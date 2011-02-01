% zero_skew: bool. True if you want to force the skew to zero. Note that this is
% not necessarily a good idea. Even if the sensor is perfectly orthogonal, you
% can still have nonzero skew if the sensor is not perfectly orthogonal to the
% optical axis, which could happen with a slightly misaligned lens.
% if centre_pp is true, the principle point will be initialised to the centre of
% the image. if square_aspect is on, then fx and fy in A will be forced to be
% equal.
% Note that centre_pp does not enforce that constraint during the optimisation -
% it only affects the initial guess. If false, then an estimate will be computed
% via linear algebra.
% first_order: bool. True if you want to use a first order radial distortion
% model (ie only using the r^2 term). If false, a second order model will be
% used (ie r^2 and r^4). Second order models more complicated lens distortion,
% however the first order term dominates in most decent quality lenses, and so
% including the second order term only leads to numerical instability. You can
% see this by checking the confidence_interval after calibrating - the interval
% corresponding to the second order interval could be much larger (as a
% percentage of the magnitude of the actual value for that parameter) than the
% others. If you set this to false, the k2 (r^4) term will be fixed to zero.
%
% The n* and delta* parameters are particular to the glass scanning
% calibration protocol. They describe the number of 'free' calib images,
% the number of images captured using a translation stage moving in deltaT
% mm increments, and the number of images captured using the rotation stage
% moving in deltaR degree increments. The images must be captured in that
% order. deltaT and deltaR should always be positive, regardless of the
% direction in which the planes are moved, because it refers to the
% absolute delta between consecutive planes.
function Calibration = zhang_calibrate_init( calib, zero_skew, ...
                                             centre_pp, square_aspect, ...
                                             first_order, nFree, nTStage, ...
                                             nRStage, deltaT, deltaR )

Calibration = zhang_load( calib );
Calibration.centre_pp = centre_pp;
Calibration.zero_skew = zero_skew;
Calibration.square_aspect = square_aspect;
Calibration.first_order = first_order;
Calibration.nFree = nFree;
Calibration.nTStage = nTStage;
Calibration.nRStage = nRStage;
Calibration.deltaT = deltaT;
Calibration.deltaR = deltaR;


% helper function: build V matrix from H
vij = @(i,j,H) [ H(1,i)*H(1,j)
                 H(1,i)*H(2,j) + H(2,i)*H(1,j)
                 H(2,i)*H(2,j)
                 H(3,i)*H(1,j) + H(1,i)*H(3,j)
                 H(3,i)*H(2,j) + H(2,i)*H(3,j)
                 H(3,i)*H(3,j) ];


V = [];
for i = 1:length( Calibration.Images )
    if Calibration.Images(i).Active
        % append G matrix to V
        H = Calibration.Images(i).H;
        G = [ vij(1,2,H)'; (vij(1,1,H)-vij(2,2,H))' ];
        V = [ V; G ];
    end
end


% solving Vb=0 directly gives trivial b==0 so instead we find the eigenvector
% associated with the smallest eigenvalue of V'V. Actually we use svd
% instead because it automatically sorts singular vals in descending order
[~,~,RIGHT] = svd( V );
b = RIGHT(:,end);

% eigenvalue approach, gives same results as above svd method
% [eigenvecs,eigenvals] = eig( V'*V );
% [mineigen,idx] = min( diag(eigenvals) );
% b = eigenvecs(:,idx);


% build A
% note that Teixeria's document had an error in the v0 and lambda formulae
v0 = ( b(2)*b(4)-b(1)*b(5) ) / ( b(1)*b(3)-b(2)^2 );
lambda = b(6) - ( b(4)^2 + v0*(b(2)*b(4)-b(1)*b(5)) ) / b(1);
% [10 Dec 2010] this warning placed here in response to bug report from Carsten
% Fries. Basically, if all the images are from roughly the same viewpoint, the
% lambda will be negative, and then you get imaginary numbers in the A matrix
% due to sqrt of negative number. Need to capture better images. In the
% meanwhile, we simply clamp lambda and beta to near 0 (since it's just an
% initial guess) and then hope that the optimisation will still work. It
% actually did work on Carsten's data, but for lambda=0.001 it returned -ve
% focal distances. This suggests that we should actually be using a bounded
% optimisation in zhang_calibrate_optimise
if lambda <=0
    disp( ['Warning: negative lambda. Calibration may be unreliable, ' ...
           'please capture more images at more varied angles.'] );
    lambda = 0.1;
end
alpha = sqrt( lambda / b(1) );
beta = sqrt( lambda*b(1) / (b(1)*b(3)-b(2)^2) );
if imag(beta) ~= 0
    beta = 0.1;
end
gamma = -b(2)*alpha^2*beta / lambda;
u0 = gamma*v0 / beta - b(4)*alpha^2 / lambda;


% form initial guess for A
if centre_pp
    u0 = 0;
    v0 = 0;
end
if zero_skew
    gamma = 0;
end
A = [ alpha  gamma  u0
        0    beta   v0
        0      0    1   ];
%N = Calibration.N;
% want NDC not pixel coordinates, so disable the normalisation undo here
%A = N\A;
if square_aspect
    A(1,1) = ( A(1,1)+A(2,2) ) / 2;
    A(2,2) = A(1,1);
end



% initial guess for radial distortion
Calibration.k1 = 0;
Calibration.k2 = 0;

% get initial guesses for R & T matrices for each image
for i = 1:length( Calibration.Images )
    if Calibration.Images(i).Active
        name = Calibration.Images(i).Name;
        H    = Calibration.Images(i).H;
        wPt  = Calibration.Images(i).wPt;
        disp( ['calculating A,R,T for ',name] );
        % should already have homogeneous_w == 1, but just in case do it again
        H = H / H(3,3);
        h1 = H(:,1);
        h2 = H(:,2);
        h3 = H(:,3);
        % note this deviates from Zhang's paper a bit because ||h1|| != ||h2||
        % note also that A\h1 is equivalent (but superior) to inv(A)*h1
        lambda1 = 1 / norm( A\h1 );
        lambda2 = 1 / norm( A\h2 );
        lambdat = (lambda1 + lambda2) / 2;
        r1 = lambda1 * (A\h1);
        r2 = lambda2 * (A\h2);
        r3 = cross( r1, r2 );
        t = lambdat * (A\h3);
        R = [ r1 r2 r3 ];
        T = [ t ];
        % now R isn't a valid rotation matrix, so reorthogonalise it
        % (this is not in Mert's code)
        R = reorthogonalise( R );
        % convert 3x3 rotation matrix to 3x1 rodrigues vector
        Rv = rodrigues_mat2vec( R );
        % reproject (assuming no lens distortion, and in NDC)
        pPt = unhomo( A*[R,T]*wPt );
        % output
        Calibration.Images(i).T   = T;
        Calibration.Images(i).R   = R;
        Calibration.Images(i).Rv  = Rv;
        Calibration.Images(i).pPt = pPt;
    end
end

% output
Calibration.A0 = A;
Calibration.A  = A;
Calibration.x0 = zhang_paramvector( Calibration );
Calibration.x  = Calibration.x0;
disp( 'Initial guess for A:' );
disp( A );
if ischar( calib )
    disp( 'saving to disk' );
    save( calib, 'Calibration', '-append' );
    disp( 'done' );
end


