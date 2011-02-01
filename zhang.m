
% Zhang's calibration algorithm
% as per "Zhang's camera calibration step by step" by Teixeira et al
function A = zhang()


% %% generate point correspondences
% % load caltag
% addpath( '~/imager/CALTag-1.0' );
% datafile = '~/imager/CALTag-1.0/CALTagData.mat';
% % divide by homogeneous coordinate to make w==1
% unhomo = @(x) x ./ repmat( x(3,:), [3,1] );
% % loop through JPG files in directory
% files = dir( '*.JPG' );
% for i = 1:length( files )
%     f = files(i).name;
%     I = imread( f );
%     % build normalisation matrix
%     [height,width,depth] = size( I );
%     N = [ 2/width     0      -1
%              0     2/height  -1
%              0        0       1 ];
%     disp( ['Running caltag on ', f] );
%     [wPt,iPt] = caltag( I, datafile, false );
%     if ~isempty( iPt )
%         disp( 'Detection successful' );
%         % want each point as a column
%         wPt = wPt';
%         iPt = iPt';
%         % caltag returns [row;col] and we want [x;y]
%         iPt = flipud( iPt );
%         wPt = flipud( wPt );
%         % add homogeneous coordinates
%         iPt(3,:) = 1;
%         wPt(3,:) = 1;
%         % normalise image coordinates to [-1,1]
%         NiPt = N * iPt;
%         % compute homography per image
%         H = homography2d( wPt, NiPt );
%         H = H / H(3,3);
%         % due to abuse of notation we have wPt as [X;Y;1] but now we are going to
%         % be dealing with 4x4 matrices, so we have to convert back to proper
%         % homogeneous 3D coordinates where Z==0 and W==1
%         wPt(3,:) = 0;
%         wPt(4,:) = 1;
%         % save results. because we store a separate N matrix for each image
%         % we could potentially input images of different resolutions
%         [pathstr,name] = fileparts( f );
%         outputfile = [name,'.mat'];
%         save( outputfile, 'name', 'wPt', 'iPt', 'NiPt', 'N', 'H' );
%     end
% end
% clear i;


%% compute intrinsic calibration guess

% divide by homogeneous coordinate to make w==1
unhomo = @(x) x ./ repmat( x(3,:), [3,1] );

disp( 'Initialising intrinsic parameters' );
vij = @(i,j,H) [ H(1,i)*H(1,j)
                 H(1,i)*H(2,j) + H(2,i)*H(1,j)
                 H(2,i)*H(2,j)
                 H(3,i)*H(1,j) + H(1,i)*H(3,j)
                 H(3,i)*H(2,j) + H(2,i)*H(3,j)
                 H(3,i)*H(3,j) ];
files = dir( '*.mat' );
L = length( files );
V = [];
for j = 1:length( files )
    f = files(j).name;
    load( f, 'wPt', 'iPt', 'NiPt', 'N', 'H' );
    % append G matrix to V
    G = [ vij(1,2,H)'; (vij(1,1,H)-vij(2,2,H))' ];
    V = [ V; G ];
end

% solving Vb=0 directly gives b==0 so instead we find the eigenvector
% associated with the smallest eigenvalue of V'V. Actually we use svd
% instead because it automatically sorts singular vals in descending order
[LEFT,SIGMA,RIGHT] = svd( V );
b = RIGHT(:,end);

% eigenvalue approach, gives same results as above svd method
% [eigenvecs,eigenvals] = eig( V'*V );
% [mineigen,idx] = min( diag(eigenvals) );
% b = eigenvecs(:,idx);

% build A
%v0 = ( b(2)*b(4)-b(1)*b(5) ) / ( b(4)*b(3)-b(2)^2 ); % from Teixeria (with error)
v0 = ( b(2)*b(4)-b(1)*b(5) ) / ( b(1)*b(3)-b(2)^2 ); % from Mert: correct
%lambda = b(6) - ( b(2)^2 + v0*(b(2)*b(4)-b(1)*b(5)) ) / b(1); % from Teixeira (with error)
lambda = b(6) - ( b(4)^2 + v0*(b(2)*b(4)-b(1)*b(5)) ) / b(1); % from Mert: correct
alpha = sqrt( lambda / b(1) );
beta = sqrt( lambda*b(1) / (b(1)*b(3)-b(2)^2) );
gamma = -b(2)*alpha^2*beta / lambda;
u0 = gamma*v0 / beta - b(4)*alpha^2 / lambda;

% form initial guess for A
A = [ alpha  gamma  u0
      0      beta   v0
      0      0      1   ];
% calling N\A is equiv to inv(N)*A
A = N \ A;
% optional
disp( ['Forcing skew to zero (from ',num2str(A(1,2)),')'] );
A(1,2) = 0;


%% get initial guesses for R & T matrices for each image
for j = 1:length( files )
    f = files(j).name;
    disp( ['Calculating RT for ',f] );
    load( f, 'H', 'N' );
    % undo the normalisation
    H = N\H;
    % should already have homogeneous_w == 1, but just in case do it again
    H = H / H(3,3);
    h1 = H(:,1);
    h2 = H(:,2);
    h3 = H(:,3);
    % note this deviates from Zhang's paper a bit because ||h1|| != ||h2||
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
    % debug: plotting projected points like this works
    % pPt = unhomo(A*([R,T]*wPt)); plot( pPt(1,:),pPt(2,:),'r+' );
    % write to disk
    save( f, '-append', 'A', 'R', 'T' );
end


% Gordon's experiments with the translation stage show that you can get
% much better results by enforcing constraints on the position of the camera
% or plane, and also by running the optimiser for longer. So the questions is:
% did he get better results just via running longer or did the extra
% constraints really help. If they did then it must be because you're adding
% known 3D distances in the Z direction, whereas with a plane you only know
% the true distances in XY. So if you have a translation stage you can get
% the 3rd dimension, but if you don't then you have to make sure to rotate
% the plane to a very steep angle. But then your saddle point locator becomes
% less accurate. So can we design a better saddle finder that takes the local
% grid distortion into account? And also, can we enforce the constraint that
% the intrinsic matrix should be invariant to the image set used to calibrate?
% ie instead of minimising just reprojection error across all images, we do
% that in an inner loop for a set of N-1 images and then in the outer loop
% minimise the distance between the intrinsic matrices for each of the N image
% subsets. Is that somehow equivalent to having a known 3rd dimension
% measurement in that we resolve the ambiguity between Z distance and zoom?


  
%% levenberg-marquardt to minimise reprojection error

% parameter vector is (ignoring skew): 
%    [ (focalx; focaly);  
%      (principlex; principley);
%      (R1x;R1y;R1z); (T1x;T1y;T1z);
%      (R2x;R2y;R2z); (T2x;T2y;T2z);
%      ...; 
%      (Rnx;Rny;Rnz); (Tnx;Tny;Tnz) ];
% The rotation parameters are rodrigues vectors

x0 = [ A(1,1); A(2,2); A(1,3); A(2,3) ];
mappings = {};
files = dir( '*.mat' );
for i = 1:length( files )
    f = files(i).name;
    load( f, 'R', 'T', 'wPt', 'iPt' );
    mappings{i} = struct( 'wPt',wPt, 'iPt',iPt, 'R',R, 'T',T );
    Rv = rodrigues_mat2vec( R );
    x0 = [ x0; Rv; T ];
end
    
options = optimset( 'Algorithm',{'levenberg-marquardt',0.01}, ...
                    'Diagnostics','on', ...
                    'PlotFcns', {@optimplotfunccount,@optimplotresnorm}, ...
                    'TolX', 1e-12, ...
                    'MaxFunEvals', 10000*length(x0) );
%x = lsqnonlin( @reprojection_error, x0, [],[], options )


% TODO: create function to plot all images with the reprojected points from
% a given x vector
% TODO: alternate between this optimisation and adjusting the radial
% distortion
% TODO: there is still something wrong though, mean reprojection error is
% over 100 pixels currently I think


    %% NB: because this is an inner function it has access to the variables
    % in the caller. This can lead to some very hard-to-trace bugs if you
    % reuse a variable name inside here. Hence the clumsy obj_ namespace
    % hack.
    function y = reprojection_error( x )
        
        y = 0;
        
        % unpack x vector
        obj_A = [ x(1)  0     x(3)
                  0     x(2)  x(4)
                  0     0     1    ];
        Nplanes = (length(x)-4) / 6;
        for obj_i = 1:Nplanes
            first = (obj_i-1)*6 + 5;
            obj_R = rodrigues_vec2mat( x(first:first+2) );
            obj_T = x(first+3:first+5);
            obj_wPt = mappings{obj_i}.wPt;
            obj_iPt = mappings{obj_i}.iPt;
            obj_pPt = unhomo( obj_A*[obj_R,obj_T]*obj_wPt );
            obj_err = sum( (obj_iPt - obj_pPt).^2 );
            y = y + sum( obj_err );
        end
        
       % y = sum( (iPt - unhomo(N\(A*[R,T]*wPt))).^2 );
    end


end



%% plotting reprojection error

i = 1;

[filepath,filebase,fileext] = fileparts( files(i).name );
filename = [filebase,'.JPG'];
I = imread( filename );
imshow( I );
hold on;

% plotting
Qi = unhomo( mappings{i}.iPt );
Qw = mappings{i}.wPt;
plot( Qi(1,:), Qi(2,:), 'r+' );
mappedQw = unhomo( A * [mappings{i}.R mappings{i}.T] * Qw );
plot( mappedQw(1,:), mappedQw(2,:), 'gO' );
hold off;







