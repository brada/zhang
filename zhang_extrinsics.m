% Given an existing calibration and a new image containing a grid, find the
% extrinsic parameters
% There is significant duplication of code here with detectcorners and init and
% optimise, so it could certainly be refactored for better maintainability
%
% R,T are returned rotation matrix (orthogonalised) and translation vector
% RE_* are the mean reprojection errors, in NDC and pixels
%
% UNTESTED
function [R,T,RE_NDC,RE_pixels] = zhang_extrinsics( calib, image, CALTag_datafile, debug )

if ~ismember( ndims(image), [2,3] )
    error( 'Input image must be 2 or 3 dimensional' );
end
if ~exist( CALTag_datafile, 'file' )
    error( [CALTag_datafile,' does not exist'] );
end
if nargin < 4
    debug = false;
end


Calibration = zhang_load( calib );
A  = Calibration.A;
N  = zhang_ndc_matrix( image );
k1 = Calibration.k1;
k2 = Calibration.k2;

disp( 'running CALTag' );
[wPt,iPt] = caltag( image, CALTag_datafile, false );
if ~isempty( iPt )
    % want each point as a column
    iPt = iPt';
    wPt = wPt';
    % caltag returns [row;col] and we want [x;y]
    iPt = flipud( iPt );
    wPt = flipud( wPt );
    % add homogeneous coordinates
    iPt(3,:) = 1;
    wPt(3,:) = 1;
    % normalise image coordinates to [-1,1]
    NiPt = N * iPt;
    % compute homography
    H = homography2d( wPt, NiPt );
    H = H / H(3,3);
    % due to abuse of notation we have wPt as [X;Y;1] but now we are going to
    % be dealing with 4x4 matrices, so we have to convert back to proper
    % homogeneous 3D coordinates where Z==0 and W==1
    wPt(3,:) = 0;
    wPt(4,:) = 1;
    
    % linear approximation for initial guess %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    % convert 3x3 rotation matrix to 3x1 rodrigues vector
    Rv = rodrigues_mat2vec( R );
    % reproject (assuming no lens distortion yet)
    pPt = unhomo( A*[R,T]*wPt );
    
    % nonlinear optimisation for R,T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x0 = [Rv; T];
    if debug
        plotfuncs   = {@optimplotx,@optimplotfirstorderopt,@optimplotresnorm,@optimplotstepsize};
        display     = 'iter';
        diagnostics = 'on';
    else
        diagnostics = 'off';
        plotfuncs   = [];
        display     = 'final';
    end
    JacobianPattern = ones( size(wPt,2), 6 );
    options = optimset( 'Algorithm','trust-region-reflective', ...
                        'JacobPattern',JacobianPattern, ...
                        'Diagnostics',diagnostics, ...
                        'PlotFcns', plotfuncs, ...
                        'TolX', 1e-12, ...
                        'TolFun', 1e-12, ...
                        'MaxIter', 100, ...
                        'MaxFunEvals', 1000*length(x0), ...
                        'Display',display, ...
                        'TypicalX',x0 );
    obj_fcn = @(x) single_image_reprojection_error( wPt, NiPt, A, k1, k2, x );
    disp( 'Running nonlinear optimisation' );
    [x,resnorm,residual,exitflag,output,~,jacobian] = lsqnonlin( obj_fcn, x0, [],[], options );
    if debug
        disp( 'resnorm:' );
        disp( resnorm );
        disp( 'exitflag:' );
        disp( exitflag );
        disp( 'output:' );
        disp( output );
    end
    % set output variables
    R = rodrigues_vec2mat( x(1:3) );
    T = x(4:6);
    RE_NDC = mean( residual );
    RE_pixels = mean( residual ) / N(1,1);

else
    disp( 'detection failed, cannot compute extrinsics' );
    R = [];
    T = [];
    RE_NDC = [];
    RE_pixels = [];
end


    % x is current estimate vector of extrinsics [Rx,Ry,Rz,Tx,Ty,Tz]
    function y = single_image_reprojection_error( wPt, NiPt, A, k1, k2, x )
        R = rodrigues_vec2mat( x(1:3) );
        T = x(4:6);
        
        pPt = zhang_projectpoints_single_image( wPt, A, R, T, k1, k2 );
        
        delta = pPt - NiPt;
        y = sqrt( sum(delta.^2) )';
    end

end
