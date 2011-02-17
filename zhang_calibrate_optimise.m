% calls optimisation toolbox to improve upon the initial guess x0 of the
% calibration parameters, by minimising the reprojection error
%
% debug [bool]: True to output residual at each iteration. 
%
% start_from_cur_soln [bool]: True to begin optimising from the most
% recent x stored in the calibration object (otherwise it will begin at x0,
% the initial guess)
%
% maxiter [int]: maximum number of iterations to perform
%
function [Calibration,x] = zhang_calibrate_optimise( calib, debug, start_from_cur_soln, maxiter )

if nargin < 4
    maxiter = 100;
end
if nargin < 3
    start_from_cur_soln = false;
end
if nargin < 2
    debug = false;
end

Calibration = zhang_load( calib );

if start_from_cur_soln
    x0 = Calibration.x;
else
    x0 = Calibration.x0;
end

if debug
    plotfuncs   = {@optimplotx,@optimplotfirstorderopt,@optimplotresnorm,@optimplotstepsize};
    display     = 'iter';
    diagnostics = 'on';
else
    diagnostics = 'off';
    plotfuncs   = [];
    display     = 'off';
end

% form the sparsity pattern for the (nPoints x nParams) jacobian matrix. 
% because the intrinsic parameters affect all points we have ones all down
% the left hand side (for each column corresponding to an intrinsic
% variable). then because the extrinsic parameters for each image only
% affect that image, the columns of J corresponding to those extrinsic
% parameters have a block diagonal form (ie changes in the, say, x-rotation
% of the 2nd plane will not affect the reprojection error of points in the
% 3rd image. Hence J(point_i in 3rd image,col 1 of 2nd plan) = 0. Supplying
% this sparsity information allows lsqnonlin to avoid computing those zero
% entries via finite differences (only for the trust-region-reflective
% algorithm though) which gives huge speed boost since it avoids many
% function evaluations.
active = [Calibration.Images.Active];
Images = Calibration.Images(active);
blocks = cell( length(Images), 1 );
nIntrinsics = 7;
if Calibration.square_aspect
    nIntrinsics = nIntrinsics - 1;
end
if Calibration.zero_skew
    nIntrinsics = nIntrinsics - 1;
end
if Calibration.first_order
    nIntrinsics = nIntrinsics - 1;
end
nExtrinsics = 6;
for i = 1:length( Images )
    blocks{i} = sparse( ones( Images(i).numPoints, nExtrinsics ) );
end
J = blkdiag( blocks{:} );
J = [ ones(size(J,1),nIntrinsics), J ];
Calibration.JacobianPattern = J;


options = optimset( 'Algorithm','trust-region-reflective', ...
                    'JacobPattern',Calibration.JacobianPattern, ...
                    'Diagnostics',diagnostics, ...
                    'PlotFcns', plotfuncs, ...
                    'TolX', 1e-12, ...
                    'TolFun', 1e-12, ...
                    'MaxIter', maxiter, ...
                    'MaxFunEvals', 1000*length(x0), ...
                    'Display',display, ...
                    'TypicalX',x0 );
                              
obj_fcn = @(x) zhang_reprojectionerror(Calibration,x);

disp( 'Running nonlinear optimisation' );
[x,resnorm,residual,exitflag,output,~,jacobian] = lsqnonlin( obj_fcn, x0, [],[], options );

Calibration.x = x;
Calibration.Jacobian = jacobian;
Calibration.output = output;
Calibration.exitflag = exitflag;
Calibration.residual = residual;

% get the 95% confidence interval
Calibration.confidence_interval = nlparci( x, residual, 'Jacobian',jacobian );

% apply the current solution to all the image extrinsic parameters and
% project the points according to this solution
Calibration = zhang_unpackparamvector( Calibration, x );
Calibration = zhang_projectpoints( Calibration );

% get the mean reprojection error across all images
mean_re_ndc = [Calibration.Images(:).mean_reprojection_error_NDC];
Calibration.mean_reprojection_error_NDC = mean_re_ndc;
Calibration.mean_reprojection_error_pixels = mean_re_ndc / Calibration.N(1,1);


if debug
    disp( 'resnorm:' );
    disp( resnorm );
    disp( 'exitflag:' );
    disp( exitflag );
    disp( 'output:' );
    disp( output );
end


if ischar( calib )
    disp( 'saving to disk' );
    save( calib, 'Calibration', '-append' );
    disp( 'done' );
end

