% return the mean reprojection error [in pixels] of the calibration structure.
% If a solution has been found then the returned error is that of the solution.
% Otherwise return the error of the initial guess solution
function y = zhang_meanreprojectionerror( calib )

Calibration = zhang_load( calib );
x = Calibration.x;

active = [Calibration.Images.Active];
nPoints = size( [Calibration.Images(active).iPt], 2 );
y = sqrt( zhang_reprojectionerror( Calibration, x ) / nPoints );
