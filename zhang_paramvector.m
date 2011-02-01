% take all the essential parameters for optimisation out of the Calibration
% structure and place them into a single 1D vector for use in the optimisation
% routine. The order of the vector is:
% [ fx; fy;   // focal length divided by sensor size [or just f if squareaspect] 
%   skew;     // skew of sensor axes [optional]
%   px; py;   // principle point
%   k1;k2;    // radial distortion coefficients [or just k if firstorder]
%   Rx;Ry;Rz; // rodrigues rotation vector of first plane
%   Tx;Ty;Tz; // translation vector of first plane
%   ...
%   Rx;Ry;Rz; // rodrigues rotation vector of last plane
%   Tx;Ty;Tz; // translation vector of last plane
% ]
function x = zhang_paramvector( calib )

Calibration = zhang_load( calib );

% extract parameters
A  = Calibration.A;
k1 = Calibration.k1;
k2 = Calibration.k2;
active = [Calibration.Images.Active];
Rv = [Calibration.Images(active).Rv];
T  = [Calibration.Images(active).T];
RvT = [ Rv; T ];

% compose parameters into 1D vector
if Calibration.square_aspect
    focal = A(1,1);
else
    focal = [ A(1,1); A(2,2) ];
end
if Calibration.zero_skew
    skew = [];
else
    skew = A(1,2);
end
if Calibration.first_order
    k = k1;
else
    k = [k1; k2];
end
principle = A([1,2],3);
x = [ focal; skew; principle; k; RvT(:) ];
