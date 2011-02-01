% return a vector of values that should be zero if the proposed parameter
% solution x satisfies the constraints. this corresponds to the ceq
% function in fmincon. so if you want to enforce a constraint that the
% first two parameters should sum to 8 and the 3rd and 4th parameters
% should sum to 9, then you would return y = [ x(1)+x(2)-8; x(3)+x(4)-9 ];
% note that this function is entirely dependent on the particular
% calibration setup: if you know that the second last and last photos are
% exactly 10mm apart (because you used a translation stage to capture them)
% then you could enforce that constraint on the T parameter.
%
% For the GlassScan project, we take Nfree images of the LCD monitor being
% manually moved around, then Ntstage images of it being translated from its
% current position, away from the camera 10mm each step, and then Nrstage images
% of the grid on the rotation stage. deltaT is the distance between successive
% nTStage grids (in mm) and deltaR is the rotation angle between successive
% nRStage grids (in degrees)
function [c,ceq] = zhang_constraint_eq( calib, x )

Calibration = zhang_load( calib );
nFree = Calibration.nFree;
nRStage = Calibration.nRStage;
nTStage = Calibration.nTStage;
deltaT = Calibration.deltaT;
deltaR = Calibration.deltaR;

c = [];

% EXAMPLE:
%
% This function currently assumes that the final 3 images are taken with a
% translation stage exactly 1in apart. Hence the R parameters for all three
% images should be equal. Als, image A and B are 1in apart,
% and B and C are 1in apart, and A and C are 2in apart. This final
% transitive constraint implicitly constrains the 3 grids to undergo pure
% translation along a 1D line.
%
% Calibration = zhang_unpackparamvector( Calibration, x );
% active = [Calibration.Images.Active];
% Images = Calibration.Images(active);
% T3 = Images(end).T;
% R3 = Images(end).Rv;
% T2 = Images(end-1).T;
% R2 = Images(end-1).Rv;
% T1 = Images(end-2).T;
% R1 = Images(end-2).Rv;
% 
% enforcedDistance = 25.4;
% ceq = [R1-R2; R2-R3; norm(T1-T2)-enforcedDistance];
%
% END OF EXAMPLE


Calibration = zhang_unpackparamvector( Calibration, x );
Images = Calibration.Images;
ImagesTStage = Images(nFree+1:nFree+nTStage);
ImagesRStage = Images(nFree+nTStage+1:nFree+nTStage+nRStage);


% extract T and R vectors for the translation stage images
T_TStage = [ImagesTStage(:).T];
R_TStage = [ImagesTStage(:).Rv];
% find distance between pairs (final element is distance first to last)
% setting distance first to last to be (n-1)*delta ensures that the points
% are all collinear
D = sqrt( sum( (T_TStage - circshift(T_TStage,[0,-1])).^2 ) );
% deltaT mm separation, in inches
delta = deltaT / 25.4;
cnstrnt_TStage_T  = [ abs(D(1:end-1)-delta), abs(D(end)-delta*(nTStage-1)) ];
% same rotation vector for each image
cnstrnt_TStage_Rv = sqrt( sum( (R_TStage - circshift(R_TStage,[0,-1])).^2 ) );

% plane normals (Z-axes) for RStage grids
% alas, "Z=[ImagesRStage(:).R(:,3)]" is invalid matlab syntax
Z = [ImagesRStage(:).R];
Z = Z(:,3:3:end);
% rotation angle between successive grids
D = acos( dot( Z, circshift(Z,[0,-1]) ) );
% deltaR deg separation, in radians
delta = deg2rad(deltaR);
cnstrnt_RStage_R = [ abs(D(1:end-1)-delta), abs(D(end)-delta*(nRStage-1)) ];


% assemble all constraints
ceq = [cnstrnt_TStage_T, cnstrnt_TStage_Rv, cnstrnt_RStage_R];

