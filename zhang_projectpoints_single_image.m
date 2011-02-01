% helper function for zhange_projectpoints and zhang_extrinsics
function pPt = zhang_projectpoints_single_image( wPt, A, R, T, k1, k2 )


% convert from world coordinates to camera coordinates
cPt = [R,T] * wPt;
% project into pinhole camera's normalised image plane at Z=1
nPt = cPt(1:2,:) ./ repmat( cPt(3,:), [2,1] );
% distance of point from optical centre
r2 = sum( nPt.^2 );
% apply radial distortion
dPt = nPt .* repmat( 1 + k1*r2 + k2*r2.^2, [2,1] );
% convert to homogeneous coords
dPt(3,:) = 1;
% apply intrinsic transformation
pPt = unhomo( A * dPt );