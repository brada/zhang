% display the idx'th image along with detected points (red crosses) and
% projected points (green circles)
function zhang_plotpoints( calib, idx )

Calibration = zhang_load( calib );

if idx < 1 || idx > length(Calibration.Images)
    error( 'Invalid index' );
end
if ~Calibration.Images(idx).Active
    error( ['Image ',num2str(idx),' is not active'] );
end

hold off;
imshow( Calibration.Images(idx).Image );
hold on;
iPt = unhomo( Calibration.Images(idx).iPt );
pPt = unhomo( Calibration.Images(idx).pPt );
pPt = calib.N \ pPt; % bugfix reported by Ton ten Kate, 17 Feb 2011
plot( iPt(1,:), iPt(2,:), 'r+' );
plot( pPt(1,:), pPt(2,:), 'gO' );
legend( 'detected image points', 'reprojected points' );
    
    