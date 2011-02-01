function zhang_plotplane( calib, idx )

Calibration = zhang_load( calib );

if Calibration.Images(idx).Active
    % get the corners of the detected points
    % assume they all lie in a constant Z plane (usually Z=0)
    wPt = Calibration.Images.wPt;  
    minW = min( wPt, [], 2 );
    maxW = max( wPt, [], 2 );
    A = minW;
    B = [ minW(1); maxW(2); 0; 1 ];
    C = maxW;
    D = [ maxW(1); minW(2); 0; 1 ];
    cnr = [A,B,C,D];
    
    R = Calibration.Images(idx).R;
    T = Calibration.Images(idx).T;
    cnr = [R,T] * cnr;
    colour = [1;1;1;1];
    
    patch( cnr(1,:), cnr(2,:), cnr(3,:), colour, 'FaceAlpha',0.3 );
end
