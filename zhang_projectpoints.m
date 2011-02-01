% given the current estimates of intrinsic and extrinsic parameters, project all
% the world points in each image onto their respective image planes
function Calibration = zhang_projectpoints( calib )

Calibration = zhang_load( calib );
A  = Calibration.A;
N  = Calibration.N;
k1 = Calibration.k1;
k2 = Calibration.k2;

for i = 1:length( Calibration.Images )
    if Calibration.Images(i).Active
        wPt = Calibration.Images(i).wPt;
        R = Calibration.Images(i).R;
        T = Calibration.Images(i).T;
        
        pPt = zhang_projectpoints_single_image( wPt, A, R, T, k1, k2 );
              
        Calibration.Images(i).pPt = pPt;
        re = zhang_reprojectionerror( Calibration, i );
        Calibration.Images(i).reprojection_error = re;
        Calibration.Images(i).mean_reprojection_error_NDC = mean( re );
        Calibration.Images(i).mean_reprojection_error_pixels = mean(re)/N(1,1);
    else
        Calibration.Images(i).pPt = [];
        Calibration.Images(i).reprojection_error = [];
        Calibration.Images(i).mean_reprojection_error_NDC = [];
        Calibration.Images(i).mean_reprojection_error_pixels = [];
    end
end
        
