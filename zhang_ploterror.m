
% calib can either be a Calibration structure or a filename (incl .mat
% suffix) containing the Calibration object on disk
function zhang_ploterror( calib )

Calibration = zhang_load( calib );

hold off;
for i = 1:length( Calibration.Images )    
    if Calibration.Images(i).Active
        iPt = unhomo( Calibration.Images(i).iPt );
        %pPt = unhomo( Calibration.Images(i).pPt );
        % bugfix reported by Ton ten Kate, 17 Feb 2011
        pPt = unhomo( calib.N \ Calibration.Images(i).pPt );
        % get 2D difference vector
        err = pPt - iPt;
        plot( err(1,:), err(2,:), '+' );
        hold all;
    end
end
title( 'Reprojection error across all images' );
hold off;

        
        
        