% if x is a parameter vector then this
% returns a vector of the reprojection error per point [in pixels]. Note
% that it does NOT return the square of the error, or the sum of the
% squares. Other function (like fminunc) expect the sum of squares, but
% lsqnonlin (which we use) expects a vector of unsquared function values.
%
% if x is a scalar, it should refer to the index of an image to compute the
% reprojection error vector for all points in that image.
function y = zhang_reprojectionerror( calib, x )


Calibration = zhang_load( calib );


if isscalar( x )
    
    % want to return the reprojection error vector just for the xth image
    y = [];
    if x > 0 && x <= length( Calibration.Images )
        if Calibration.Images(x).Active
            pPt = Calibration.Images(x).pPt;
            NiPt = Calibration.Images(x).NiPt;
            if isempty( NiPt )
                warning( 'Zhang:nopoints', 'No points in that image' );
            else
                delta = pPt - NiPt;
                y = sqrt( sum(delta.^2) )';
                % if you want the error in pixels, then scale it by W/2 i.e.
                % inverse of the top left entry of the NDC matrix
            end
        end
    else
        warning( 'Zhang:arg', 'Trying to access nonexistant image' );
    end
    
else

    % want to return the reprojection error vector for all active images

    Calibration = zhang_unpackparamvector( Calibration, x );
    Calibration = zhang_projectpoints( Calibration );

    active = [Calibration.Images.Active];
    Images = Calibration.Images;
    y = cell( length(Images), 1 );
    for i = 1:length( Images )
        y{i} = zhang_reprojectionerror( Calibration, i );
    end
    y = cell2mat( y(active) );

end

