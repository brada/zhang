% returns the calibration structure unchanged if it is a calibration
% structure. But if the argument is a filename, load that file from disk
% and return the contained structure
function Calibration = zhang_load( calib )

if isstruct( calib )
    Calibration = calib;
else
    if ~exist( calib, 'file' )
        error( [calib,' does not exist'] );
    else
        Calibration = [];
        disp( 'loading calibration data into memory' );
        load( calib, 'Calibration' );
    end
end
