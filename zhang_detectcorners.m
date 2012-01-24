% run caltag on all the images inside the calibration structure, saving the
% detected image points and corresponding world points in each image structure.
% You may pass in a Calibration structure or else a filename of a .mat file
% containing that structure as the calib argument. CALTag_datafile is the .mat
% file containing the description on the caltag pattern being used.
% If reading the Calibration structure from file it will be saved to file upon
% exit, in addition to being returned as an argument.
function Calibration = zhang_detectcorners( calib, CALTag_datafile )

if ~exist( CALTag_datafile, 'file' )
    error( [CALTag_datafile,' does not exist'] );
end
Calibration = zhang_load( calib );

for i = 1:length( Calibration.Images )
    
    if Calibration.Images(i).Active
        I = Calibration.Images(i).Image;
        disp( ['running CALTag on ',Calibration.Images(i).Name] );
        [wPt,iPt] = caltag( I, CALTag_datafile, false );
        if size( iPt, 1 ) > 4
            % want each point as a column
            iPt = iPt';
            wPt = wPt';
            % caltag returns [row;col] and we want [x;y]
            iPt = flipud( iPt );
            wPt = flipud( wPt );
            % add homogeneous coordinates
            iPt(3,:) = 1;
            wPt(3,:) = 1;
            % normalise image coordinates to NDC, i.e. [-1,1] horizontally
            NiPt = Calibration.N * iPt;
            % compute homography per image
            H = homography2d( wPt, NiPt );
            H = H / H(3,3);
            % due to abuse of notation we have wPt as [X;Y;1] but now we are going to
            % be dealing with 4x4 matrices, so we have to convert back to proper
            % homogeneous 3D coordinates where Z==0 and W==1
            wPt(3,:) = 0;
            wPt(4,:) = 1;
            % save to structure
            nPoints = size( iPt, 2 );
            disp( [num2str(nPoints),' points detected'] );
            Calibration.Images(i).iPt  = iPt;
            Calibration.Images(i).NiPt = NiPt;
            Calibration.Images(i).wPt  = wPt;
            Calibration.Images(i).H    = H;  
            Calibration.Images(i).numPoints = size( iPt, 2 );
        else
            disp( 'detection failed, disabling image' );
            Calibration.Images(i).Active = false;
            Calibration.Images(i).iPt  = [];
            Calibration.Images(i).NiPt = [];
            Calibration.Images(i).wPt  = [];
            Calibration.Images(i).H    = [];
            Calibration.Images(i).numPoints = 0;
        end
    else
        Calibration.Images(i).iPt  = [];
        Calibration.Images(i).NiPt = [];
        Calibration.Images(i).wPt  = [];
        Calibration.Images(i).H    = [];
        Calibration.Images(i).numPoints = 0;
    end
     
end


% output
if ischar( calib )
    disp( 'saving to disk' );
    save( calib, 'Calibration', '-append' );
    disp( 'done' );
end

