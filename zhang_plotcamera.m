function zhang_plotcamera( calib, idx )

Calibration = zhang_load( calib );

if Calibration.Images(idx).Active
    
    % plot single image plane for static grid
    wPt = Calibration.Images(idx).wPt;
    minW = min( wPt, [], 2 );
    maxW = max( wPt, [], 2 );
    A = minW;
    B = [ minW(1); maxW(2); 0; 1 ];
    C = maxW;
    D = [ maxW(1); minW(2); 0; 1 ];
    cnr = [A,B,C,D];
    colour = [1;1;1;1]; 
    patch( cnr(1,:), cnr(2,:), cnr(3,:), colour );
    gridSize = max( maxW - minW ) / 2;
      
    % inverse transformation
    R = Calibration.Images(idx).R;
    T = Calibration.Images(idx).T;
    M = [R, T; [0,0,0,1]];
   
    % camera coordinates of the arrows we want to draw
    cop = [0;0;0;1]; % centre of projection
    oa  = [0;0;norm(T);1]; % optical axis (use norm(T) for length)
    up  = [0;norm(T)/10;0;1]; % unit up vector
    rt  = [norm(T)/10;0;0;1]; % right vector
    
    % a sphere at the camera's position
    [Sx,Sy,Sz] = sphere( 10 ); 
    [nrows,ncols] = size(Sx);
    S = [ [Sx(:), Sy(:), Sz(:)]*gridSize/2, ones(numel(Sx),1) ];
    S = (M \ S')';
    Sx = reshape( S(:,1), [nrows,ncols] );
    Sy = reshape( S(:,2), [nrows,ncols] );
    Sz = reshape( S(:,3), [nrows,ncols] );
      
    % lines representing the positive camera axes (longest one is the
    % optical axis)
    cop = M\cop;
    oa = M\oa;
    up = M\up;
    rt = M\rt;
      
    % draw the sphere and lines
    surf( Sx, Sy, Sz );
    daspect([1,1,1]);
    hold on
    line_oa = [ cop, oa ];
    line( line_oa(1,:), line_oa(2,:), line_oa(3,:) );
    line_up = [ cop, up ];
    line( line_up(1,:), line_up(2,:), line_up(3,:) );
    line_rt = [ cop, rt ];
    line( line_rt(1,:), line_rt(2,:), line_rt(3,:) );     
    
end