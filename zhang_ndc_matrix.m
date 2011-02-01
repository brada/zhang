% compute the matrix to map pixel coordinates in an image to "normalised device
% coordinates" which lie in [-1,1] horizontally and [-H/W,H/W] vertically where
% H and W are height and width of input image
% image must be a valid image matrix
function N = zhang_ndc_matrix( image )

if ~ismember( ndims(image), [2,3] )
    error( 'Input image must be 2 or 3 dimensional' );
end


dims = size( image );
height = dims(1);
width = dims(2);

% old normalisation matrix
%N = [ 2/width     0      -1
%         0     2/height  -1
%         0        0       1 ];
% new matrix, after email discussion started 21 August 2010 with Wolfgang
% and Derek. The problem with using different W and H scale factors is that
% you transform into a unit square so you lose the aspect ratio, so that
% radial undistortion would apply elliptical correction rather than circular
% This matrix also accounts for the (aspect corrected) half pixel shift so
% that -1/2 pixels maps to -1 "normalised pixels" and W-1+1/2 maps to +1
% because the pixel centres lie at integer coordinates. See my notebook,
% page 1st September 2010
N = [ 2/width     0      -1+1/width
         0     2/width   -height/width+1/width
         0        0       1                     ];


% NOTE however that MDA uses "y-flip" which means that the Y axis points up,
% rather than down, as in most image formats. To accommodate that, we would have
% to negate the signs of the 2nd row of the above matrix. Doing so will
% also ensure we get a nice, consistent right hand coordinate system
% (otherwise we end up with weird y-flips in the 3D world)
% UPDATE: this does NOT produce the desired effect. Even with/without
% setting all the image point y-coords to "imgheight-y" does not work.
%N = [ 2/width     0         -1+1/width
%         0     -(2/width)   -(-height/width+1/width)
%         0        0          1                       ];
