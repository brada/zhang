% init calibration structure. This will copy the images into the .mat file
% so be aware that it grow to a large size if you have many images
% The output_filename is optional
function Calibration = zhang_init( image_pattern, output_filename )


if ~ischar( image_pattern )
    error( 'image_pattern should be a filename pattern, eg *.JPG' );
end
if nargin > 1 && exist( output_filename, 'file' )
    warning( 'Zhang:arg', [output_filename, ' will be overwritten'] );
end


files = dir( image_pattern );
image_path = fileparts( image_pattern );
if isempty( files )
    error( 'No files found matching pattern' );
end

images = [];
N_0    = [];
W_0    = [];
H_0    = [];
for i = 1:length( files )

    name = files(i).name;
    if ~isempty( image_path )
        name = [image_path, '/', name];
    end
    disp( ['loading ',name] );
    I = imread( name );
    if ndims( I ) == 3
        I = rgb2gray( I );
    end
    I = im2uint8( I );

    [height,width] = size( I );
    N = zhang_ndc_matrix( I );
    % all images will be resized to match the first one
    if isempty( N_0 )
        N_0 = N;
        W_0 = width;
        H_0 = height;
    else
        if ~isequal( N, N_0 )
            % note this does not accommodate for changes in aspect ratio of the
            % other images, only uniform scaling
            warning( 'Zhang:res', [name,' has incompatible resolution: resizing'] );
            I = imresize( I, [H_0,W_0] );
        end
    end

    img = struct( 'Name',name, 'Image',I, 'Active',true );
    images = [images;img];
end


% output
Calibration = struct( 'A',[], 'N',N_0, 'Images',images );
if nargin > 1
    disp( 'saving to disk' ); 
    if exist( output_filename, 'file' )
        save( output_filename, 'Calibration', '-append' );
    else
        save( output_filename, 'Calibration' );
    end
    disp( 'done' );
end

