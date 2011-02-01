% divide by homogeneous coordinate to make W == 1
% assumes x is a 3xN or 4xN matrix
% returns x with the W==1 coordinate included
function x = unhomo( x )

if size( x, 1 ) == 3
    x = x ./ repmat( x(3,:), [3,1] );
elseif size( x, 1 ) == 4
    x = x ./ repmat( x(4,:), [4,1] );
end
