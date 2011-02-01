function R = reorthogonalise( R )
% if R isn't a valid 3x3 rotation matrix, return the closest valid rotation
% matrix to R, in terms of the Frobenius norm

if ~isequal( size(R), [3,3] )
    error( 'R must be 3x3' );
end

[U,S,V] = svd( R );

% the SVD can potentially cause a reflection so we check to see if this
% happened so we can correct it later
S = diag( [1 1 det(U*V')] );

% Zhang just says R=U*V' but other references say the S accommodates for
% reflections
R = U*S*V';