function R = rodrigues_vec2mat( v )
% copied and simplified from Bouguet toolbox

if length( v ) ~= 3
    error( 'v must have 3 elements' );
end
theta = norm( v );
if theta < eps
    R = eye( 3 );
else
    v = v(:);
    omega = v / theta;
    alpha = cos( theta );
    beta = sin( theta );
    gamma = 1 - cos( theta );
    omegav = [ 0          -omega(3)  omega(2)
               omega(3)   0          -omega(1)
               -omega(2)  omega(1)   0         ];
    A = omega * omega';
    R = eye(3)*alpha + omegav*beta + A*gamma;
end