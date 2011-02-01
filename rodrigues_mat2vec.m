function v = rodrigues_mat2vec( R )
% copied and simplified from Bouguet toolbox

[M,N] = size( R );
if M ~= N && N ~= 3
    error( 'R must be a 3x3 matrix' );
end


[U,S,V] = svd(R);
R = U*V';

tr = (trace(R)-1)/2;
theta = real( acos(tr) );

if sin(theta) >= 1e-5,
    
    vth = 1/(2*sin(theta)); 
    om1 = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';
    om = vth*om1;
    v = om*theta;
   
else
    % case norm(om)=0; 
    if tr > 0;
        v = [0 0 0]';
    % case norm(om)=pi;
    else
        v = theta * (sqrt((diag(R)+1)/2).*[1;2*(R(1,2:3)>=0)'-1]);
    end;
end;

