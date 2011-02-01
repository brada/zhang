% rotation matrix for rotation by theta degrees around arbitrary axis v
function M = zhang_rotationmatrix( v, theta )

v = v / norm(v);
theta = deg2rad(theta);
c = cos(theta);
s = sin(theta);
t = 1 - cos(theta);
X = v(1);
Y = v(2);
Z = v(3);

M = [ t*X^2 + c,      t*X*Y + s*Z,    t*X*Z - s*Y;
      t*X*Y - s*Z,    t*Y^2 + c,      t*Y*Z + s*X;
      t*X*Y + s*Y,    t*Y*Z - s*X,    t*Z^2 + c ];
