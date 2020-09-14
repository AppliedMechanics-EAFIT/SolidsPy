/*
Square membrane.

Units in m.
*/

side = Pi;

// Points
Point(1) = {0, 0, 0, 1};
Point(2) = {side, 0, 0, 1};
Point(3) = {side, side, 0, 1};
Point(4) = {0, side, 0, 1};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Physical groups
Physical Curve(1) = {1, 2, 3, 4};  // Top line
Physical Surface(2) = {1};  // Surface

// Mesh parameters
ndiv = 10;
Transfinite Curve {4, 2, 3, 1} = ndiv Using Progression 1;
Transfinite Surface {1};
Mesh.ElementOrder = 2;
