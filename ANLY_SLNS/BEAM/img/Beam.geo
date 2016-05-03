// Gmsh project created on Tue May 03 12:45:55 2016
Point(1) = {0, 0, 0, 1.0};
Point(2) = {24, 0, 0, 1.0};
Point(3) = {24, 8, 0, 1.0};
Point(4) = {0, 8, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Transfinite Line {3, 1} = 4 Using Progression 1;
Transfinite Line {4, 2} = 2 Using Progression 1;
Transfinite Surface {6};
Recombine Surface {6};
