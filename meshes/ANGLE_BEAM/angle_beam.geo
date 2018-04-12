// Gmsh project created on Thu Mar 16 11:18:30 2017

// Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {10, 4, 0, 1.0};
Point(3) = {10, 2, 0, 1.0};
Point(4) = {0, 4, 0, 1.0};

// Lines
Line(1) = {1, 3};
Line(2) = {3, 2};
Line(3) = {2, 4};
Line(4) = {4, 1};

// Surface 
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Transfinite Line {3, 1} = 21 Using Progression 1;
Transfinite Line {4, 2} = 11 Using Progression 1;
Transfinite Surface {6};
Recombine Surface {6};
