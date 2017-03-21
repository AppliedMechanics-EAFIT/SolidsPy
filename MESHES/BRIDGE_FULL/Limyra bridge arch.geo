/* Bridge near Limyra

References
----------
https://en.wikipedia.org/wiki/Bridge_near_Limyra
https://en.wikipedia.org/wiki/File:Limyra_Bridge_Arch.svg
*/

// Parameters of the bridge
radius = 8.08891;
span = 12.75;
clear_span = 10.65;
rise = 2;
thickness = 1.25;

// Parameters of the mesh
npts_arc1 = 10;
npts_arc2 = 8;
npts_rad = 6;

// Points
Point(1) = {0, rise  - radius, 0};
Point(2) = {0, rise + thickness, 0};
Point(3) = {span/2, rise + thickness, 0};
Point(4) = {span/2, 0, 0};
Point(5) = {clear_span/2, 0, 0};
Point(6) = {0, rise, 0};
Point(7) = {0.5*radius, rise + (Sqrt(3)/2 - 1)*radius, 0};
Point(8) = {-span/2, rise + thickness, 0};
Point(9) = {-span/2, 0, 0};
Point(10) = {-clear_span/2, 0, 0};
Point(11) = {-0.5*radius, rise + (Sqrt(3)/2 - 1)*radius, 0};

// Lines
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {6, 2};
Line(5) = {7, 3};
Circle(6) = {5, 1, 7};
Circle(7) = {7, 1, 6};
Line(8) = {2, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {11, 8};
Circle(12) = {10, 1, 11};
Circle(13) = {11, 1, 6};

// Surfaces
Line Loop(14) = {-1, 5, -7, -4};
Plane Surface(15) = {14};
Line Loop(16) = {-6, -5, -2, -3};
Plane Surface(17) = {16};
Line Loop(18) = {8, -11, 13, 4};
Plane Surface(19) = {18};
Line Loop(20) = {11, 9, 10, 12};
Plane Surface(21) = {20};

// Meshing properties
Transfinite Line {1, 7, 8, 13} = npts_arc1 Using Progression 1;
Transfinite Line {2, 6, 9, 12} = npts_arc2 Using Progression 1;
Transfinite Line {3, 4, 5, 10, 11} = npts_rad Using Progression 1;
Transfinite Surface {15};
Transfinite Surface {17};
Transfinite Surface {19};
Transfinite Surface {21};
Recombine Surface {15};
Recombine Surface {17};
Recombine Surface {19};
Recombine Surface {21};

// Physical groups
Physical Line("Top") = {1, 8};
Physical Line("Right") = {2};
Physical Line("Left") = {9};
Physical Line("Bottom") = {3, 10};
Physical Surface("Body") = {15, 17, 19, 21};
