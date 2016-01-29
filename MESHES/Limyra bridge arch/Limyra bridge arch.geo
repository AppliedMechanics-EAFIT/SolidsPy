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
Point(1) = {0, rise + thickness, 0};
Point(2) = {span/2, rise + thickness, 0};
Point(3) = {span/2, 0, 0};
Point(4) = {clear_span/2, 0, 0};
Point(5) = {0, rise, 0};
Point(6) = {0, rise  - radius, 0};
Point(7) = {0.5*radius, rise + (Sqrt(3)/2 - 1)*radius, 0};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {5, 1};
Line(5) = {7, 2};
Circle(6) = {4, 6, 7};
Circle(7) = {7, 6, 5};

// Surfaces
Line Loop(8) = {1, -5, 7, 4};
Plane Surface(9) = {8};
Line Loop(10) = {6, 5, 2, 3};
Plane Surface(11) = {10};
Transfinite Line {7, 1} = npts_arc1 Using Progression 1;
Transfinite Line {6, 2} = npts_arc2 Using Progression 1;
Transfinite Line {3, 4, 5} = npts_rad Using Progression 1;
Transfinite Surface {9};
Transfinite Surface {11};
Recombine Surface {9};
Recombine Surface {11};
