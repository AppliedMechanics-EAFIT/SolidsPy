size = 1;
Point(1) = {-75, 0, 0, size};
Point(2) = {-25, 0, 0, size};
Point(3) = {-25, -5, 0, size};
Point(4) = {25, -5, 0, size};
Point(5) = {25, 0, 0, size};
Point(6) = {75, 0, 0, size};
Point(7) = {-75, -5, 0, size};
Point(8) = {75, -5, 0, size};
Point(9) = {-75, -11, 0, size};
Point(10) = {75, -11, 0, size};
Point(11) = {-75, -15, 0, size};
Point(12) = {75, -15, 0, size};
Point(13) = {-75, -18, 0, size};
Point(14) = {75, -18, 0, size};
Point(15) = {-75, -33, 0, size};
Point(16) = {75, -33, 0, size};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {7, 8};
Line(6) = {9, 10};
Line(7) = {11, 12};
Line(8) = {13, 14};
Line(9) = {15, 16};
Line(10) = {1, 7};
Line(11) = {6, 8};
Line(12) = {7, 9};
Line(13) = {8, 10};
Line(14) = {9, 11};
Line(15) = {10, 12};
Line(16) = {11, 13};
Line(17) = {12, 14};
Line(18) = {13, 15};
Line(19) = {14, 16};
Line Loop(1) = {-5, 12, 6, -13};
Line Loop(2) = {-6, 14, 7, -15};
Line Loop(3) = {-7, 16, 8, -17};
Line Loop(4) = {-8, 18, 9, -19};
//Plane Surface
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
// Physical surface
Physical Surface(100) = {1};
Physical Surface(200) = {2};
Physical Surface(300) = {3};
Physical Surface(400) = {4};





