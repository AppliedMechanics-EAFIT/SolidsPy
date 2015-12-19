// Input .geo for dam
// author: Juan Gomez

c = 5.0; 						// for size elements

// Define vertex points 						
Point(1) = {-50.0,-150.0, 0, c};		        // {x,y,z, size}
Point(2) = {100.0,-150.0, 0, c};
Point(3) = {100.0, 0.000, 0, c};
Point(4) = {50.00, 0.000, 0, c};
Point(5) = {0.000, 0.000, 0, c};
Point(6) = {-50.0, 0.000, 0, c};
Point(7) = { 50.0, 100.0, 0, c};
Point(8) = { 25.0, 100.0, 0, c};


// Define boundary lines
Line(1) = {1, 2};					// {Initial_point, end_point}
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {4, 7};
Line(8) = {7, 8};
Line(9) = {8, 5};

// Joint Lines
Line Loop(1) = {1, 2, 3, 4, 5, 6};			 // {Id_line1,id_line2, ... }
Line Loop(2) = {-4, 7, 8, 9};

// surface for mesh 					// {Id_Loop}
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// For Mesh 4 nodes
//Recombine Surface {1};				// {Id_Surface}

// "Structure" mesh 
//Transfinite Surface {1,2};				// {Id_Surface}

// Physical surface. Two material 
Physical Surface(100) = {1};
Physical Surface(200) = {2};

//Physical line. Boundary 
//Physical Line(1000) = {1};
//Physical Line(2000) = {2};
//Physical Line(3000) = {4};
//Physical Line(4000) = {6};
//Physical Line(5000) = {7};
