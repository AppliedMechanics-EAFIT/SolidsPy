// Input .geo for bimaterial
// author: Juan Gomez

c = 5.0; 						// for size elements

// Define vertex points 						
Point(1) = {0.0, 0.0, 0, c};				// {x,y,z, size}
Point(2) = {2.0, 0.0, 0, c};
Point(3) = {2.0, 1.0, 0, c};
Point(4) = {2.0, 2.0, 0, c};
Point(5) = {0.0, 2.0, 0, c};
Point(6) = {0.0, 1.0, 0, c};


// Define boundary lines
Line(1) = {1, 2};					// {Initial_point, end_point}
Line(2) = {2, 3};
Line(3) = {3, 6};
Line(4) = {6, 1};
Line(5) = {3, 4};
Line(6) = {4, 5};
Line(7) = {5, 6};

// Joint Lines
Line Loop(1) = {1, 2, 3, 4};			        // {Id_line1,id_line2, ... }
Line Loop(2) = {-3, 5, 6, 7};

// surface for mesh 					// {Id_Loop}
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// For Mesh 4 nodes
//Recombine Surface {1};				// {Id_Surface}

// "Structure" mesh 
Transfinite Surface {1,2};				// {Id_Surface}

Physical Surface(100) = {1};
Physical Surface(200) = {2};
