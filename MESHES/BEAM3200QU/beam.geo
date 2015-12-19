// Input .geo for beam
// author: Juan Gomez

c = 0.25; 									// for size elements

// Define vertex points 						
Point(1) = {0, -4, 0, c};				// {x,y,z, size}
Point(2) = {24, -4, 0, c};
Point(3) = {24,  4, 0, c};
Point(4) = {0, 4, 0, c};

// Define boundary lines
Line(1) = {1, 2};						// {Initial_point, end_point}
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Joint Lines
Line Loop(1) = {1, 2, 3, 4};			// {Id_line1,id_line2, ... }

// surface for mesh 					// {Id_Loop}
Plane Surface(1) = {1};

// For Mesh 4 nodes
Recombine Surface {1};					// {Id_Surface}

// "Structure" mesh 
Transfinite Surface {1};				// {Id_Surface}

Physical Surface(100) = {1};
