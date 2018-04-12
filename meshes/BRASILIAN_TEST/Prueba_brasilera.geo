L = 0.1;

// Puntos
Point(1) = {0, 0, 0, L};
Point(2) = {1, 0, 0, L};
Point(3) = {0, 1, 0, L};

// Líneas
Line(1) = {3, 1};
Line(2) = {1, 2};
Circle(3) = {2, 1, 3};

// Superficie
Line Loop(1) = {2, 3, 1};
Plane Surface(1) = {1};

// Grupos Físicos
Physical Line(100) = {1};
Physical Line(200) = {2};
Physical Point(3) = {3};
Physical Surface(1000) = {1};
