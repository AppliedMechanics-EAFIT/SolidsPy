/*
Solucion de Flamant-semiespacio a 45 grados con carga puntual
Autor: Juan Gomez
Fecha: Mayo 15, 2017
*/
lado = 1.0;  // lado caras inclinadas

// Puntos
Point(1) = {0, -lado, 0, 1.0};  // Centro
Point(2) = {0, 0, 0, 1.0};      // Vertice de la cunia
Point(3) = {-lado, -lado, 0, 1.0};
Point(4) = {lado, -lado, 0, 1.0};

//Lineas
Line(1) = {4, 2};
Line(2) = {2, 3};
Circle(3) = {3, 1, 4};

// Loops y superficies
Line Loop(4) = {1, 2, 3};
Plane Surface(5) = {4};

//Subdivision de superficies
Physical Surface(10000) = {5};
Transfinite Surface {5};

// Linea fisica para aplicar condiciones de frontera
Physical Line(100) = {3};