//TUNEL DOBLE 
c=2.5;
z=1.75;
//PUNTOS Y COORDENADAS

Point(1) = {0,0,0,c};
Point(2) ={0,100,0,c};
Point(3) ={200,100,0,c};
Point(4) ={200,0,0,c};
Point(5) ={105,50,0,c};
Point(6) ={109,50,0,z};
Point(7) ={124,50,0,c};
Point(8) ={57,50,0,c};
Point(9) ={61,50,0,z};
Point(10) ={76,50,0,c};

//LINEAS

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Circle(5) = {6,7,6};
Circle(6) = {5,7,5};
Circle(7) = {9,10,9};
Circle(8) = {8,10,8};

//LOOPS

Line loop(1) = {1,2,3,4};
Line loop(2) = {5};
Line loop(3) = {6};
Line loop(4) = {7};
Line loop(5) = {8};

//SURFACES

Plane Surface(1) = {1,3,5};
Plane Surface(2) = {3,2};
Plane Surface(3) = {5,4};

//PHYSICAL SURFACES
Physical Surface(10) = {1};
Physical Surface(20) = {2,3};