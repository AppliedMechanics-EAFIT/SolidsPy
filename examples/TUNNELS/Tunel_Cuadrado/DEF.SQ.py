

"""
Se crea un codigo para que modifique el archivo llamado nodes.txt creado en GMESH 
y compilado por medio de mesher.for que contiene los nodos, sus coordenadas 
en X y Y, y sus reacciones(inicialmente, no hay reacciones cada punto esta
representado por "0"); el cual ubica las reacciones en cada punto colocandole un
"-1". Adicional a esto, se crea un archivo que contiene los nudos con las cargas en X y/o Y,
y un archivo con informacion del material.
"""

import numpy as np
H=int(input('ingrese longitud del suelo: '))
a=np.loadtxt("nodes.txt")
Nelem=len(a[:,1])
Cx= np.zeros(Nelem)
Cy= np.zeros(Nelem)

for i in range(Nelem):
    Cx[i]=a[i,1]
    Cy[i]= a[i,2]
    
print ("1: Carga Distribuida")
print ("2: Carga Puntual")
d= int (input('Tipo de Carga:'))
w= int(input ('Ingrese magnitud de la carga: '))


## REACCIONES
Rx=np.zeros(Nelem)
Ry=np.zeros(Nelem)
m=0
M=0

#print 'Reacciones x'
while M+1 <=Nelem:
    if Cx[M]==H or Cx[M]==0:
        i=-1
            
    else:
        i=0

    Rx[m]=i
    m=m+1
    M=M+1
#print Rx


#print 'Reacciones y'
M=0
m=0
while M+1<=Nelem:
    if Cy[M]==0:
        i=-1
    else:
        i=0
        
    Ry[m]=i
    m=m+1
    M=M+1
#print Ry

#CARGAS
CargaX=np.zeros(Nelem)
CargaY=np.zeros (Nelem)
Nodo= np.zeros(Nelem)

    
if d == 2:
    for i in range (Nelem):
        
        if Cx[i]==H/2 and Cy[i]==H/2:
            CargaY[i]= w
        else:
            CargaY [i]=0
            #print CargaY
            
        Nodo[i]=i
        #print Nodo


if d==1: 
    i=0
    for i in range (Nelem):
        
        if Cy[i]==H/2:
            CargaY[i]= w
            
        else:
            CargaY[i]=0
            #print CargaY
           
        Nodo[i]= i   
        #print Nodo 
        
#Matriz nodes 
i=0
j=0
M= np.zeros([Nelem,5])
V= np.zeros((5*Nelem))
V[0:Nelem]= Nodo[0:Nelem]
V[Nelem:(2*Nelem)]= Cx[0: Nelem]
V[(2*Nelem):(3*Nelem)]= Cy[0:Nelem]
V[(3*Nelem):(4*Nelem)] = Rx[0:Nelem]
V[(4*Nelem):(5*Nelem)]= Ry[0:Nelem]

m=0
for j in range (5):
    for i in range (Nelem):
        M[i,j]= V[m]
        m=m+1
np.savetxt ("nodes.txt",M,'%3i %8.3f %8.3f %3i %3i \n')

#Matriz cargas 
i=0
j=0

C= np.zeros([Nelem,3])
V1= np.zeros((3*Nelem))
V1[0:Nelem]=Nodo[0:Nelem]
V1[Nelem:(2*Nelem)]= CargaX[0:Nelem]
V1[(2*Nelem):(3*Nelem)]=CargaY[0:Nelem]
m=0
for j in range (3):
    for i in range(Nelem):
        C[i,j]= V1[m]
        m=m+1
np.savetxt("loads.txt",C,'%3i %3i %3i')

"""
PROGRAMA QUE PERMITE ENCONTRAR LOS MODULOS DE YOUNG Y EL COEFICIENTE DE 
POISSON VIRTUALES A PARTIR DE LOS REALES
"""
import numpy as np

N=int(input('Numero de Materiales:')) 
if N==2:
    
    E1=int(input('Modulo de Young Material 1:'))
    print('TENGA EN CUENTA --> Coeficiente de Poisson entre 0.0 y 0.5')
    v1=float(input('Coeficiente de Poisson Material 1:'))
    if v1>=0 and v1<=1:
       E2=int(input('Modulo de Young Material 2:'))
       v2=float(input('Coeficiente de Poisson Material 2:'))
       if v2>=0 and v2<=1:
         EV1=(E1/(1-(v1)**2))
         EV2=(E2/(1-(v2)**2))
         vv1=(v1/(1-v1))
         vv2=(v2/(1-v2))
         print (EV1,vv1)
         print (EV2,vv2)
         M=np.zeros([2,2])
         A=[EV1,EV2,vv1,vv2]
         m=0
         for j in range (2):
           for i in range (2):
               M[i,j]=A[m]
               m=m+1             
         np.savetxt('mater.txt',M,'%8.3f %8.3f')
       else: 
            print('Coeficiente de Poisson debe estar entre 0 y 0.5')
    else:
        print('Coeficiente de Poisson debe estar entre 0 y 0.5')

else:
   
    E1=int(input('Modulo de Young Material:'))
    print('TENGA EN CUENTA --> Coeficiente de Poisson entre 0.0 y 0.5')
    v1=float(input('Coeficiente de Poisson Material:'))
    if v1>=0 and v1<=1:
        EV1=(E1/(1-(v1)**2))
        vv1=(v1/(1-v1))
        print (EV1,vv1)
        print (EV1,vv1)
        N=np.zeros([2,2])
        B=[EV1,EV1,vv1,vv1]
        n=0
        for j in range (2):
            for i in range (2):
               N[i,j]=B[n]
               n=n+1
        np.savetxt('mater.txt',N,'%8.3f %8.3f')
    else:
        print('Coeficiente de Poisson debe estar entre 0 y 0.5')
    
#    SE IMPRIME DOS VECES EL MISMO PORQUE EN EL ARCHIVO DE mater.txt DEBEN
#    HABER MINIMO 2 FILAS.  



