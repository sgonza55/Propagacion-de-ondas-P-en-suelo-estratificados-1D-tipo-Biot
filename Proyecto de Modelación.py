#Matrices de masa, amortiguamiento y rigidez
import numpy as np
import matplotlib.pyplot as plt

#n = int(input("Digite numero de estratos: ")) #Numero de estratos
n =20

if (n<4):
    print ("No se puede usar el programa")

r=n+1

#Matrices locales de masa

MatMlocal = np.zeros((n,2,2))#Matriz de masa local
MatPlocal = np.zeros((n,2,2))#Matrices locales de rigidez

g=9.81#m/s^2
rho=1800/g
Vp=100#m/s

H=0.0
LongElastico=0.0
cont=0.0

for i in range (n):
    if (i<(n-3)):
        #l=float(input("Digite longitud del estrato:"))
        l=1
        H=H+l
        mv=(1/(rho*(Vp**2)))
        krig=(1/(mv*l))
        coef=rho*l
        for j in range (2):
            for k in range (2):
                    if (j == k): 
                        MatPlocal[i,j,k]= krig
                        MatMlocal[i,j,k]= coef/3
                    else:
                        MatPlocal[i,j,k]= -krig
                        MatMlocal[i,j,k]= coef/6
        mvarr=mv
    else:
        if (cont==0):
            cont=cont+1
            #l=float(input("Digite longitud del estrato:"))
            le=l
            LongElastico= LongElastico + l
            mv=(1/(rho*(Vp**2)))
            krig=(1/(mv*l))
            coef=rho*l
            for j in range (2):
                for k in range (2):
                        if (j == k): 
                            MatPlocal[i,j,k]= krig
                            MatMlocal[i,j,k]= coef/3
                        else:
                            MatPlocal[i,j,k]= -krig
                            MatMlocal[i,j,k]= coef/6
            mvaba=mv
        else:
            #l=float(input("Digite longitud del estrato:"))
            LongElastico= LongElastico + l
            mv=(1/(rho*(Vp**2)))
            krig=(1/(mv*l))
            coef=rho*l
            for j in range (2):
                for k in range (2):
                        if (j == k): 
                            MatPlocal[i,j,k]= krig
                            MatMlocal[i,j,k]= coef/3
                        else:
                            MatPlocal[i,j,k]= -krig
                            MatMlocal[i,j,k]= coef/6
            
        
#print (MatMlocal)
#print (MatPlocal)
print ("La longitud total del suelo saturado es:", H)
print ("La longitud total del suelo elastico es:", LongElastico)

#---------------------------------------------------------------------------
#Ensamlbe general de masa
MatM = np.zeros((r,r))
l=0
jant=0
kant=0
tempant=0
for i in range (0,n): 
    for j in range (0,2):
        for k in range(0,2):
            if (j==k):
                temp1= MatMlocal[i,j,k]
                if (jant==j+l and kant==k+l):
                    MatM[j+l,k+l]=temp1+tempant
            else:
                temp2= MatMlocal[i,j,k]
                MatM[j+l,k+l]= temp2           
    if (jant==0):
        jant= jant+ j +l
        kant= kant+k+l
        tempant=temp1
        l=l+1
    elif (jant!=0):
        jant=1+l
        kant=1+l
        tempant=temp1
        l=l+1
MatM[j+l-1,k+l-1]=MatMlocal[i,j,k]
print ("Matriz de masa",MatM)
#-----------------------------------------------------------------------------
#Ensamlbe general de rigidez
MatP = np.zeros((r,r))
l=0
jant=0
kant=0
tempant=0
for i in range (0,n): 
    for j in range (0,2):
        for k in range(0,2):
            if (j==k):
                temp1= MatPlocal[i,j,k]
                if (jant==j+l and kant==k+l):
                    MatP[j+l,k+l]=temp1+tempant
            else:
                temp2= MatPlocal[i,j,k]
                MatP[j+l,k+l]= temp2           
    if (jant==0):
        jant= jant+ j +l
        kant= kant+k+l
        tempant=temp1
        l=l+1
    elif (jant!=0):
        jant=1+l
        kant=1+l
        tempant=temp1
        l=l+1
MatP[j+l-1,k+l-1]=MatPlocal[i,j,k]
print ("Matriz de rigidez" , MatP)

#----------------------------------------------------------------------------
#Matriz de amortiguamiento


MatC=np.zeros((r,r))
delta=rho*Vp
for i in range (r):
    for j in range (r):
        if (i==r-1):
            if (j==r-1):
                MatC[i,j]= delta
  
print ("Matriz de amortiguamiento",MatC)
#------------------------------------------------------------------------------
#Constantes del metodo de integracion
beta=0.7
if (beta >= 0.5):
    alfa=0.25*((0.5+beta)**2)
    print ("beta vale:",beta)
    print ("alfa vale:",alfa)
else:
    print ("No se puede ejecutar la integracion, beta supera el rango")

den =rho
#El FC puede tomar el valor de 2 0 4, hay que preguntar bien, que valor ha de tomar
fc = 4
Ts = 1
Vs = Vp
Iter=int(input("Digite numero de iteraciones: "))

#Inicilizar vectores de desplazmiento, de velocidad y aceleracion

fuerza=np.zeros(r)
U1=np.zeros(r)
UVel1=np.zeros(r)
UAce1=np.zeros(r)
Uprim=np.zeros(Iter)
deltacum=np.zeros(Iter)
y2acum=np.zeros(Iter)

#le es el valor de la longitud de ese estrato donde agg el valor de las fuerzas
#Constantes de integracion

delt=0.01

a0= 1/(alfa*((delt)**2))
a1=beta/(alfa*delt)
a2=1/(alfa*delt)
a3=(1/((2*alfa)))-1
a4= (beta/alfa)-1
a5=(delt/2)*((beta/alfa)-2)
a6= delt*(1-beta)
a7=beta*delt

# Matriz de rigidez efectiva

per1= (a0*MatM)
per2= (a1*MatC)

Keff= MatP + per1 + per2

#Inversa de la matriz de rigidez

Kefft=np.linalg.inv(Keff)
 
#ful=alpha**2
#col1=  ((2*(np.exp(-(ful*(delt-z1)**2))))*((-8*ful*(delt**2))+(8*ful*delt*z1)-(4*ful*(delt**2))+(16*ful*delt*z1)-(12*ful*(z1**2))+(6*alpha)))
#col2=  ((4*ful*(delt**3))-(8*ful*z1*(delt**2))+(4*ful*delt*(z1**2))-(6*alpha*delt)-(4*ful*z1*(delt**2))+(8*ful*(z1**2)*delt)-(4*ful*(z1**3))+(6*alpha*z1))
#DerRicker1= col1 +col2
     
delt=0

for i in range (Iter):
    
    z1 = H+le
    z2 = H
    
#Calculo las fuerza por pulso de Ricker en cada delta de tiempo
    alfa11 = ((np.pi*fc)**2) * ((delt-(z1/Vs)-Ts)**2)
    alfa12 = ((np.pi*fc)**2) * ((delt+(z1/Vs)-Ts)**2)
    
    y1= ((((2*alfa11)-1) * (np.exp(-alfa11))) + (((2*alfa12)-1) * (np.exp(-alfa12)))) * ( 1/(mvaba*le))
    
    alfa21 = ((np.pi*fc)**2) * ((delt-(z2/Vs)-Ts)**2)
    alfa22 = ((np.pi*fc)**2) * ((delt+(z2/Vs)-Ts)**2)
    
    
    y2= ((((2*alfa21)-1) * (np.exp(-alfa21))) + (((2*alfa22)-1) * (np.exp(-alfa22)))) * (-1/(mvaba*le))
    
    fuerza[r-4] = y1
    fuerza[r-3] = y2
    
    y2acum[i]=y2
        
#Matriz de fuerza efectiva

    sum11= (a0*U1)+(a2*UVel1)+(a3*UAce1)
    sum12= (a1*U1)+(a4*UVel1)+(a5*UAce1)

    colt1= (np.dot(MatM, sum11))
    colt2= (np.dot(MatC, sum12))
    
    Feff1= fuerza + colt1 + colt2
    
    
    Udel1= np.dot(Kefft,Feff1)
    Uprim[i]=Udel1[0]
    
    UAcedel1= ((a0*(Udel1-U1))-(a2*UVel1)-(a3*UAce1))
    UVeldel1= UVel1+(a6*UAce1)+(a7*UAcedel1)

    U1=Udel1
    UVel1=UVeldel1
    UAce1=UAcedel1
    
    deltacum[i]=delt
    delt=delt + 0.01
    

uprim = open('uprim.txt', 'w')
for i in range(Iter):
    uprim.write(" %f %f %f " % (Uprim[i],y2acum[i],deltacum[i]))
    uprim.write("\n")    
uprim.close()

plt.figure()
plt.grid()
ax = plt.axes(xlim=(0, 10), ylim=(-60,40))
plt.plot(deltacum,Uprim)
plt.xlabel('$Tiempo$',fontsize=15)
plt.ylabel('$Desplazamientos-del-primer-nodo$',fontsize=15)
plt.show()
    
