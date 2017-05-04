#Matrices de masa, amortiguamiento y rigidez
import numpy as np

n = int(input("Digite numero de estratos: ")) #Numero de estratos

if (n<5):
    print ("No se puede usar el programa")

r=n+1

#Matrices locales de masa

MatMlocal = np.zeros((n,2,2))#Matriz de masa local
MatPlocal = np.zeros((n,2,2))#Matrices locales de rigidez

g=9.8#m/s^2
rho=18000/g
Vp=100#m/s

H=0.0
LongElastico=0.0
cont=0.0

for i in range (n):
    if (i<(n-4)):
        l=float(input("Digite longitud del estrato:"))
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
            l=float(input("Digite longitud del estrato:"))
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
            l=float(input("Digite longitud del estrato:"))
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
#Ensamlbe general
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
fc = 4
Ts = 1
Vs = Vp
alpha= ((np.pi*fc)/(Vs-Ts))**2

#Inicilizar vectores de desplazmiento, de velocidad y aceleracion

U1=np.zeros(r)
UVel1=np.zeros(r)
UAce1=np.zeros(r)
U2=np.zeros(r)
UVel2=np.zeros(r)
UAce2=np.zeros(r)


#le es el valor de la longitud de ese estrato donde agg el valor de las fuerzas

delt=0
for i in range (10):
    delt=delt + 0.1
    fuerza=np.zeros(r)
    z1 = H+le
    z2 = H
    rest= z1
    col1=((4*alpha)*(np.exp(-(alpha*-(rest**2)))))*((((2*alpha)*-rest)**2)-(2*alpha))
    DerRicker1= ((4*alpha)*(np.exp(-(alpha*(rest**2)))))*((((2*alpha)*rest)**2)-(2*alpha))+col1
    y1 = (((2*((np.pi*fc)**2)*(((z1)/(Vs-Ts))**2))-1))*(np.exp(-((np.pi*fc)**2)*(((z1)/(Vs-Ts))**2))) + (((2*((np.pi*fc)**2)*(((-z1)/(Vs-Ts))**2))-1))*(np.exp(-((np.pi*fc)**2)*(((-z1)/(Vs-Ts))**2)))*(1/(mvarr*le))-(((rho*le)/6)*DerRicker1)
    rest= z2
    col2=((4*alpha)*(np.exp(-(alpha*-(rest**2)))))*((((2*alpha)*-rest)**2)-(2*alpha))
    DerRicker2= ((4*alpha)*(np.exp(-(alpha*(rest**2)))))*((((2*alpha)*rest)**2)-(2*alpha))+col2     
    y2 = (((2*((np.pi*fc)**2)*(((z2)/(Vs-Ts))**2))-1))*(np.exp(-((np.pi*fc)**2)*(((z2)/(Vs-Ts))**2))) + (((2*((np.pi*fc)**2)*(((-z2)/(Vs-Ts))**2))-1))*(np.exp(-((np.pi*fc)**2)*(((-z2)/(Vs-Ts))**2)))*(1/(mvaba*le))+(((rho*le)/6)* DerRicker2)
    fuerza[r-5] = y2
    fuerza[r-4] = y1
          
#Constantes de integracion
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
    print ("Matriz de rigidez efectiva", Keff, "en el tiempo:" , delt, "segundos")
    
#Matriz de fuerza efectiva de los nodos r-5 y r-4
    sum11= (a0*U1)+(a2*UVel1)+(a3*UAce1)
    sum12= (a1*U1)+(a4*UVel1)+(a5*UAce1)
    sum21= (a0*U2)+(a2*UVel2)+(a3*UAce2)
    sum22= (a1*U2)+(a4*UVel2)+(a5*UAce2)
    print (sum11)
    
    colt= (np.dot(MatM,sum11))
    Feff1=fuerza  + colt + (np.dot(MatC, sum12))
    Feff2=fuerza  + (np.dot(MatM, sum21)) + (np.dot(MatC, sum22))
    print (Feff1)
    
    Kefft=np.transpose(Keff)
    Udel1= np.dot(Kefft,Feff1)
    print (Udel1)
    
    Udel2= np.dot(Kefft,Feff2)
    UAcedel1= ((a0*(Udel1-U1))-(a2*UVel1)-(a3*UAce1))
    UAcedel2= ((a0*(Udel2-U2))-(a2*UVel2)-(a3*UAce2))
    UVeldel1= UVel2+(a6*UAce1)+(a7*UAcedel1)
    UVeldel2= UVel2+(a6*UAce2)+(a7*UAcedel2)
    
    U1=Udel1
    UVel1=UVeldel1
    UAce1=UAcedel1
    U2=Udel2
    UVel2=UVeldel2
    UAce2=UAcedel2