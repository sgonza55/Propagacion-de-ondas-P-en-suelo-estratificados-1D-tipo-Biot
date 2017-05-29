#%%#Librerias importadas necesarias para el codigo
#plt.gca().invert_yaxis()
import numpy as np
import matplotlib.pyplot as plt

#%%#Al multplicar por el delta de tiempo, nos da el tiempo total del sistema
Iter= 200 #int(input("Digite numero de iteraciones: "))

#%%#Numero de estratos del sistema
n = 6 #int(input("Digite numero de estratos: ")) 

#%%#Número de estratos de solo la zona elastica
ZoneElas= 3 #int(input("Digite numero de estratos de la zona elastica: ")) 

if (n<ZoneElas):
    print ("No se puede usar el programa")
    
#%%# Numero de divisiones
div=3 #int(input("Digite en cuanto desea partir el estrato")) 

#%%Número de nodos del suelo poroso
poros = div*(n-ZoneElas)

#%% Punto presion poros
soil = poros-1

#%%# ORDEN MATRICES 
C=(div*(n-ZoneElas))+ZoneElas
A=C+1 

mat = poros+A

#%%#Constantes del sistema
g=9.81#m/s^2
PesoEspecificoAgua=9800 #N/m^3

#%%#Vector donde se guardaran TODAS las velocidades de propagacion, Ss, e y l del modelo
Vpall=np.zeros(n)
Ssall=np.zeros(n)
eall=np.zeros(n)
lall=np.zeros(n)

#%%#Iniciazion de las matrices locales de masa y rigidez locales
MatMlocal = np.zeros((C,2,2))#Matriz de masa local
MatPlocal = np.zeros((C,2,2))#Matrices locales de rigidez
MatHlocal = np.zeros((poros,2,2))

#%%#Algoritmo que forma las matrices de masa y rigidez
H=0.0  #Longitud suelo saturado
LongElastico=0.0  #Longitud suelo elastico
cont=0.0
conta=0

for i in range (n-ZoneElas):
    Ss=2.7 #float(input("Digite peso especifico realtivo de los solidos: "))
    e= 0.9 #float(input("Digite relacion de vacios: "))
    Vp=100 #int(input("Digite velocidad de propagacion del estrato: "))
    l=9 #int(input("Digite longitud del estrato:"))
    Vpall[conta]=Vp
    Ssall[conta]=Ss
    eall[conta]=e
    lall[conta]=l
    conta=conta +1
ele=l/div
    
for i in range (ZoneElas):
    Ss=2.7 #float(input("Digite peso especifico realtivo de los solidos: "))
    e= 0.9 #float(input("Digite relacion de vacios: "))
    Vp=100 #int(input("Digite velocidad de propagacion del estrato: "))
    Vpall[conta]=Vp
    Ssall[conta]=Ss
    eall[conta]=e
    lall[conta]=ele
    conta=conta +1

u=0                   
for i in range (n):
    if (i<(n-ZoneElas)):
        Ss=Ssall[i]
        e=eall[i]
        Vp=Vpall[i]
        l=lall[i]
        GamaSat=PesoEspecificoAgua*((Ss+e)/(1+e))
        rho=GamaSat/g
        H=H+l
        for p in range (div):
            lediv=l/div
            mv=(1/(rho*(Vp**2)))
            krig=(1/(mv*lediv))
            coef=rho*lediv
            ache=krig/(PesoEspecificoAgua*lediv)
            for j in range (2):
                for k in range (2):
                        if (j == k): 
                            MatPlocal[u,j,k]= krig
                            MatMlocal[u,j,k]= coef/3
                            MatHlocal[u,j,k]= ache         
                        else:
                            MatPlocal[u,j,k]= -krig
                            MatMlocal[u,j,k]= coef/6
                            MatHlocal[u,j,k]= -ache         
            u=u+1
    else:
        Ss=Ssall[i]
        e=eall[i]
        Vp=Vpall[i]
        l=lall[i]
        GamaSat=PesoEspecificoAgua*((Ss+e)/(1+e))
        rho=GamaSat/g
        LongElastico= LongElastico + l
        mv=(1/(rho*(Vp**2)))
        krig=(1/(mv*l))
        coef=rho*l
        le=l
        for j in range (2):
            for k in range (2):
                    if (j == k): 
                        MatPlocal[u,j,k]= krig
                        MatMlocal[u,j,k]= coef/3
                    else:
                        MatPlocal[u,j,k]= -krig
                        MatMlocal[u,j,k]= coef/6
        u=u+1                         
        mvele=mv
        le=l

print ("La longitud total del suelo saturado es:", H)
print ("La longitud total del suelo elastico es:", LongElastico)

#---------------------------------------------------------------------------
#%%#Ensamlbe general de la matriz de masa
MatM = np.zeros((A,A))
l=0
jant=0
kant=0
tempant=0
for i in range (0,C): 
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

#%%#Ensamlbe GLOBAL de la matriz de masa
Mtot = np.zeros((mat,mat))

g= poros
for i in range (0,A):
    for j in range (0,A):
        if (i==j):
            Mtot[g,g]= MatM[i,j]
            g = g+1
t=0
for i in range (0,C):
    Mtot [poros+t, poros+t+1] = MatM[i,i+1]
    t = t+1
t=0
for i in range (0,C):
    Mtot [poros+t+1, poros+t] = MatM[i+1,i]
    t = t+1
    
#%%Matriz de amortiguamiento
MatC=np.zeros((A,A))
Ctot = np.zeros((mat,mat))
VP1=Vpall[n-1] 
delta=rho*VP1
for i in range (A):
    for j in range (A):
        if (i==C):
            if (j==C):
                MatC[i,j]= delta
  
print ("Matriz de amortiguamiento",MatC)

#%%Matriz Q
MatQ=np.zeros((poros+1,poros+1))

t=0
for i in range (0,poros):
    MatQ [t, t+1] = -0.5
    t = t+1
t=0
for i in range (0,poros):
    MatQ [t+1, t] = 0.5
    t = t+1    
    
for i in range (0,poros):
    for j in range (0, poros):
        if (i==j):
            if (i==0 and j==0):
                MatQ[i,j]= -0.5
            else:
                MatQ[i,j]=0
                                   
MatQ[poros,poros]=0.5 

#%% Matriz Q1    
MatQ1=np.zeros((poros+1,poros))
    
t=0
for i in range (0,poros):
    MatQ1 [t,t] = MatQ[t,t+1]
    t = t+1

t=0
for i in range (0,soil):
    MatQ1 [t+2,t] = MatQ[t+2,t+1]
    t = t+1
    
MatQ1[poros,soil]=MatQ[poros,poros]
   
#Transpuesta de la matriz Q     
MatQT=np.transpose(MatQ1)

#%%Ensamlbe GLOBAL de la matriz de amortiguamiento
Ctot = np.zeros((mat,mat))
            
t=0
for i in range (0,poros):
    Ctot [t, poros+t] = MatQT[t,t]
    t = t+1

t=0
for i in range (0,soil):
    Ctot [t, poros+t+2] = MatQT[t,t+2]
    t = t+1

Ctot [soil, (2*poros)] = MatQT[soil,poros]
    
g= poros
for i in range (0,A):
    for j in range (0,A):
        if (i==j):
            Ctot[g,g]= MatC[i,j]
            g = g+1
t=0
for i in range (0,C):
    Ctot [poros+t, poros+t+1] = MatC[i,i+1]
    t = t+1
t=0
for i in range (0,C):
    Ctot [poros+t+1, poros+t] = MatC[i+1,i]
    t = t+1
    
#-----------------------------------------------------------------------------
#%%#Ensamlbe general de la matriz de rigidez
MatP = np.zeros((A,A))
l=0
jant=0
kant=0
tempant=0
for i in range (0,C): 
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

#%%#Ensamble de la matriz H
MatH = np.zeros((poros+1,poros+1))
l=0
jant=0
kant=0
tempant=0
for i in range (0,poros): 
    for j in range (0,2):
        for k in range(0,2):
            if (j==k):
                temp1= MatHlocal[i,j,k]
                if (jant==j+l and kant==k+l):
                    MatH[j+l,k+l]=temp1+tempant
            else:
                temp2= MatHlocal[i,j,k]
                MatH[j+l,k+l]= temp2           
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
MatH[j+l-1,k+l-1]=MatHlocal[i,j,k]

#%%#Ensamlbe GLOBAL de la matriz de rigidez
Ptot = np.zeros((mat,mat))

for i in range (0,poros):
    for j in range (0,poros):
        if (i==j):
            Ptot[i,i]= MatH[i+1,j+1]
t=0
for i in range (0,soil):
    Ptot [t,t+1] = MatH[t+1,t+2]
    t = t+1

t=0
for i in range (0,soil):
    Ptot [t+1, t] = MatH[t+2,t+1]
    t = t+1
    
t=0
for i in range (0,soil):
    Ptot [poros+2+t, t] = -MatQ[t+2,t+1]
    t = t+1
    
t=0
for i in range (0,poros):
    Ptot [poros+t, t] = -MatQ[t,t+1]
    t = t+1

Ptot [(poros*2),soil]= -MatQ[poros,poros]

g= poros
for i in range (0,A):
    for j in range (0,A):
        if (i==j):
            Ptot[g,g]= MatP[i,j]
            g = g+1
t=0
for i in range (0,C):
    Ptot [poros+t, poros+t+1] = MatP[i,i+1]
    t = t+1
t=0
for i in range (0,C):
    Ptot [poros+t+1, poros+t] = MatP[i+1,i]
    t = t+1
#------------------------------------------------------------------------------
#%%#PROCESO DE METODO DE ELEMENTOS FINITOS
#------------------------------------------------------------------------------
#%%#Constantes del metodo de integracion
beta= 0.7 #float(input("Digite constante beta de integración: "))
if (beta >= 0.5):
    alfa=0.25*((0.5+beta)**2)
    print ("beta vale:",beta)
    print ("alfa vale:",alfa)
else:
    print ("No se puede ejecutar la integracion, beta supera el rango")

#Delta de tiempo del sistema (incremento constante del tiempo)
delt=0.01

#%%#Constantes de integracion
a0= 1/(alfa*((delt)**2))
a1=beta/(alfa*delt)
a2=1/(alfa*delt)
a3=(1/((2*alfa)))-1
a4= (beta/alfa)-1
a5=(delt/2)*((beta/alfa)-2)
a6= delt*(1-beta)
a7=beta*delt

#%%#Inicilizar vectores de fuerza, presion de poros, tiempo, desplazmiento, de velocidad y aceleracion del sistema
fuerza=np.zeros(mat)
deltacum=np.zeros(Iter)
U1=np.zeros(mat)
UVel1=np.zeros(mat)
UAce1=np.zeros(mat)
Uprim=np.zeros((mat-poros+1,Iter))

#%%# Matriz de rigidez efectiva
per1= (a0*Mtot)
per2= (a1*Ctot)
Keff= Ptot + per1 + per2

#%%#Inversa de la matriz de rigidez
Kefft=np.linalg.inv(Keff)

#%%#Datos iniciales para empezar el proceso de iteracion del FEM     
delt=0 #Tiempo inicial
fc = 4 #Frecuancia
Ts = 1 #Punto maximo del pulso de Ricker
 
#%%#Proceso de iteración del sistema
for i in range (Iter):
    
    z1 = H+le
    z2 = H
    T= 2*poros
    K= T+1
    
#%%Calculo las fuerza por pulso de Ricker en cada delta de tiempo
    Vs=Vpall[(ZoneElas-1)]
    alfa11 = ((np.pi*fc)**2) * ((delt-(z1/Vs)-Ts)**2)
    alfa12 = ((np.pi*fc)**2) * ((delt+(z1/Vs)-Ts)**2)
    y1= ((((2*alfa11)-1) * (np.exp(-alfa11))) + (((2*alfa12)-1) * (np.exp(-alfa12)))) * ( 1/(mvele*le))
    
    Vs=Vpall[ZoneElas-1]
    alfa21 = ((np.pi*fc)**2) * ((delt-(z2/Vs)-Ts)**2)
    alfa22 = ((np.pi*fc)**2) * ((delt+(z2/Vs)-Ts)**2)
    y2= ((((2*alfa21)-1) * (np.exp(-alfa21))) + (((2*alfa22)-1) * (np.exp(-alfa22)))) * (-1/(mvele*le))

#%%Ubicar la aplicacion de las fuerzas en los nodos correspondientes    
    fuerza[T] = y1
    fuerza[K] = y2

#%%#Matriz de fuerza efectiva
    sum11= (a0*U1)+(a2*UVel1)+(a3*UAce1)
    sum12= (a1*U1)+(a4*UVel1)+(a5*UAce1)
    colt1= (np.dot(Mtot, sum11))
    colt2= (np.dot(Ctot, sum12))
    Feff1= fuerza + colt1 + colt2

#%%#Calculo del vector de desplazamientos en el nuevo tiempo   
    Udel1= np.dot(Kefft,Feff1)
    
#%%#Calculo de los vectores de velocidad y aceleracion en el nuevo tiempo
    UAcedel1= ((a0*(Udel1-U1))-(a2*UVel1)-(a3*UAce1))
    UVeldel1= UVel1+(a6*UAce1)+(a7*UAcedel1)

#%%#Vectores que registran el desplzamiento de un nodo y de cada tiempo respectivamente
    
#%%#Vectores que registran el desplzamiento de un nodo y de cada tiempo respectivamente
    for j in range(0,mat-poros):
        Uprim[j,i]=Udel1[j+poros]
    
    deltacum[i]=delt


#%%#Iniicalizar nuevamente las variables para la siguiente iteracion
    U1=Udel1
    UVel1=UVeldel1
    UAce1=UAcedel1
    
#%%#Incremento de tiempo
    delt=delt + 0.01
 
#%%#FIN DEL PROCESO DE FEM
#------------------------------------------------------------------------------
#%%#Algoritmo que guarda desplazamientos de un nodo en formato .txt
uprim = open('uprim.txt', 'w')
for i in range(Iter):
    uprim.write(" %f %f " % (Uprim[i],deltacum[i]))
    uprim.write("\n")    
uprim.close()



#%%#Algoritmo que grafica los desplazamientos de cada nodo Vs el tiempo
dx=0
plt.figure()
plt.grid(True)
plt.grid(linewidth = 0.8)
for j in range (0,mat-poros-1):
    plt.plot(deltacum,(5*Uprim[j,:] + dx ), 'k',label="Desplazamientos")
    dx=dx-0.05
plt.xlabel('Tiempo',fontsize=15)
plt.ylabel('Desplazamientos',fontsize=15)
plt.show()

