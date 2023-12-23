#113786

import math
import numpy as np
import matplotlib.pyplot as plt

m = 2000 # massa do corpo em kg
P = 40000 # potência em W
x0 = 0 # posição inicial
v0 = 1 # velocidade inicial
u = 0.04 # coeficiente de resistência do alcatrão
cres = 0.25 # coeficiente de resistência do ar
A = 2 # Area frontal
g=9.8

ti = 0
tf = 100
dt = 0.001
n = int((tf-ti)/dt)

t = np.linspace(ti,tf,n+1)

def planoInclinado_res_1D(x0,v0,n,dt,cres,u,A,m,P,ang):
    x=np.empty(n+1)
    vx=np.empty(n+1)
    ax=np.empty(n+1)
    
    p_ar=1.225
    
    x[0]=x0
    vx[0]=v0
    ax[0]=0
    
    for i in range(n):
        vv=np.abs(vx[i])
        ax[i]=-g*np.sin(ang) -u*g*np.cos(ang) -(0.5*cres*A*p_ar*vx[i]*vv)/m + P/(m*vx[i])
        vx[i+1]=vx[i]+ax[i]*dt
        x[i+1]=x[i]+vx[i]*dt
    return x,vx,ax

values = planoInclinado_res_1D(x0,v0,n,dt,cres,u,A,m,P,ang=5)
x = values[0]
v = values[1]

# 2a) Gráfico de x em função de t
plt.xlabel("t (s)")
plt.ylabel("x (m/s)")
plt.plot(t,x)
plt.grid()
plt.title("2a) Evolução da posição em função do tempo")
plt.show()

# 2a) Gráfico de v em função de t
plt.xlabel("t (s)")
plt.ylabel("v (m/s)")
plt.plot(t,v)
plt.grid()
plt.title("2a) Evolução da velocidade em função do tempo")
plt.show()

# 2b) Tempo para atingir 2000 m
# 2c) Trabalho realizado pelo motor
for i in range(n):
    if (2000-0.05)<x[i]<(2000+0.05):
        print("b)Tempo depois de 2km:", t[i], "s")
        print("c)Trabalho realizado pelo motor:", P*t[i], "J")
        break
