#113786

import numpy as np
import matplotlib.pyplot as plt

dt=0.0001
tf=100
n=int(tf/dt+0.1)

t=np.linspace(0,tf,n)

#Dados do problema
xi = 0 #m
vi = 0 #m/s
k = 1 #N/m
m = 1 #kg
xeq = 0
b = 0.02 #kg/s
f0 = 7.5 #N
wf = 1 #rad/s

v = np.empty(n)
x = np.empty(n)
a = np.empty(n)

v[0] = vi
x[0] = xi

for i in range(n-1):
    a[i] = (-k/m)*x[i] - (b/m)*v[i] + (f0/m)*np.cos(wf*t[i])    #Acelaração
    v[i+1] = v[i] + a[i] * dt   #Velocidade
    x[i+1] = x[i] + v[i+1]*dt   #Posição
    
plt.figure()
plt.plot(t,x)
plt.title('2a)')
plt.ylabel('Posição(m)')
plt.xlabel( 't (s)' )
plt.grid()

plt.show()