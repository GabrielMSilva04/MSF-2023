import numpy as np
import matplotlib.pyplot as plt

t0 = 0
tf = 20
dt = 0.001
n = int((tf-t0)/dt)
g = 9.8

t = np.linspace(t0, tf, n)

a = np.empty(n)
w = np.empty(n)
teta = np.empty(n)

teta[0] = np.radians(1)
w[0] = 0

tempos = []
for i in range(n-1):
    a[i] = (-g/2)*np.sin(teta[i])
    w[i+1] = w[i] + a[i]*dt
    teta[i+1] = teta[i] + w[i+1]*dt

    if teta[i-1] < teta[i] > teta[i+1]:
        tempos.append(t[i])


periodos = []
for i in range(len(tempos)-1):
    periodos.append(tempos[i+1]-tempos[i])

periodo = sum(periodos)/len(periodos)
print("Período: ", periodo)

plt.plot(t, teta)
plt.grid()
plt.show()
