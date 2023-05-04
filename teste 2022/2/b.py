import numpy as np
import matplotlib.pyplot as plt

y0 = 800
g = 9.8
vt = 60
tf = 100
dt = 0.001

n = int(tf/dt)

t = np.empty(n)
vy = np.empty(n)
ay = np.empty(n)
y = np.empty(n)

t[0] = 0
vy[0] = 0
y[0] = y0

for i in range(n-1):
    t[i+1] = t[i]+dt
    if t[i] > 10:
        vt = 5
        vv = np.abs(vy[i])
        dres = g/vt**2
        ay[i] = -g-dres*vv*vy[i]
        vy[i+1] = vy[i]+ay[i]*dt
    else:
        vv = np.abs(vy[i])
        dres = g/vt**2
        ay[i] = -g-dres*vv*vy[i]
        vy[i+1] = vy[i]+ay[i]*dt
    y[i+1] = y[i]+vy[i]*dt

for i in range(n-1):
    if y[i] < 0+0.01 and y[i+1] > 0-0.01:
        print(
            f'Chegada ao solo - tempo: {t[i+1]}; posição: {y[i+1]}; velocidade (m/s): {vy[i+1]};')


plt.plot(t, y)
plt.xlabel('Tempo (s)')
plt.ylabel('Posição (m)')
plt.title('Queda do paraquedista')
plt.grid()
plt.show()