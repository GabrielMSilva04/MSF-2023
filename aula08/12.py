import numpy as np
import matplotlib.pyplot as plt

P = 0.4 * 735.4975
vi = 1
m = 75
Cres = 0.9
A = 0.3
p = 1.225
u = 0.004
g = 9.8

ti = 0
tf = 200
dt = 0.0001

n = int((tf - ti) / dt)
t = np.linspace(ti, tf, n)

x = np.empty(n)
x[0] = 0

vx = np.empty(n)
vx[0] = vi

v = np.empty(n)
v[0] = vi

ax = np.empty(n)
ax[0] = 0

for i in range(n - 1):
    v[i] = np.sqrt(vx[i] ** 2)
    ax[i] = -u * g - (Cres / (2 * m)) * A * p * v[i] * vx[i] + (P / (m * v[i]))
    vx[i + 1] = vx[i] + ax[i] * dt
    x[i + 1] = x[i] + vx[i] * dt
    if vx[i] - vx[i - 1] < 0.0000001:
        vt = vx[i]
    if 1999 < x[i] < 2000:
        t2000 = t[i]
print("Tempo aos 2km:", t2000)

print("Velocidade terminal: ", vt)

v[n - 1] = np.sqrt(vx[n - 1] ** 2)

for i in range(n - 1):
    if (0.90 * vt - dt) < v[i] < (0.90 * vt + dt):
        print("Tempo: ", t[i])
        break

plt.plot(t, v)
plt.xlabel("Tempo (s)")
plt.ylabel("Velocidade (m/s)")
plt.grid()
plt.show()
