import numpy as np
import matplotlib.pyplot as plt

m = 1
k = 1
w = np.sqrt(k/m)
t0 = 0
tf = 100
dt = 0.0001
n = int((tf-t0)/dt)
xeq = 1.5
Em = 0.75

t = np.linspace(t0, tf, n+1)

a = np.empty(n+1)

x = np.empty(n+1)
x[0] = np.sqrt(np.sqrt((2*Em)/k) + xeq**2)

v = np.empty(n+1)
v[0] = 0

for i in range(n):
    a[i] = (-2*k*(x[i]**2-xeq**2)*x[i])/m

    v[i+1] = v[i] + a[i]*dt

    x[i+1] = x[i] + v[i+1]*dt

plt.plot(t, x)
plt.show()
