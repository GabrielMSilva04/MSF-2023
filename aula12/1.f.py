import numpy as np
import matplotlib.pyplot as plt

t0 = 0
tf = 400
dt = 0.01
n = int((tf-t0)/dt)
k = 1
b = 0.05
F0 = 7.5
m = 1
w = 1
countMax = 0

t = np.linspace(t0, tf, n)

a = np.empty(n)

v = np.empty(n)
v[0] = -4

x = np.empty(n)
x[0] = -2

ind = np.transpose([0 for i in range(1000)])

afo = np.zeros(15)
bfo = np.zeros(15)


def abfourier(tp, xp, it0, it1, nf):
    dt = tp[1]-tp[0]
    per = tp[it1]-tp[it0]
    ome = 2*np.pi/per

    s1 = xp[it0]*np.cos(nf*ome*tp[it0])
    s2 = xp[it1]*np.cos(nf*ome*tp[it1])
    st = xp[it0+1:it1]*np.cos(nf*ome*tp[it0+1:it1])
    soma = np.sum(st)

    q1 = xp[it0]*np.sin(nf*ome*tp[it0])
    q2 = xp[it1]*np.sin(nf*ome*tp[it1])
    qt = xp[it0+1:it1]*np.sin(nf*ome*tp[it0+1:it1])
    somq = np.sum(qt)

    intega = ((s1+s2)/2+soma)*dt
    af = 2/per*intega
    integq = ((q1+q2)/2+somq)*dt
    bf = 2/per*integq
    return af, bf


for i in range(n-1):
    a[i] = (-k*x[i])/m - (b*v[i])/m + (F0*np.cos(w*t[i]))/m

    v[i+1] = v[i] + a[i]*dt

    x[i+1] = x[i] + v[i+1]*dt

    if t[i] > 200 and x[i-1] < x[i] and x[i+1] < x[i]:
        countMax = countMax + 1
        ind[countMax] = int(i)

t0 = ind[countMax-1]
t1 = ind[countMax]

for i in range(15):
    af, bf = abfourier(t, x, t0, t1, i)
    afo[i] = af
    bfo[i] = bf

ii = np.linspace(0, 14, 15)
plt.figure()
plt.ylabel("[ an ]")
plt.xlabel("n")
plt.bar(ii, np.abs(afo))
plt.grid()
plt.show()
