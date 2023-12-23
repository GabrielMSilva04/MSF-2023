import numpy as np
import matplotlib.pyplot as plt

k = 1
m = 1
b = 0.05
F0 = 7.5
Wf = 1.4
alpha = 0.001

def oscilador_quartico(k,m,b,F0,Wf,alpha,dt, tf):
    ti = 0
    n = int((tf-ti)/dt+0.1)
    tempo = np.empty(n+1)
    x = np.empty(n+1)
    vx = np.empty(n+1)
    a = np.empty(n+1)
    Em = np.empty(n+1)

    t0 = 0.
    x0 = 2.0
    vx0 = 4

    tempo[0] = t0
    vx[0] = vx0
    x[0] = x0

    k = 1
    m = 1
    b = 0.05
    F0 = 7.5
    Wf = 1.4
    alpha = 0.001
    ampl = 0
    countMax = 0
    tMax = []
    periodo = []
    for i in range(n):
        tempo[i+1] = tempo[i]+dt
        a[i] = -(k/m)*x[i]*(1+2*alpha*x[i]**2) - \
            (b/m)*vx[i]+(F0/m)*np.cos(Wf*tempo[i])
        vx[i+1] = vx[i]+a[i]*dt
        x[i+1] = x[i]+vx[i+1]*dt

    return a, vx, x, tempo


dt = 0.001
tf = 300
a_1, vx_1, x_1, t_1 = oscilador_quartico(k,m,b,F0,Wf,alpha,dt, tf)
dt = 0.01
a_2, vx_2, x_2, t_2 = oscilador_quartico(k,m,b,F0,Wf,alpha,dt, tf)

x_temp = x_1[t_1 > 150]
t_temp = t_1[t_1 > 150]
maximos_x = x_temp[:-2][np.diff(np.sign(np.diff(x_temp))) == -2]
maximos_t = t_temp[:-2][np.diff(np.sign(np.diff(x_temp))) == -2]
print("Amplitude:", np.round(np.mean(maximos_x), 3), "m")
print("Período:",  np.round(np.mean(np.diff(maximos_t)), 3), "s")