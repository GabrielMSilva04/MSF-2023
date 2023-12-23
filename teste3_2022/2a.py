import numpy as np
import matplotlib.pyplot as plt

x0 = 0
b = 0.05
k =1
m = 1
#Fx=F0*cos(w_f*t)
w_f = 1.4
F0 = 7.5
x0 = 2
v0 = 4
ti = 0
tf = 400
dt = 0.01
n = int((tf-ti)/dt)

def oscSimpFA_1D(x0,v0,k,m,t,b,F0,w_f,n,dt):
    #oscilador Simples sujeito forçado e amortecido
    x=np.empty(n+1)
    v=np.empty(n+1)
    a=np.empty(n+1)
    x[0]=x0
    v[0]=v0
    for i in range(n):
        a[i]=-k/m*x[i]-b/m*v[i]+F0/m*np.cos(w_f*t[i])
        v[i+1]=v[i]+a[i]*dt
        x[i+1]=x[i]+v[i+1]*dt
    return x,v,a

def amp_per_comp(x, t, n, reg_est):
    #amplitde, periodo e comprimeno de onda com regime estacionario
    ind_max=[i for i in range(1,n-1) if x[i-1]<=x[i]>=x[i+1] if t[i]>reg_est]
    x_max=[x[i] for i in ind_max]
    t_max=[t[i] for i in ind_max]
    A=np.average(x_max)

    T_lst=[t_max[i+1]-t_max[i] for i in range(len(t_max)-1)]
    T=np.average(T_lst)

    lmbd_lst=[x_max[i+1]-x_max[i] for i in range(len(x_max)-1)]
    lmbd=np.average(lmbd_lst)

    return A, T, lmbd

t = np.linspace(ti,tf,n+1)
values = oscSimpFA_1D(x0,v0,k,m,t,b,F0,w_f,n,dt)
x = values[0]

reg_est = 250
values1 = amp_per_comp(x, t, n, reg_est)
A = values1[0]

print("Amplitude: ", A)
plt.plot(t, x)
plt.grid()
plt.show()