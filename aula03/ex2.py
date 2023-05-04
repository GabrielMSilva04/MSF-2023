import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

x, y, t = sp.symbols('x, y, t')

vt = 6.8 #m/s
g = 9.8 #m/s**2

y = (vt**2/g)*sp.log(sp.cosh(g*t/vt))
v = sp.diff(y)
a = sp.diff(v)

tempo1 = sp.nsolve(y-20, t, 0) #y - 20= 0 -> y= 20

tempo2 = sp.nsolve(0.5*(g*(t**2))-20, t, 0) # 1/2 * at**2

y = sp.lambdify(t, y, 'numpy')
v = sp.lambdify(t, v, 'numpy')
a = sp.lambdify(t, a, 'numpy')

print(v(float(tempo1)))

plt.plot(np.linspace(0,4,1000),y(np.linspace(0,4,1000))) #1000 pontos em 0<x<4
plt.plot(np.linspace(0,4,1000),v(np.linspace(0,4,1000)))
plt.plot(np.linspace(0,4,1000),a(np.linspace(0,4,1000)))


print("com resistência do ar: "+str(tempo1))
print("velocidade final: "+ str(v(float(tempo1))))
print("aceleração final: "+ str(a(float(tempo1))))
print("sem resistência do ar: "+str(tempo2))
print("velocidade final: "+ str(g*tempo2)) #v = a*t
print("aceleração final: "+ str(g))

plt.show()