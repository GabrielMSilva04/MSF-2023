import matplotlib.pyplot as plt
import numpy as np

t = np.array([0.5+i for i in range(0,11)]) #0.5 , 1.5 ... 9.5 , 10.5
s = np.array([0.1, 1.4, 1.7, 6.5, 7.7, 10.4, 19.5, 26.1, 26.5, 45.9, 52.5])

x = t
y = s

npontos = x.size

xy = x*y  # element by element product
x2 = x*x
y2 = y*y

sx = x.sum()
sy = y.sum()
sxy = xy.sum()
sxx = x2.sum()
syy = y2.sum()


n = npontos
rn = n*sxy-sx*sy
rd = (n*sxx-sx**2)*(n*syy-sy**2)
r2 = rn**2/rd
r = np.sqrt(r2)
m = (n*sxy-sx*sy)/(n*sxx-sx**2)
dm = abs(m)*np.sqrt((1/r**2-1)/(n-2))
bn = sxx*sy-sx*sxy
bd = n*sxx-sx**2
b = bn/bd
db = dm*np.sqrt(sxx/n)

print('m +/-dm= ', m, "+/-", dm)
print('b +/-db= ', b, "+/-", db)
print('r2= ', r2)

reta = np.polyfit(x, y, 1)

plt.scatter(x, y, color='orange')

plt.plot(t, reta[0]*t + reta[1], color='green', label='Reta de ajuste')
plt.title("Gráfico")
plt.legend()
plt.xlabel("Tempo (s)")
plt.ylabel("Posição (m)")
plt.grid()
plt.show()
