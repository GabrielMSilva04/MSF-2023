import math
import numpy as np
import matplotlib.pyplot as plt

def funcao(x,y):
    N = len(x)
    sumx = sum(xi for xi in x)
    sumy = sum(yi for yi in y)
    sumxy = sum(x[i]*y[i] for i in range(len(x)))
    sumx2 = sum(xi**2 for xi in x)
    sumy2 = sum(yi**2 for yi in y)


    m = (N*(sumxy) - (sumx * sumy))/(N*(sumx2)-sumx**2)
    b = (sumx2*sumy - sumx*sumxy)/((N*sumx2)-sumx**2)
    r2 = (N*sumxy-(sumx*sumy))**2 / ((N*sumx2-sumx**2)*(N*sumy2-sumy**2)) #quanto mais perto de 1, melhor o ajuste é

    delta_m = abs(m) * math.sqrt((1/r2-1)/(N-2))
    delta_b = delta_m * math.sqrt(sumx2/N)

    return m,b,r2,delta_m,delta_b

x = np.array([222, 207.5, 194, 171.5, 153, 133, 113, 92])
y = np.array([3, 2.2, 2, 1.8, 1.6, 1.4, 1.2, 1])

f = funcao(x, y)
m, b, r2, delta_m, delta_b = f

print("m =", m)
print("b =", b)
print("r2 =", r2)
print("erro de m =", delta_m)
print("erro de b =", delta_b)


plt.scatter(x, y, color='blue')
plt.plot(x, y, 'o')
plt.plot(x, m*x + b, color='green', label='Reta de ajuste')

plt.show()