import numpy as np
import matplotlib.pyplot as plt
from funcao import funcao

f = funcao([222, 207.5, 194, 171.5, 153, 133, 113, 92], [2.3, 2.2, 2, 1.8, 1.6, 1.4, 1.2, 1])



x = [222, 207.5, 194, 171.5, 153, 133, 113, 92]
y = [2.3, 2.2, 2, 1.8, 1.6, 1.4, 1.2, 1]
plt.plot(x, y, 'o')

x2 = x
y2 = [f[0]*xi+f[1] for xi in x]
plt.plot(x2, y2)

plt.show()

print("L = 165.0cm , X = " + str(f[0]*165+f[1]))