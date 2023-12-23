#113786

import numpy as np
from matplotlib import pyplot as plt

m = 1 #kg
xeq = 0
k = 1 #N/m
alpha = 0.05 #N/m^2
x = np.arange(-8, 4, 0.1)
Ep = 0.5*k*(x**2)+alpha*(x**3)

plt.plot(x,Ep)
plt.grid()

plt.title('1a) Diagrama de Energia Potencial')
plt.xlabel('x(m)')
plt.ylabel('Energia Potencial(J)')

plt.show()
