import numpy as np

import matplotlib.pyplot as plt

x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
y = np.array([2, 4, 5, 4, 5, 6, 3, 8, 5, 4])

z = np.polyfit(x, y, len(x)-1)

p = np.poly1d(z)


xp = np.linspace(min(x), max(x), 1000)

plt.plot(x, y, '.', xp, p(xp), '-')
plt.show()
