import numpy as np
import math
import matplotlib.pyplot as plt

# parametros
v0 = 100 / 3.6
vy0 = v0 * math.sin(math.pi / 18)
vx0 = v0 * math.cos(math.pi / 18)

vterminal = 100 / 3.6
y0 = 0
x0 = 0

g = 9.8
t0 = 0


# inic
def euler(tf, deltat):  # deltat = tempo para cada intervalo na integração
    n = np.int_((tf - t0) / deltat)

    t = np.zeros(n + 1)  # array com tamanho igual ao número de divisões na integração
    x = np.zeros(n + 1)
    y = np.zeros(n + 1)
    vy = np.zeros(n + 1)
    vx = np.zeros(n + 1)
    ay = np.zeros(n + 1)
    ax = np.zeros(n + 1)

    xr = np.zeros(n + 1)
    yr = np.zeros(n + 1)
    vyr = np.zeros(n + 1)
    vxr = np.zeros(n + 1)
    ayr = np.zeros(n + 1)
    axr = np.zeros(n + 1)

    vy[0] = vy0
    vx[0] = vx0
    vyr[0] = vy0
    vxr[0] = vx0
    ayr[0] = -9.8
    t[0] = t0

    # euler method

    for i in range(n):
        # s/ resistência
        ax[i] = 0
        ay[i] = -g  # queda livre
        # (em geral pode ser qualquer função de x[i] e vx[i])

        y[i + 1] = y[i] + vy[i] * deltat + (1 / 2 * ay[i] * deltat ** 2)
        vy[i + 1] = vy[i] + ay[i] * deltat  # atualizar velocidade sabendo aceleração

        x[i + 1] = x[i] + vx[i] * deltat + (1 / 2 * ax[i] * deltat ** 2)
        vx[i + 1] = vx[i] + ax[i] * deltat

        # c/ resistência
        axr[i + 1] = ax[i] - ((g / vterminal ** 2) * math.sqrt(vxr[i] ** 2 + vyr[i] ** 2) * vxr[i])
        ayr[i + 1] = ay[i] - ((g / vterminal ** 2) * math.sqrt(vyr[i] ** 2 + vyr[i] ** 2) * vyr[i])

        yr[i + 1] = yr[i] + vyr[i] * deltat + (1 / 2 * ayr[i] * deltat ** 2)
        vyr[i + 1] = vyr[i] + ayr[i] * deltat

        xr[i + 1] = xr[i] + vxr[i] * deltat + (1 / 2 * axr[i] * deltat ** 2)
        vxr[i + 1] = vxr[i] + axr[i] * deltat

        t[i + 1] = t[i] + deltat
    return t, x, y, xr, yr


x = np.array(euler(10, 0.01)[1])
y = np.array(euler(10, 0.01)[2])
xr = np.array(euler(10, 0.01)[3])
yr = np.array(euler(10, 0.01)[4])

plt.plot(x, y, label='s/ resistencia')
plt.plot(xr, yr, label='c/ resistencia')

plt.ylim([0, 1.5])
plt.xlim([0, 30])
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.legend()
plt.show()
