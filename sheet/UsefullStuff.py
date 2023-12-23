import numpy as np
from matplotlib import pyplot as plt
import sympy as sym
from mpl_toolkits import mplot3d
from numpy import linalg as LA


# funções:

# auxiliares:


def prodExt(a, b):
    # calcula o produto externo de dois vetores tridimensionais
    return (a[1] * b[2] - b[1] * a[2], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def aproxRet(f, n, dx):
    # realiza a aproximação regular de uma função
    return dx * np.sum(f[0:n])


def aproxTrap(f, n, dx):
    # realiza a aproximação trapezoidal de uma função
    return dx * (0.5 * (f[0] + f[n]) + np.sum(f[1:n]))


def maxminv(xm1, xm2, xm3, ym1, ym2, ym3):
    xab = xm1 - xm2
    xac = xm1 - xm3
    xbc = xm2 - xm3
    a = ym1 / (xab * xac)
    b = -ym2 / (xab * xbc)
    c = ym3 / (xac * xbc)
    xmla = (b + c) * xm1 + (a + c) * xm2 + (a + b) * xm3
    xm = 0.5 * xmla / (a + b + c)
    xta = xm - xm1
    xtb = xm - xm2
    xtc = xm - xm3
    ymax = a * xtb * xtc + b * xta * xtc + c * xta * xtb
    return xm, ymax


def lim(t, x, n):
    indMax = []
    for i in range(n):
        if (x[i - 1] < x[i] and x[i] > x[i + 1] and t[i] > 0):
            indMax.append(i)

    indMin = []
    for i in range(n):
        if (x[i - 1] > x[i] and x[i] < x[i + 1] and t[i] > 0):
            indMin.append(i)

    tmax = np.zeros(len(indMax))
    xmax = np.zeros(len(indMax))
    tmin = np.zeros(len(indMin))
    xmin = np.zeros(len(indMin))
    c = 0
    for i in indMax:
        tmax[c], xmax[c] = maxminv(t[i - 1], t[i], t[i + 1], x[i - 1], x[i], x[i + 1])
        c += 1

    j = 0
    for i in indMin:
        tmin[j], xmin[j] = maxminv(t[i - 1], t[i], t[i + 1], x[i - 1], x[i], x[i + 1])
        j += 1

    amplitude = np.mean(xmax)
    print("Amplitude (Máximo):", amplitude)
    minimo = np.mean(xmin)
    print("Mínimo:", minimo)
    periodo = tmax[1] - tmax[0]
    print("Período:", periodo)
    freq = 1 / periodo
    print("Frequência:", freq)


def euler(t, x0, v0, dt, n):
    x = np.empty(n + 1)
    x[0] = x0
    v = np.empty(n + 1)
    v[0] = v0
    for i in range(n):
        a = acel(t[i], x[i], v[i])
        v[i + 1] = v[i] + a * dt
        x[i + 1] = x[i] + v[i] * dt
    return x, v


def fourier(x0, v0, m, k, b, F0, Wf, A, t, n, ti_estac, dt):
    # faz a preparação necessária ao cálculo dos coeficientes de Fourier e apresenta-os num gŕafico de barras
    # atenção: alterar a equação da aceleração se necessário
    x = np.empty(n + 1)
    v = np.empty(n + 1)
    a = np.empty(n + 1)
    E = np.empty(n + 1)

    v[0] = v0
    x[0] = x0
    cntMax = 0
    ind = np.transpose([0 for i in range(1000)])
    afo = np.zeros(15)
    bfo = np.zeros(15)

    for i in range(n):
        a[i] = -(k / m) * x[i] - (b / m) * v[i] + (F0 / m) * np.cos(Wf * t[i])
        v[i + 1] = v[i] + a[i] * dt
        x[i + 1] = x[i] + v[i + 1] * dt
        E[i] = 0.5 * m * v[i] ** 2 + 0.5 * k * x[i] ** 2
        if t[i] > ti_estac and x[i - 1] < x[i] and x[i + 1] < x[i]:
            cntMax += 1
            ind[cntMax] = int(i)

    t0 = ind[cntMax - 1]
    t1 = ind[cntMax]
    for i in range(15):
        af, bf = fourier_calc(t, x, t0, t1, dt, i)
        afo[i] = af
        bfo[i] = bf

    ii = np.linspace(0, 14, 15)
    plt.figure()
    plt.ylabel('| a_n |')
    plt.xlabel('n')
    plt.bar(ii, np.abs(afo))
    plt.grid()
    plt.show()

    ii = np.linspace(0, 14, 15)
    plt.figure()
    plt.ylabel('| b_n |')
    plt.xlabel('n')
    plt.bar(ii, np.abs(bfo))
    plt.grid()
    plt.show()


def fourier_calc(t, x, t0, t1, dt, nf):
    # faz os cálculos dos coeficientes de Fourier. Parte de 'fourier'
    T = t[t1] - t[t0]
    ome = 2 * np.pi / T

    s1 = x[t0] * np.cos(nf * ome * t[t0])
    s2 = x[t1] * np.cos(nf * ome * t[t1])
    st = x[t0 + 1:t1] * np.cos(nf * ome * t[t0 + 1:t1])
    soma = np.sum(st)

    q1 = x[t0] * np.sin(nf * ome * t[t0])
    q2 = x[t1] * np.sin(nf * ome * t[t1])
    qt = x[t0 + 1:t1] * np.sin(nf * ome * t[t0 + 1:t1])
    somq = np.sum(qt)

    intega = ((s1 + s2) / 2 + soma) * dt
    af = 2 / T * intega
    integq = ((q1 + q2) / 2 + somq) * dt
    bf = 2 / T * integq
    return af, bf


# 1 dimensão:


def regLin(x, y):
    # calcula a regressão linear de uma função
    x2 = x ** 2
    y2 = y ** 2
    xy = x * y
    n = x.size

    sx = x.sum()
    sy = y.sum()
    sxy = xy.sum()
    sx2 = x2.sum()
    sy2 = y2.sum()

    m = (n * sxy - sx * sy) / (n * sx2 - (sx ** 2))
    b = (sx2 * sy - sx * sxy) / (n * sx2 - (sx ** 2))

    r2n = n * sxy - sx * sy
    r2d1 = n * sx2 - (sx ** 2)
    r2d2 = n * sy2 - (sy ** 2)
    r2 = (r2n ** 2) / (r2d1 * r2d2)

    varM = abs(m) * np.sqrt(((1 / r2) - 1) / (n - 2))
    varB = varM * np.sqrt(sx2 / n)
    return m, b, r2, varM, varB


def getR_1D(r0, v0, dt, n, a):
    # obtém a posição e velocidade a uma dimensão com aceleração constante
    v = np.empty(n + 1)
    v[0] = v0
    r = np.empty(n + 1)
    r[0] = r0
    for i in range(n):
        v[i + 1] = v[i] + a * dt
        r[i + 1] = r[i] + v[i] * dt
    return r, v


def planoInclinado_simp_1D(x0, v0, n, dt, m, pot=0, teta=0):
    # calcula posição, velocidade e aceleração para um corppo a (tentar) subir sem forças resistentes
    x = np.empty(n + 1)
    vx = np.empty(n + 1)
    ax = np.empty(n + 1)
    g = 9.8

    x[0] = x0
    vx[0] = v0
    ax[0] = 0

    for i in range(n):
        vv = np.abs(vx[i])
        ax[i] = -g * np.sin(teta) + pot / (m * vx[i])
        vx[i + 1] = vx[i] + ax[i] * dt
        x[i + 1] = x[i] + vx[i] * dt
    return x, vx, ax


def planoInclinado_atr_1D(x0, v0, n, dt, c_atr, A, m, pot, teta=0):
    # calcula posição, velocidade e aceleração para um corppo a (tentar) subir com atrito
    x = np.empty(n + 1)
    vx = np.empty(n + 1)
    ax = np.empty(n + 1)
    g = 9.8

    x[0] = x0
    vx[0] = v0
    ax[0] = 0

    for i in range(n):
        vv = np.abs(vx[i])
        f_cic = pot / vx[i]
        ax[i] = -g * np.sin(teta) - c_atr * g * np.cos(teta) + pot / (m * vx[i])
        vx[i + 1] = vx[i] + ax[i] * dt
        x[i + 1] = x[i] + vx[i] * dt
    return x, vx, ax


def planoInclinado_res_1D(x0, v0, n, dt, cres, u, A, m, P, ang=0):
    # calcula posição, velocidade e aceleração para um corppo a (tentar) subir com atrito e resistência do ar
    x = np.empty(n + 1)
    vx = np.empty(n + 1)
    ax = np.empty(n + 1)

    p_ar = 1.225
    g = 9.8

    x[0] = x0
    vx[0] = v0
    ax[0] = 0

    for i in range(n):
        vv = np.abs(vx[i])
        ax[i] = -g * np.sin(ang) - u * g * np.cos(ang) - (0.5 * cres * A * p_ar * vx[i] * vv) / m + P / (m * vx[i])
        vx[i + 1] = vx[i] + ax[i] * dt
        x[i + 1] = x[i] + vx[i] * dt
    return x, vx, ax


def amp_per_comp(x, t, n):
    # amplitde, periodo e comprimeno de onda
    ind_max = [i for i in range(1, n - 1) if x[i - 1] <= x[i] >= x[i + 1]]
    x_max = [x[i] for i in ind_max]
    t_max = [t[i] for i in ind_max]
    A = np.average(x_max)

    T_lst = [t_max[i + 1] - t_max[i] for i in range(len(t_max) - 1)]
    T = np.average(T_lst)

    lmbd_lst = [x_max[i + 1] - x_max[i] for i in range(len(x_max) - 1)]
    lmbd = np.average(lmbd_lst)

    return A, T, lmbd


def amp_per_comp(x, t, n, reg_est):
    # amplitde, periodo e comprimeno de onda com regime estacionario
    ind_max = [i for i in range(1, n - 1) if x[i - 1] <= x[i] >= x[i + 1] if t[i] > reg_est]
    x_max = [x[i] for i in ind_max]
    t_max = [t[i] for i in ind_max]
    A = np.average(x_max)

    T_lst = [t_max[i + 1] - t_max[i] for i in range(len(t_max) - 1)]
    T = np.average(T_lst)

    lmbd_lst = [x_max[i + 1] - x_max[i] for i in range(len(x_max) - 1)]
    lmbd = np.average(lmbd_lst)

    return A, T, lmbd


def amp_per_comp(x, t, n):
    # x minimo e máximo, periodo e comprimento de onda
    ind_max = [i for i in range(1, n - 1) if x[i - 1] <= x[i] >= x[i + 1]]
    ind_min = [i for i in range(1, n - 1) if x[i - 1] >= x[i] <= x[i + 1]]
    x_max = [x[i] for i in ind_max]
    x_min = [x[i] for i in ind_min]
    t_max = [t[i] for i in ind_max]
    x_max1 = np.average(x_max)
    x_min1 = np.average(x_min)

    T_lst = [t_max[i + 1] - t_max[i] for i in range(len(t_max) - 1)]
    T = np.average(T_lst)

    lmbd_lst = [x_max[i + 1] - x_max[i] for i in range(len(x_max) - 1)]
    lmbd = np.average(lmbd_lst)

    return x_max1, x_min1, T, lmbd


def amp_per_comp(x, t, n, reg_est):
    # x minimo e máximo, periodo e comprimento de onda com regime estacionario
    ind_max = [i for i in range(1, n - 1) if x[i - 1] <= x[i] >= x[i + 1] if t[i] > reg_est]
    ind_min = [i for i in range(1, n - 1) if x[i - 1] >= x[i] <= x[i + 1] if t[i] > reg_est]
    x_max = [x[i] for i in ind_max]
    x_min = [x[i] for i in ind_min]
    t_max = [t[i] for i in ind_max]
    x_max1 = np.average(x_max)
    x_min1 = np.average(x_min)

    T_lst = [t_max[i + 1] - t_max[i] for i in range(len(t_max) - 1)]
    T = np.average(T_lst)

    lmbd_lst = [x_max[i + 1] - x_max[i] for i in range(len(x_max) - 1)]
    lmbd = np.average(lmbd_lst)

    return x_max1, x_min1, T, lmbd


def oscHarmSimp_1D(x0, v0, k, m, n, dt):
    # oscilador Harmónico simples
    x = np.empty(n + 1)
    v = np.empty(n + 1)
    a = np.empty(n + 1)
    x[0] = x0
    v[0] = v0
    for i in range(n):
        a[i] = -k / m * x[i]
        v[i + 1] = v[i] + a[i] * dt
        x[i + 1] = x[i] + v[i + 1] * dt
    return x, v, a


def oscSimpFA_1D(x0, v0, k, m, t, b, F0, w_f, n, dt):
    # oscilador Simples sujeito forçado e amortecido
    x = np.empty(n + 1)
    v = np.empty(n + 1)
    a = np.empty(n + 1)
    x[0] = x0
    v[0] = v0
    for i in range(n):
        a[i] = -k / m * x[i] + (-b * v[i] + F0 * np.cos(w_f * t[i])) / m
        v[i + 1] = v[i] + a[i] * dt
        x[i + 1] = x[i] + v[i + 1] * dt
    return x, v, a


global g, vt


def rk4_1D(t, x0, v0, dt, n):
    # método de Runge-Kutta de 4ª ordem que retorna a posição e a velocidade, requer que se defina à parte a função acel para a aceleração
    x = np.empty(n + 1)
    x[0] = x0
    v = np.empty(n + 1)
    v[0] = v0
    for i in range(n):
        ax1 = acel(t[i], x[i], v[i])
        c1v = ax1 * dt
        c1x = v[i] * dt

        ax2 = acel(t[i] + dt / 2, x[i] + c1x / 2, v[i] + c1v / 2)
        c2v = ax2 * dt
        c2x = (v[i] + c1v / 2) * dt

        ax3 = acel(t[i] + dt / 2, x[i] + c2x / 2, v[i] + c2v / 2)
        c3v = ax3 * dt
        c3x = (v[i] + c2v / 2) * dt

        ax4 = acel(t[i] + dt / 2, x[i] + c2x / 2, v[i] + c2v / 2)
        c4v = ax4 * dt
        c4x = (v[i] + c3v) * dt

        x[i + 1] = x[i] + (c1x + 2 * c2x + 2 * c3x + c4x) / 6
        v[i + 1] = v[i] + (c1v + 2 * c2v + 2 * c3v + c4v) / 6
    return x, v


def oscAcop_EulerCromer_1D(x0, v0, x_eq, k, m, n, dt):
    # oscilador acoplado harmónico de dois corpos por Euler-Cromer
    mA = m[0]
    mB = m[1]

    xA = np.empty(n + 1)
    xA[0] = x0[0]
    xB = np.empty(n + 1)
    xB[0] = x0[1]

    xA_eq = x_eq[0]
    xB_eq = x_eq[1]

    vA = np.empty(n + 1)
    vA[0] = v0[0]
    vB = np.empty(n + 1)
    vB[0] = v0[1]

    aA = np.empty(n + 1)
    aB = np.empty(n + 1)

    k_ext = k[0]
    k_int = k[1]

    for i in range(n):
        aA[i] = (-k_ext * (xA[i] - xA_eq) - k_int * ((xA[i] - xA_eq) - (xB[i] - xB_eq))) / mA
        aB[i] = (-k_ext * (xB[i] - xB_eq) + k_int * ((xA[i] - xA_eq) - (xB[i] - xB_eq))) / mB

        vA[i + 1] = vA[i] + aA[i] * dt
        vB[i + 1] = vB[i] + aB[i] * dt

        xA[i + 1] = xA[i] + vA[i + 1] * dt
        xB[i + 1] = xB[i] + vB[i + 1] * dt
    return (xA, xB), (vA, vB), (aA, aB)


def oscAcop_ModosNormais_1D(x_eq, A, w, t, phi, n, dt):
    # oscilador acoplado harmónico de dois corpos por sobreposição de Modos Normais
    xA = np.empty(n + 1)
    xB = np.empty(n + 1)

    xA_eq = x_eq[0]
    xB_eq = x_eq[1]

    A1 = A[0]
    A2 = A[1]
    w1 = w[0]
    w2 = w[1]
    phi1 = phi[0]
    phi2 = phi[1]

    xA = xA_eq + A1 * np.cos(w1 * t + phi1) + A2 * np.cos(w2 * t + phi2)
    xB = xB_eq + A1 * np.cos(w1 * t + phi1) - A2 * np.cos(w2 * t + phi2)
    return xA, xB


def oscAcop_EulerCromer_FA_1D(x0, v0, x_eq, k, m, t, b, F0, w_f, n, dt):
    # oscilador acoplado forçado e amortecido de dois corpos por Euler-Cromer
    xA = np.empty(n + 1)
    xA[0] = x0[0]
    xB = np.empty(n + 1)
    xB[0] = x0[1]

    mA = m[0]
    mB = m[1]
    # se não forem aplicads forças num corpo, F0 e w_f são zero nesse corpo
    F0A = F0[0]
    F0B = F0[1]
    w_f_A = w_f[0]
    w_f_B = w_f[1]

    xA_eq = x_eq[0]
    xB_eq = x_eq[1]

    vA = np.empty(n + 1)
    vA[0] = v0[0]
    vB = np.empty(n + 1)
    vB[0] = v0[1]

    aA = np.empty(n + 1)
    aB = np.empty(n + 1)

    k_ext = k[0]
    k_int = k[1]

    for i in range(n):
        aA[i] = (-k_ext * (xA[i] - xA_eq) - k_int * ((xA[i] - xA_eq) - (xB[i] - xB_eq)) - b * vA[i] + F0A * np.cos(
            w_f_A * t[i])) / mA
        aB[i] = (-k_ext * (xB[i] - xB_eq) + k_int * ((xA[i] - xA_eq) - (xB[i] - xB_eq)) - b * vB[i] + F0B * np.cos(
            w_f_B * t[i])) / mB

        vA[i + 1] = vA[i] + aA[i] * dt
        vB[i + 1] = vB[i] + aB[i] * dt

        xA[i + 1] = xA[i] + vA[i + 1] * dt
        xB[i + 1] = xB[i] + vB[i + 1] * dt
    return (xA, xB), (vA, vB), (aA, aB)


# as duas seguintes são mais complicados, não aconselho ver a menos que se perceba da coisa
def oscAcop_EulCrom_Mult_Simp_1D(x0, v0, x_eq, k, m, numCorp, n, dt):
    x = []
    v = []
    a = []
    # k são as constantes da molas, da esquerda para a direita

    # criar posição .velocidade e aceleração para numCorp objetos
    for i in range(numCorp):
        x[i] = np.empty(n + 1)
        v[i] = np.empty(n + 1)
        a[i] = np.empty(n + 1)
        # atribuir valores iniciais:
        x[i] = x0[i]
        v[i] = v0[i]

    for i in range(n):
        # aceleração dos dois corpos conectados à parede:
        a[0][i] = (-k[0] * (x[0][i] - x_eq[0]) - k[1] * ((x[0][i] - x_eq[0]) - (x[1][i] - x[1][i]))) / m[0]
        a[-1][i] = (-k[-1] * (x[-1][i] - x_eq[-1]) - k[-2]((x[-1][i] - x_eq[-1]) - (x[-2][i] - x_eq[-2]))) / m[-1]
        # aceleração dos corpos no meio:
        for iC in range(1, numCorp - 1):
            a[iC][i] = (-k[iC] * ((x[iC] - x_eq[iC]) - (x[iC - 1] - x_eq[iC - 1])) - k[iC + 1] * (
                        (x[iC] - x_eq[iC]) - (x[iC + 1] - x_eq[iC + 1]))) / m[iC]

        # velocidade e posição:
        for i2 in range(numCorp):
            v[i2][i + 1] = v[i2][i] + a[i2][i] * dt
            x[i2][i + 1] = x[i2][i] + v[i2][i + 1] * dt
        return x, v, a


def oscAcop_EulCrom_Mult_FA_1D(x0, v0, x_eq, k, m, numCorp, t, b, F0, w_f, n, dt):
    x = []
    v = []
    a = []
    # k são as constantes da molas, da esquerda para a direita

    # criar posição, velocidade e aceleração para numCorp objetos
    for i in range(numCorp):
        x[i] = np.empty(n + 1)
        v[i] = np.empty(n + 1)
        a[i] = np.empty(n + 1)
        # atribuir valores iniciais:
        x[i] = x0[i]
        v[i] = v0[i]

    for i in range(n):
        # aceleração dos dois corpos conectados à parede:
        a[0][i] = (-k[0] * (x[0][i] - x_eq[0]) - k[1] * ((x[0][i] - x_eq[0]) - (x[1][i] - x[1][i])) + F0[0] * np.cos(
            w_f[0] * t[i])) / m[0]
        a[-1][i] = (-k[-1] * (x[-1][i] - x_eq[-1]) - k[-2]((x[-1][i] - x_eq[-1]) - (x[-2][i] - x_eq[-2])) + F0[
            -1] * np.cos(w_f[-1] * t[i])) / m[-1]
        # aceleração dos corpos no meio:
        for iC in range(1, numCorp - 1):
            a[iC][i] = (-k[iC] * ((x[iC] - x_eq[iC]) - (x[iC - 1] - x_eq[iC - 1])) - k[iC + 1] * (
                        (x[iC] - x_eq[iC]) - (x[iC + 1] - x_eq[iC + 1])) + F0[iC] * np.cos(w_f[iC] * t[i])) / m[iC]

        # velocidade e posição:
        for i2 in range(numCorp):
            v[i2][i + 1] = v[i2][i] + a[i2][i] * dt
            x[i2][i + 1] = x[i2][i] + v[i2][i + 1] * dt
        return x, v, a


# 2 dimensões:


def getR_2D(r0, v0, a0, dt, n):
    # corpo num plano com aceleração constante
    x = np.empty(n + 1)
    y = np.empty(n + 1)
    vx = np.empty(n + 1)
    vy = np.empty(n + 1)

    x[0] = r0[0]
    y[0] = r0[1]
    vx[0] = v0[0]
    vy[0] = v0[1]
    ax = a0[0]
    ay = a0[1]

    for i in range(n):
        vx[i + 1] = vx[i] + ax * dt
        vy[i + 1] = vy[i] + ay * dt
        x[i + 1] = x[i] + vx[i] * dt
        y[i + 1] = y[i] + vy[i] * dt
    return (x, y), (vx, vy)


def resAr_2D(r0, v0, n, dt, vt):
    # corpo num plano com resistência do ar
    x = np.empty(n + 1)
    y = np.empty(n + 1)
    vx = np.empty(n + 1)
    vy = np.empty(n + 1)
    ax = np.empty(n + 1)
    ay = np.empty(n + 1)

    g = 9.80
    x[0] = r0[0]
    y[0] = r0[1]
    vx[0] = v0[0]
    vy[0] = v0[1]
    ax[0] = 0
    ay[0] = -g
    dres = g / vt ** 2

    for i in range(n):
        vv = np.sqrt(vx[i] ** 2 + vy[i] ** 2)

        ax[i] = -dres * vv * vx[i]
        ay[i] = -g - dres * vv * vy[i]

        vx[i + 1] = vx[i] + ax[i] * dt
        vy[i + 1] = vy[i] + ay[i] * dt

        x[i + 1] = x[i] + vx[i] * dt
        y[i + 1] = y[i] + vy[i] * dt
    return (x, y), (vx, vy), (ax, ay)


# 3 dimensões:


def getR_3D(r0, v0, a0, dt, n):
    # corpo no espaço com aceleração constante
    x = np.empty(n + 1)
    y = np.empty(n + 1)
    z = np.empty(n + 1)
    vx = np.empty(n + 1)
    vy = np.empty(n + 1)
    vz = np.empty(n + 1)

    x[0] = r0[0]
    y[0] = r0[1]
    z[0] = r0[2]
    vx[0] = v0[0]
    vy[0] = v0[1]
    vz[0] = v0[2]
    ax = a0[0]
    ay = a0[1]
    az = a0[2]

    for i in range(n):
        vx[i + 1] = vx[i] + ax * dt
        vy[i + 1] = vy[i] + ay * dt
        x[i + 1] = x[i] + vx[i] * dt
        y[i + 1] = y[i] + vy[i] * dt
    return (x, y, z), (vx, vy, vz)


def airRes_3D(r0, v0, a0, n, dt, vt):
    # corpo no espaço com resistência do ar
    x = np.empty(n + 1)
    y = np.empty(n + 1)
    z = np.empty(n + 1)
    vx = np.empty(n + 1)
    vy = np.empty(n + 1)
    vz = np.empty(n + 1)
    ax = np.empty(n + 1)
    ay = np.empty(n + 1)
    az = np.empty(n + 1)

    x[0] = r0[0]
    y[0] = r0[1]
    z[0] = r0[2]
    vx[0] = v0[0]
    vy[0] = v0[1]
    vz[0] = v0[2]
    ax[0] = a0[0]
    ay[0] = a0[1]
    az[0] = a0[2]
    g = 9.80
    dres = g / vt ** 2

    for i in range(n):
        vv = np.sqrt(vx[i] ** 2 + vy[i] ** 2)

        ax[i] = a0[0] - dres * vv * vx[i]
        ay[i] = a0[1] - dres * vv * vy[i]
        az[i] = a0[2] - dres * vv * vz[i]

        vx[i + 1] = vx[i] + ax[i] * dt
        vy[i + 1] = vy[i] + ay[i] * dt
        vz[i + 1] = vz[i] + az[i] * dt

        x[i + 1] = x[i] + vx[i] * dt
        y[i + 1] = y[i] + vy[i] * dt
        z[i + 1] = z[i] + vz[i] * dt
    return (x, y, z), (vx, vy, vz), (ax, ay, az)


def magnus_3D(r0, v0, a0, rot, p_ar, r, n, dt, vt, m):
    # corpo no espaço com resistência do ar e rotação
    g = 9.80
    A = np.pi * r ** 2
    apr = 0.5 * p_ar * A * r

    x = np.empty(n + 1)
    y = np.empty(n + 1)
    z = np.empty(n + 1)

    vx = np.empty(n + 1)
    vy = np.empty(n + 1)
    vz = np.empty(n + 1)

    ax = np.empty(n + 1)
    ay = np.empty(n + 1)
    az = np.empty(n + 1)

    x[0] = r0[0]
    y[0] = r0[1]
    z[0] = r0[2]

    vx[0] = v0[0]
    vy[0] = v0[1]
    vz[0] = v0[2]

    ax[0] = a0[0]
    ay[0] = a0[1]
    az[0] = a0[2]

    dres = g / vt ** 2
    for i in range(n):
        vv = np.sqrt(vx[i] ** 2 + vy[i] ** 2 + vz[i] ** 2)
        rot_v = prodExt(rot, (vx[i], vy[i], vz[i]))

        mag_x = apr * rot_v[0] / m
        mag_y = apr * rot_v[1] / m
        mag_z = apr * rot_v[2] / m

        ax[i] = a0[0] - dres * vv * vx[i] + mag_x
        ay[i] = a0[1] - dres * vv * vy[i] + mag_y
        az[i] = a0[2] - dres * vv * vz[i] + mag_z

        vx[i + 1] = vx[i] + ax[i] * dt
        vy[i + 1] = vy[i] + ay[i] * dt
        vz[i + 1] = vz[i] + az[i] * dt

        x[i + 1] = x[i] + vx[i] * dt
        y[i + 1] = y[i] + vy[i] * dt
        z[i + 1] = z[i] + vz[i] * dt
    return (x, y, z), (vx, vy, vz), (ax, ay, az)
