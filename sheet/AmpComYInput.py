import numpy as np


def intlaginvv(yinp, xm1, xm2, xm3, ym1, ym2, ym3):
    # interpolação inversa usando o polinómio de Lagrange
    # Dados (input): yinp, (x0,y0), (x1,y1), (x2,y2)
    # Resultados (output): xout, yout
    xab = xm1 - xm2
    xac = xm1 - xm3
    xbc = xm2 - xm3
    a = ym1 / (xab * xac)
    b = -ym2 / (xab * xbc)
    c = ym3 / (xac * xbc)
    am = a + b + c
    bm = a * (xm2 + xm3) + b * (xm1 + xm3) + c * (xm1 + xm2)
    cm = yinp + a * xm2 * xm3 + b * xm1 * xm3 + c * xm1 * xm2
    xout = (bm + np.sqrt(bm * bm - 4 * am * cm)) / (2 * am)

    if xm3 > xm1 and (xout < xm1 or xout > xm3):
        xout = (bm - np.sqrt(bm * bm - 4 * am * cm)) / (2 * am)
    if xm1 > xm3 and (xout < xm3 or xout > xm1):
        xout = (bm - np.sqrt(bm * bm - 4 * am * cm)) / (2 * am)

    xta = xout - xm1
    xtb = xout - xm2
    xtc = xout - xm3
    yout = a * xtb * xtc + b * xta * xtc + c * xta * xtb
    return xout, yout
