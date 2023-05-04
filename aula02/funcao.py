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
    r2 = (N*sumxy-(sumx*sumy))**2 / ((N*sumx2-sumx**2)*(N*sumy2-sumy**2))

    dm = abs(m) * math.sqrt((1/r2-1)/(N-2))
    db = dm * math.sqrt(sumx2/N)

    return m,b,r2,dm,db