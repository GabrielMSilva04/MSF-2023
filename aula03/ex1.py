import numpy as np
import matplotlib.pyplot as plt

vA = 70 / (3600/1000)

def f(time):
    lstA = []
    lstB = []
    x = []
    dA = 0
    dB = 0
    for i in range(0,time+1):
        dA = vA * i
        lstA.append(dA)
        dB = i**2
        lstB.append(dB)
        x.append(i)
    return lstA,lstB,x

ftuple = f(30) #tempo

yA = ftuple[0]
yB = ftuple[1]
x = ftuple[2]

plt.plot(x, yA)
plt.plot(x, yB)

plt.xlabel('t (s)')
plt.ylabel('x (m)')

#ponto
#vt=t**2 -> v=t
plt.plot(vA, vA*vA, 'o')

plt.show()

print("instante: "+str(vA))
print("distância percorrida: "+str(vA*vA))

