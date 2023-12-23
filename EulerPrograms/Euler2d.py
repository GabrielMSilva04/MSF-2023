import numpy as np
import matplotlib.pyplot as plt
import math

def graph(x,y,title = "Gráfico", xaxis = "Eixo X", yaxis = "Eixo Y"):
    plt.plot(x,y)

    plt.xlabel(xaxis)
    plt.ylabel(yaxis)

    plt.grid()

    plt.title(title)

# máximo pelo polinómio de Lagrange
def maximo(xm1,xm2,xm3,ym1,ym2,ym3):  
    xab=xm1-xm2
    xac=xm1-xm3
    xbc=xm2-xm3

    a=ym1/(xab*xac)
    b=-ym2/(xab*xbc)
    c=ym3/(xac*xbc)

    xmla=(b+c)*xm1+(a+c)*xm2+(a+b)*xm3
    xmax=0.5*xmla/(a+b+c)

    xta=xmax-xm1
    xtb=xmax-xm2
    xtc=xmax-xm3

    ymax=a*xtb*xtc+b*xta*xtc+c*xta*xtb
    return xmax, ymax

# raiz pelo polinómio de Lagrange
def zerosv(xm1,xm2,xm3,ym1,ym2,ym3):  
    xab=xm1-xm2
    xac=xm1-xm3
    xbc=xm2-xm3

    a=ym1/(xab*xac)
    b=-ym2/(xab*xbc)
    c=ym3/(xac*xbc)

    am=a+b+c
    bm=a*(xm2+xm3)+b*(xm1+xm3)+c*(xm1+xm2)
    cm=a*xm2*xm3+b*xm1*xm3+c*xm1*xm2

    xzero=(bm+np.sqrt(bm*bm-4*am*cm))/(2*am)
    if xm3 > xm1 and (xzero < xm1 or xzero > xm3): 
        xzero=(bm-np.sqrt(bm*bm-4*am*cm))/(2*am)


    if xm1 > xm3 and (xzero < xm3 or xzero > xm1):
        xzero=(bm-np.sqrt(bm*bm-4*am*cm))/(2*am)

    xta=xzero-xm1
    xtb=xzero-xm2
    xtc=xzero-xm3
    yzero=a*xtb*xtc+b*xta*xtc+c*xta*xtb
    return xzero, yzero

x0 = 
y0 = 

alfa = *math.pi/180
v0 = *1000/3600##velocidade inicial (m/s)
vx0 = v0*math.cos(alfa)
vy0 = v0*math.sin(alfa)

g = 9.80 

## Passo temporal para o Método de Euler
dt = 0.001 ##passo temporal
t0 = 0.00 ##tempo inicial
tf = 5.00 ##tempo final
N = int((tf-t0)/dt) ##Nº de passos a realizar
print("Número de passos: ", N)

## Criação de arrays de velocidade, posição e tempo com N lugares, para os N passos,
#  preenchidos por 0's. Caso a velocidade ou posição inicial seja diferente de 0, temos que fazer a correção.

## NOTA: Vx é constante s/ resistência do ar (igual a Vx0)
t = np.zeros(N+1)

vy = np.zeros(N+1)
vx = vx0

x = np.zeros(N+1)
y = np.zeros(N+1)

vy[0] = vy0

## Ciclo que implementa o Método de Euler, e atualiza o array tempo para podemos fazer um plot
for i in range(N):
    t[i+1] = t[i] + dt

    vy[i+1] = vy[i] - g*dt

    x[i+1] = x[i] + vx*dt
    y[i+1] = y[i] + vy[i]*dt


graph(x,y,"Altura x Distância", "Distância (m)", "Altura (m)")
plt.show()

## Cálculo da altura máxima e do alcance com recurso aos polinómios de Lagrange

for i in range(N):
    if y[i-1] < y[i] and  y[i+1] < y[i]:
        print('Success')
        maxx, maxy = maximo(x[i-1], x[i], x[i+1], y[i-1], y[i], y[i+1])
        print("A altura máxima atingida foi {:.2f}m após {:.2f}s.".format(maxy, t[i]))
        break

for i in range(N):
    if y[i]*y[i-1] < 0:
        print('Success')
        xzero, yzero = zerosv(x[i-1], x[i], x[i+1], y[i-1], y[i], y[i+1])
        print("O alcance da bola foi {:.2f}m após {:.2f}s.".format(xzero, t[i]))
        break