import numpy as np
import matplotlib.pyplot as plt
import math

def graph(x,y,title = "Gráfico", xaxis = "Eixo X", yaxis = "Eixo Y"):
    plt.plot(x,y,'o')

    plt.xlabel(xaxis)
    plt.ylabel(yaxis)

    plt.grid()
    
    #plt.xlim(0, 25)
    #plt.ylim(0, 2)

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

xx = 100 ##velocidade em km/h

x0 = 0 ##posição inicial x
y0 = 0 ##posição inicial y

alfagraus = 10 ##angulo de lançamento (graus) INPUT
alfa = alfagraus * math.pi/180 ##angulo de lançamento (radianos)

v0 = xx*1000/3600 ##velocidade inicial (m/s)
vx0 = v0*math.cos(alfa)
vy0 = v0*math.sin(alfa)

g = 9.80
mass = 1 ##massa do objeto kg

## Passo temporal para o Método de Euler
dt = 0.001 ##passo temporal
t0 = 0.00 ##tempo inicial
tf = 5.00 ##tempo final
N = int((tf-t0)/dt) ##Nº de passos a realizar
print("Número de passos: ", N)

## Com resistência do ar
vt = xx*1000/3600
D = g/vt**2  ##Cálculo de D para o efeito da Resistência do Ar


## Criação de arrays de velocidade, posição e tempo com N lugares, para os N passos,
#  preenchidos por 0's. Caso a velocidade ou posição inicial seja diferente de 0, temos que fazer a correção.

## NOTA: Vx é constante s/ resistência do ar (igual a Vx0)
t = np.zeros(N+1)

ax = np.zeros(N+1)
ay = np.zeros(N+1)

vx = np.zeros(N+1)
vy = np.zeros(N+1)

x = np.zeros(N+1)
y = np.zeros(N+1)

## Cálculo do trabalho da Rar
fx = np.zeros(N+1)
fy = np.zeros(N+1)
ayrar = np.zeros(N+1)

vx[0] = vx0
vy[0] = vy0

x[0] = x0
y[0] = y0


## Ciclo que implementa o Método de Euler, e atualiza o array tempo para podemos fazer um plot
for i in range(N):
    t[i+1] = t[i] + dt

    vv = math.sqrt(vx[i]**2 + vy[i]**2)

    ay[i] = -( g + (D*vy[i]*vv))
    ayrar[i] = -(D*vy[i]*vv)
    ax[i] = -(D*vx[i]*vv)

    vy[i+1] = vy[i] + ay[i]*dt
    vx[i+1] = vx[i] + ax[i]*dt

    x[i+1] = x[i] + vx[i]*dt
    y[i+1] = y[i] + vy[i]*dt

    fx[i] = mass*ax[i]*vx[i] 
    fy[i] = mass*ayrar[i]*vy[i]

## Gráfico
#graph(t,y,"Trajetória da bola", "t (s)", "y (m)")
graph(x,y,"Trajetória da bola", "x (m)", "y (m)")
plt.show()



## Altura máxima e alcance pelo polinómio de Lagrange 
for i in range(N):
    if y[i-1] < y[i] and  y[i+1] < y[i]:
        print('Success')
        maxx, maxy = maximo(x[i-1], x[i], x[i+1], y[i-1], y[i], y[i+1])
        print("A altura máxima atingida foi {:.2f}m após {:.2f}s.".format(maxy, t[i]))
        break

for i in range(N):
    if y[i]*y[i+1] < 0:
        print('Success')
        xzero, yzero = zerosv(x[i-1], x[i], x[i+1], y[i-1], y[i], y[i+1])
        print("O alcance da bola foi {:.2f}m após {:.2f}s.".format(xzero, t[i]))
        break



## Calculo da Energia Mecânica e do trabalho realizado pela Rar num dado Instante xx
for i in range(N):
    if t[i+1] > 0.4-dt and  t[i+1] < 0.4+dt:
        print('A velocidade ao final de 0.4 segundos é {:.2f} m/s no eixo 0x, {:.2f} m/s no eixo 0y, e a altura é {:.2f}m.'.format(vx[i+1], vy[i+1], y[i+1]))
        vv = math.sqrt(vx[i+1]**2 + vy[i+1]**2)

        print("A velocidade tem, então, intensidade de {:.2f} m/s".format(vv))

        Ep = mass*g*y[i+1] 
        Ec = 0.5*mass*(vv**2)

        Em = Ec + Ep

        print("A energia potencial neste instante é dada por {:.2f}J, a cinética por {:.2f}J e a mecânica por {:.2f}J".format(Ep,Ec,Em))
        
        Wx = dt*((fx[0] + fx[i+1])*0.5 + np.sum(fx[1:i]))
        Wy = dt*((fy[0] + fy[i+1])*0.5 + np.sum(fy[1:i]))
        print("O trabalho realizado foi {:.2f}J.".format(Wx+Wy))

        break
