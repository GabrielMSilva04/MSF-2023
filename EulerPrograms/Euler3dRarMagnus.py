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
z0 = 


alfa = 10*math.pi/180
v0 = xx*1000/3600 ##velocidade inicial (m/s)
vx0 = v0*math.cos(alfa)
vy0 = v0*math.sin(alfa)
vz0 = 0

g = 9.80 ##aceleração gravitica

## Com resistência do ar
vt = xx*1000/3600
D = g/vt**2  ##Cálculo de D para o efeito da Resistência do Ar

## Implementando a Força de Magnus
r =  ##raio da bola (m)
dAr = 1.225
A = (math.pi)*(r**2)
Mx0 = ## Calcular com base no teorema de laplace. Não esquecer de meter fórmula completa!
My0 =
Mz0 =
m = ##massa da bola (kg)

## Passo temporal para o Método de Euler
dt = 0.001 ##passo temporal
t0 = 0.00 ##tempo inicial
tf = 5.00 ##tempo final
N = np.int((tf-t0)/dt) ##Nº de passos a realizar
print("Número de passos: ", N)

## Criação de arrays de velocidade, posição e tempo com N lugares, para os N passos,
#  preenchidos por 0's. Caso a velocidade ou posição inicial seja diferente de 0, temos que fazer a correção.

## NOTA: Vx é constante s/ resistência do ar (igual a Vx0)
t = np.zeros(N+1)

ax = np.zeros(N+1)
ay = np.zeros(N+1)
az = np.zeros(N+1)

vx = np.zeros(N+1)
vy = np.zeros(N+1)
vz = np.zeros(N+1)

x = np.zeros(N+1)
y = np.zeros(N+1)
z = np.zeros(N+1)

Mx = np.zeros(N+1)
My = np.zeros(N+1)

vz[0] = vz0
vy[0] = vy0
vx[0] = vx0

x[0] = x0
y[0] = y0
z[0] = z0

Mx[0] = Mx0
My[0] = My0

## Ciclo que implementa o Método de Euler, e atualiza o array tempo para podemos fazer um plot
for i in range(N):
    t[i+1] = t[i] + dt

    vv = math.sqrt(vx[i]**2 + vy[i]**2) ## atualizar caso tenhamos vz

    Mx[i+1] = 0.5*A*r*dAr* ##completar com formula retirada do teorema de laplace
    My[i+1] = 0.5*A*r*dAr*

    ay[i] = -( g + (D*vy[i]*vv)) + My[i]/m
    ax[i] = -D*vx[i]*vv + Mx[i]/m

    vz[i+1] = vz[i] + az[i]*dt
    vy[i+1] = vy[i] + ay[i]*dt
    vx[i+1] = vx[i] + ax[i]*dt

    x[i+1] = x[i] + vx[i]*dt
    y[i+1] = y[i] + vy[i]*dt
    z[i+1] = z[i] + vz[i]*dt

graph(x,y,"Altura x Distância", "Distância (m)", "Altura (m)")
plt.show()


## Altura Máxima e Alcance de acordo com o polinómio de Lagrange 
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
    
for i in range(n):
    if y[i]*y[i+1] < 0:
        alcance = x[i]

print('Alcance: ', alcance)