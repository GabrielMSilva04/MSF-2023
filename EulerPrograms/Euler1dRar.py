import matplotlib.pyplot as plt
import numpy as np

# Função utilizada para obter gráficos simples
def graph(x,y,title = "Gráfico", xaxis = "Eixo X", yaxis = "Eixo Y"):
    plt.plot(x,y,'o')

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
v0 = #m/s

g = 9.80 ##aceleração gravitica

## Passo temporal para o Método de Euler
dt = 0.001 ##passo temporal
t0 = 0.00 ##tempo inicial
tf = 5.00 ##tempo final
N = np.int((tf-t0)/dt) ##Nº de passos a realizar
print("Número de passos: ", N)

## Resistência do Ar
vt =  #m/s²
D = g/vt**2

## Criação de arrays de velocidade, posição e tempo com N lugares, para os N passos,
#  preenchidos por 0's. Caso a velocidade ou posição inicial seja diferente de 0, temos que fazer a correção.

t = np.zeros(N+1)

ax = np.zeros(N+1)
vx = np.zeros(N+1)
x = np.zeros(N+1)

vx[0] = v0
x[0] = x0



## Ciclo que implementa o Método de Euler, e atualiza o array tempo para podemos fazer um plot
for i in range(N):
    t[i+1] = t[i] + dt

    vv = np.abs(vx[i])
    ax[i+1] = g - D*vx[i]*vv

    vx[i+1] = vx[i] + ax[i]*dt
    x[i+1] = x[i] + vx[i]*dt

## Gráfico posição*tempo

graph(t,x,"Método de Euler","Tempo (s)","Posição (m)")
plt.show()

## Determinação de altura máxima e alcance com recurso ao polinómio de Lagrange
for i in range(N):
    if y[i-1] < y[i] and  y[i+1] < y[i]:
        print('Success')
        maxx, maxy = maximo(t[i-1], t[i], t[i+1], y[i-1], y[i], y[i+1])
        print("A altura máxima atingida foi {:.2f}m após {:.2f}s.".format(maxy,maxx))
        break

for i in range(N):
    if y[i]*y[i-1] < 0:
        print('Success')
        xzero, yzero = zerosv(t[i-1], t[i], t[i+1], y[i-1], y[i], y[i+1])
        print("A bola voltou a passar em 0m após após {:.2f}s.".format(xzero))
        break