import numpy as np
import matplotlib.pyplot as plt
import math

def graph(x,y,title = "Gráfico", xaxis = "Eixo X", yaxis = "Eixo Y"):
    plt.plot(x,y,'o')

    plt.xlabel(xaxis)
    plt.ylabel(yaxis)

    plt.grid()

    plt.title(title)

x0 = 
v0 = #m/s

g = 9.80 ##aceleração gravitica

## Passo temporal para o Método de Euler
dt = 0.001 ##passo temporal
t0 = 0.00 ##tempo inicial
tf = 5.00 ##tempo final
N = np.int((tf-t0)/dt) ##Nº de passos a realizar
print("Número de passos: ", N)

## Criação de arrays de velocidade, posição e tempo com N lugares, para os N passos,
#  preenchidos por 0's. Caso a velocidade ou posição inicial seja diferente de 0, temos que fazer a correção.

t = np.zeros(N+1)

vx = np.zeros(N+1)
x = np.zeros(N+1)

vx[0] = v0
x[0] = x0


## Ciclo que implementa o Método de Euler, e atualiza o array tempo para podemos fazer um plot
for i in range(N):
    t[i+1] = t[i] + dt

    vx[i+1] = vx[i] - g*dt

    x[i+1] = x[i] + vx[i]*dt

graph()
plt.show()

##Ifs para determinar velocidade/posição após x segundos
## Caso não seja encontrado nenhum resultado, multiplicar valor de dt por 2, 5, 10, 15.... até ser encontrado
for i in range(N):
    if t[i+1] > 3-dt and  t[i+1] < 3+dt:
        print('A velocidade ao final de 3 segundos é {:.1f} m/s, calculando recorrendo-se ao método de Euler com um passo de {:.3f}s.'.format(vx[i+1], dt))
        break

for i in range(N):
    if t[i+1] > 2-dt and  t[i+1] < 2+dt:
        print('A posição ao final de 2 segundos é {:.1f} m, calculando recorrendo-se ao método de Euler com um passo de {:.3f}s.'.format(x[i+1], dt))
        break