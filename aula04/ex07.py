import numpy as np
import sympy as sp


v0 = 10
a = -9.8
y0 = 0

t, r = sp.symbols('t, r')

#s/ resistencia
y = y0 + v0 * t + 0.5 * a*t**2
dy = sp.diff(y)

tempomax = sp.nsolve(dy-0, t, 0)
tempoy0 = sp.nsolve(y-0, t, 1)

y = sp.lambdify(t, y, 'numpy')
dy = sp.lambdify(t, dy, 'numpy')

posicaomax = y(tempomax)

print("Altura máxima: "+ str(posicaomax)) #b)
print("Tempo em que a bola passa na origem :"+ str(tempoy0)) #c)

#c/ resistencia
#y = y0 + v0 * t + 0.5 * (a-r)*t**2
#dy = sp.diff(y)

#vt = 100 / 3.6
#posicaovt = sp.nsolve(dy-vt, t, 0) #v = dy
#print(posicaovt)
