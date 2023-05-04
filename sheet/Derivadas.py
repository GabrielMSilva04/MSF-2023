import math
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from scipy import stats

v_terminal = sym.Symbol('vt', real=True, positive=True)
gravidade = sym.Symbol('g', real=True, positive=True)
temp = sym.Symbol('tempo', real=True, positive=True)
eq_velocidade = sym.Symbol('v(t)', real=True, positive=True)
eq_aceleracao = sym.Symbol('a(t)', real=True, positive=True)

print("\nVELOCIDADE")

#calcular a equação da velocidade instantanea
eq_velocidade = sym.diff((v_terminal ** 2 / gravidade) * sym.log(sym.cosh((gravidade * temp) / v_terminal)), temp, 1)
print(eq_velocidade)

#simplificar expressao da velocidade
eq_velocidade = sym.simplify(eq_velocidade)
print(eq_velocidade)


print("\nACELERAÇÃO")

#calcular a equação da aceleração instantanea
eq_aceleracao = sym.diff((v_terminal ** 2 / gravidade) * sym.log(sym.cosh((gravidade * temp) / v_terminal)), temp, 2)
print(eq_aceleracao)

#simplificar expressao da velocidade
eq_aceleracao = sym.simplify(eq_aceleracao)
print(eq_aceleracao)

#sol=sym.solve(2 + 4*t_)
