import numpy as np
import matplotlib.pyplot as plt

# Dados de entrada
x = np.array([1, 2, 3, 4, 5])
y = np.array([2, 4, 5, 4, 5])

# Calculando a média dos valores de x e y
mean_x = np.mean(x)
mean_y = np.mean(y)

# Calculando os valores para a regressão
n = len(x)
numerator = np.sum((x - mean_x) * (y - mean_y))
denominator = np.sum((x - mean_x) ** 2)
b = numerator / denominator
a = mean_y - b * mean_x

# Imprimindo os valores da regressão
print("Coeficiente angular (b):", b)
print("Coeficiente linear (a):", a)

# Imprimindo a equação da reta
print(f"Equação da reta: y = {a} + {b}*x")

# Plotando os dados e a reta de regressão
plt.scatter(x, y, color='blue')
plt.plot(x, a + b*x, color='red')
plt.title('Regressão Linear')
plt.xlabel('x')
plt.ylabel('y')
plt.show()