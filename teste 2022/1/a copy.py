import numpy as np
import matplotlib.pyplot as plt

# Definir os dados
t = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])
s = np.array([0.1, 1.4, 1.7, 6.5, 7.7, 10.4, 19.5, 26.1, 26.5, 45.9, 52.5])

# Calcular os valores do ajuste linear
n = len(t)
media_t = np.mean(t)
media_s = np.mean(s)
sum_xy = np.sum(t*s) - n*media_t*media_s
sum_xx = np.sum(t**2) - n*media_t**2
declive = sum_xy / sum_xx
ord_origem = media_s - declive*media_t

# Calcular os desvios padrão
y_pred = declive*t + ord_origem
resid = s - y_pred
SSres = np.sum(resid**2)
SStot = np.sum((s - np.mean(s))**2)
s_y = np.sqrt(SSres / (n - 2))
s_x = np.sqrt(sum_xx / (n - 1))

# Calcular os erros do ajuste linear
declive_err = s_y / s_x / np.sqrt(n - 2)
ord_origem_err = s_y * np.sqrt(1/n + media_t**2/sum_xx)

# Plotar o gráfico com a reta de ajuste e os erros
plt.plot(t, s, 'o', label='Dados experimentais')

plt.plot(t, declive*t + ord_origem, 'r', label='Ajuste linear')
plt.xlabel('Tempo (s)')
plt.ylabel('Distância percorrida (cm)')
plt.title('Gráfico de s em função de t')
plt.legend()
plt.show()

# Imprimir os valores do ajuste linear
print(f"Declive: {declive:.2f} +/- {declive_err:.2f}")
print(f"Ordenada na origem: {ord_origem:.2f} +/- {ord_origem_err:.2f}")