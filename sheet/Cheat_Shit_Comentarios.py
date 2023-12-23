import numpy as np
import matplotlib.pyplot as plt

################################ EULER ########################################
 # Parâmetros
dt = 0.01 # δt - tamanho do passo
t0 = 0 # Tempo inicial
tf = 4.0 # Tempo final
y0 = 0 # Posição inicial
vy0 = 0 # Velocidade inicial
g = 9.8 # Aceleração gravítica

# Número de passos/iterações
#
# + 0.1 para garantir que não há arrendodamentos
# para baixo
n = int((tf-t0) / dt + 0.1)

t = np.zeros(n + 1) # Tempo
y = np.zeros(n + 1) # Posição
vy = np.zeros(n + 1) # Velocidade
ay = np.zeros(n + 1) # Aceleração

# Valores inicias
vy[0] = vy0
t[0] = t0
y[0] = y0

for i in range(n):
    ay[i] = g 

    vy[i + 1] = vy[i] + ay[i] * dt
    y[i + 1] = y[i] + vy[i] * dt 
    t[i + 1] = t[i] + dt

################################ EULER CROMER ########################################

# Parâmetros
dt = 0.01 # δt - tamanho do passo
t0 = 0 # Tempo inicial
tf = 4.0 # Tempo final
y0 = 0 # Posição inicial
vy0 = 0 # Velocidade inicial
g = 9.8 # Aceleração gravítica

# Número de passos/iterações
#
# + 0.1 para garantir que não há arrendodamentos
# para baixo
n = int((tf-t0) / dt + 0.1)

t = np.zeros(n + 1) # Tempo
y = np.zeros(n + 1) # Posição
vy = np.zeros(n + 1) # Velocidade
ay = np.zeros(n + 1) # Aceleração

# Valores inicias
vy[0] = vy0
t[0] = t0
y[0] = y0

for i in range(n):
    ay[i] = g 

    vy[i + 1] = vy[i] + ay[i] * dt
    y[i + 1] = y[i] + vy[i + 1] * dt 
    t[i + 1] = t[i] + dt

################################# RESISTENCIA AR #########################################

# Este código faz a simulação da resistência do ar em um corpo em queda livre
# com D coeficiente de arrasto e m massa do corpo e velocidade inicial v0

def forcaResAr(D, v, m):
    v_norm = np.linalg.norm(v)
    v_hat = v / v_norm
    return -m * D * v_norm**2 * v_hat

################################# RESISTENCIA AR CICLISMO #########################################

# Este código faz a simulação da resistência do ar em um ciclista
# com velocidade inicial v0 e coeficiente de arrasto C e a área de seção transversal A

def forcaResAr(C_res, A, v):
    rho_ar = 1.225
    v_norm = np.linalg.norm(v)
    return -C_res/2 * A * rho_ar * v_norm * v

################################# FORÇA GRAVIDADE NEWTON (PLANETAS) #########################################

# Este código faz a simulação da força gravitacional entre dois corpos
# com m massa do corpo 1 e M massa do corpo 2 e r vetor posição do corpo 1 em relação ao corpo 2

def forcaGrav(m, M, r):
    G = 6.67259 * 10**(-11)
    r_norm = np.linalg.norm(r)
    r_hat = r / r_norm
    return G * m * M / r_norm**2 * r_hat

################################# FORÇA ELETROESTATICA #########################################

# Este código faz a simulação da força eletroestática entre dois corpos
# com q carga do corpo 1 e Q carga do corpo 2 e r vetor posição do corpo 1 em relação ao corpo 2

def forcaElet(q, Q, r):
    K = 8.987551 * 10**9
    r_norm = np.linalg.norm(r)
    r_hat = r / r_norm
    return K * q * Q / r_norm**2 * r_hat

################################# FORÇA MAGNETICA #########################################

# Este código faz a simulação da força magnética entre dois corpos
# com q carga do corpo 1 e Q carga do corpo 2 e r vetor posição do corpo 1 em relação ao corpo 2

def forcaMag(q, v, B):
    return q * np.cross(v, B)

################################# FORÇA MAGNUS #########################################

# Este código faz a simulação da força magnus entre dois corpos
# com A coeficiente de arrasto, p densidade do fluido, r raio do corpo, w velocidade angular e v velocidade linear

def forcaMagnus(A, p, r, w, v):
    return 1/2 * A * p * r * np.cross(w, v)

################################# FORÇA ATRITO #########################################

# Este código faz a simulação da força de atrito entre dois corpos
# com mu coeficiente de atrito, N vetor normal e v velocidade linear

def forcaAtrito(v, mu, N):
    v_hat = v / np.linalg.norm(v)
    return -mu * np.linalg.norm(N) * v_hat

