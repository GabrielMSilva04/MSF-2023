import numpy as np

#Funções para forças

# Força Elétrica
def forcaElet(q, Q, r): #q e Q em Coulombs(carga 1 e 2), r (raio) em metros
    K = 8.987551 * 10**9
    r_norm = np.linalg.norm(r)
    r_hat = r / r_norm
    return K * q * Q / r_norm**2 * r_hat

# Força Magnética
def forcaMag(q, v, B): #q (carga) em Coulombs, v em m/s, B em Tesla
    return q * np.cross(v, B)

# Força de Magnus (Esfera)
def forcaMagnus(A, p, r, w, v):
    return 1/2 * A * p * r * np.cross(w, v)

# Força de Atrito
def forcaAtrito(v, mu, N): #v em m/s, mu (coeficiente de atrito) e N (normal) em Newtons
    v_hat = v / np.linalg.norm(v)
    return -mu * np.linalg.norm(N) * v_hat