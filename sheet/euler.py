import numpy as np

#parametros
g = 9.8
t0 = 0
v0 = 0
vy0 = 0

#inic
def euler(tf, instante, deltat): #deltat = tempo para cada intervalo na integração
    n = np.int_((tf-t0)/deltat)

    t=np.zeros(n+1) #array com tamanho igual ao número de divisões na integração
    y=np.zeros(n+1)
    vy=np.zeros(n+1)
    ay=np.zeros(n+1)

    vy[0]=vy0
    t[0]=t0

    #euler method
    for i in range(n):
        ay[i] = g # queda livre 
        # (em geral pode ser qualquer função de x[i] e vx[i])
        y[i+1]=y[i]+vy[i]*deltat
        vy[i+1]=vy[i] + ay[i]*deltat # atualizar velocidade sabendo aceleração
        t[i+1]=t[i]+deltat
    return t[round(instante/deltat)], y[round(instante/deltat)], vy[round(instante/deltat)], ay[round(instante/deltat)]

t, y, vy, ay = 0, 1, 2, 3

print("tempo final: "+str(euler(4, 3, 0.01)[t]))
print("posição vertical final: "+str(euler(4, 3, 0.01)[y]))
print("velocidade vertical final: "+str(euler(4, 3, 0.01)[vy]))
print("aceleração vertical final: "+str(euler(4, 3, 0.01)[ay]))