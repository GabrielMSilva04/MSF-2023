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

print("alínea b(deltat=0.01): "+str(euler(4, 3, 0.01)[vy])) #ex b)
print("alínea c(deltat=0.001): "+str(euler(4, 3, 0.001)[vy])) #ex c)

print("alínea e(posição em t=2): "+str(euler(3, 2, 0.01)[y])) #ex e)
print("alínea f(passo 10x menor): "+str(euler(3, 2, 0.001)[y])) #ex f)

print("alínea h(passo 100x menor): "+str(euler(3, 2, 0.0001)[y])) #ex h)


