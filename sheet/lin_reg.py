import math
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from scipy import stats


#retona um dicionario com dados de graficos lineares (keys: m , b , r, deltam, deltab)
def dadosgraficos(x,y):
    n=len(x)
    sumx=sum(x)
    sumy=sum(y)
    xy=sum(x[i]*y[i] for i in range(n))
    xx=sum(x[i]**2 for i in range(n))
    yy=sum(y[i]**2 for i in range(n))

    m=((n*xy)-(sumx*sumy))/((n*xx)-(sumx**2))
    b=((xx*sumy)-(sumx*xy))/((n*xx)-(sumx**2))
    r2=(((n*xy)-(sumx*sumy))**2)/(((n*xx)-(sumx**2))*((n*yy)-(sumy**2)))
    deltam=abs(m)*math.sqrt(((1/r2)-1)/(n-2))
    deltab= deltam* math.sqrt(xx/n)

    return {'m':m,'b':b,'r2':r2,'deltam':deltam,'deltab':deltab}

def grafico(x, y):
    plt.plot(x, y, 'ro')
    plt.xlabel('----')
    plt.ylabel('------')
    plt.title('Exercício ---- MSF')
    plt.grid(True)
    plt.show()

x=[222.0,207.5,194.0,171.5,153.0,133.0,113.0,92.0] 
y=[2.3,2.2,2.0,1.8,1.6,1.4,1.2,1.0]
print(dadosgraficos(x,y)) #{'m ': 0.010155051683894608, 'b ': 0.05507544181394264, 'r2 ': 0.9984571397353, 'deltam ': 0.00016296903598723278, 'deltab ': 0.027130765543908482}
print("declive: ",dadosgraficos(x,y)["m"]) #declive:  0.010155051683894608
grafico(x,y)

#c1=np.polyfit (x,y,grau)
#c1[1] #declive
#c1[0] #ordenada na origem

#-----------funçoes logaritmicas-------------
m= [0.15,0.20,0.16,0.11,0.25,0.32,0.40,0.45,0.50,0.55]
T = [1.21,1.40,1.26,1.05,1.60,1.78,2.00,2.11,2.22,2.33]

logt= np.log(T)
logm= np.log(m)

#valores negativos
xmax=np.max(logm)+np.abs((np.min(logm)*0.1))
xmin=np.min(logm)-np.abs((np.min(logm)*0.1))

lin_reg = stats.linregress(logm, logt) #lista com o declive, ordenada, valor de r, ?, delta m, deltab 
print(lin_reg)
x = np.linspace(xmin, xmax, 100)#Gerar array com 100 valores 
plt.plot(x, x * lin_reg.slope + lin_reg.intercept)
plt.xlabel('Tempo(dias)')
plt.ylabel('log/Atividade (mCi)')
plt.title('Exercício 09 MSF')
plt.plot(logm, logt, "ro")
plt.show()

