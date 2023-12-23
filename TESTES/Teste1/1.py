#nMec-113786
import math
import numpy as np
import matplotlib.pyplot as plt



#a)
t = [0, 48, 96, 144, 192, 240, 288, 336, 384]
ativ = [10.03, 7.06, 4.88, 3.38, 2.26, 1.66, 1.14, 0.79, 0.58]

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
    plt.plot(x, graficoinfo["m"]*x + graficoinfo["b"])
    plt.xlabel('Tempo(horas)')
    plt.ylabel('Ativ(mBq)')
    plt.grid(True)
    plt.show()

x=np.array(t) 
y=np.array(ativ)
graficoinfo = dadosgraficos(x,y)
print("a)")
print("coeficiente de determinação: ",graficoinfo["r2"]) #0.8635533523596588
grafico(x,y)

#RESPOSTA:
#A relação entre o tempo e a atividade não é linear



#b)
logativ= np.log(ativ)

y=np.log(ativ)
graficoinfo = dadosgraficos(x,y)
print("b)")
print("coeficiente de determinação: ",graficoinfo["r2"]) #0.9995011037502572
print("declive:",graficoinfo["m"],"+/-",graficoinfo["deltam"]) #-0.007496883780133508 +/- 6.330605327615247e-05
grafico(x,y)



#c)
print("c)")
meiot = -np.log(2)/graficoinfo["m"]
delta_meiot = (graficoinfo["deltam"] * meiot)/(graficoinfo["m"])
print("meia-vida do decaimento:",meiot,"+/-",abs(delta_meiot))

