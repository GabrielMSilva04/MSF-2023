import numpy as np
import matplotlib.pyplot as plt
import math

plt.arrow(0,0,0,5,color='r',width=0.05) #(0,0)-> ponto inicial (0,5)-> ponto final
plt.arrow(0,0,5/2,(5*math.sqrt(3))/2,color='b',width=0.05)
plt.arrow(0,0,(5*math.sqrt(3))/2,-5/2,color='g',width=0.05)
plt.arrow(0,0,3.21,3.83,color='y',width=0.05)

plt.show()