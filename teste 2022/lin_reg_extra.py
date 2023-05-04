#Calcular a Regressão Linear
lin_reg = stats.linregress(x, y)

'''
m = lin_reg.slope
b = lin_reg.intercept
r²= lin_reg.rvalue
'''
print("d={:.3f}x + {:.3f}".format(lin_reg.slope, lin_reg.intercept))
print("r²={:.3f}".format(lin_reg.rvalue))