# Magnus -  Aula 7: 4.c.py
# Ondas - Aula 7: 6.b.py
# Em - Aula 8:  3.a.py ou 3.b.py(Com resistencia)
# Integrais para trabalho - Aula 8:  3.c.py
# Potencia de ciclista - Aula 9: 11.a.py

# Altura maxima
for i in range(n-1):
    if (y[i+1] < y[i]):
        print("A altura max é {:0.2f} m".format(y[i+1]))
        plt.plot(x[i+1], y[i+1], "o", label="Altura máxima")
        break

# Alcance
for i in range(n-1):
    if (y[i+1]*y[i] < 0):
        print("Alcance {:0.2f} m".format(x[i+1]))
        plt.plot(x[i+1], y[i+1], "o", label="Alcance")
        break

    # ver quando atinge um valor
    if x[i] < 2000 + dt and x[i+1] > 2000-dt:
        print("Tempo quando atinge 2km = ", t[i])

    # velocidade terminal
    if (vx[i]-vx[i-1] < tolerancia):
        vt = vx[i]
        print("Velociadade Terminal = ", vx[i])
    break

# % de velocidade terminal
for i in range(n-1):
    if ((0.90*vt-dt) < v[i] < (0.90*vt+dt)):
        print("Tempo: ", t[i])
        break
