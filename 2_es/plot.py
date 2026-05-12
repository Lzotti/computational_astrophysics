import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**2



trapezio = np.loadtxt("trap.dat")
simpson = np.loadtxt("simp.dat")

x_t = trapezio[:, 0]
trapezio_dat = trapezio[:, 1]

x_s = simpson[:, 0]
simpson_dat = simpson[:, 1]


plt.figure(figsize=(10, 6))

plt.plot(x_t, trapezio_dat, label="trapezio", color='blue')
plt.plot(x_s, simpson_dat, label="simpson", color='red')
plt.legend()
plt.xlabel("N")
plt.ylabel("Errore")
plt.xscale("log")
plt.yscale("log")
plt.title("Errore in funzione di N per i vari metodi di integrazione")
plt.show()


