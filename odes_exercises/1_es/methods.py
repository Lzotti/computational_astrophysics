import numpy as np
import scipy.optimize as opt


# scrivo una introduzione decente, pro e contro di ogni metodo e come sono implementati in scipy e commento
#problema: alcuni metodi vanno più avanti di altri

def forward_euler(t, y, h, f, args = None, extra_par=None):

    y_next = y + h * f(t, y, args)

    return y_next   

# implemento io un root finding, e uso pure newton con derivata.
def backward_euler(t, y, h, f, args = None, extra_par=None):

    y_next = opt.root(lambda y_next: y_next - y - h * f(t + h, y_next, args), y).x
    #sto attento a .x o .x[0] o .x.item
    return y_next   


def explicit_midpoint(t, y, h, f, args = None, extra_par=None):

    y_mid = y + h/2 * f(t, y, args)
    y_next = y + h * f(t + h/2, y_mid, args)

    return y_next

def Heun(t, y, h, f, args = None, extra_par=None):

    y_pred = y + h * f(t, y, args)
    y_next = y + h/2 * (f(t, y, args) + f(t + h, y_pred, args))

    return y_next

# implemento io un root finding, e uso pure newton con derivata.
def implicit_midpoint(t, y, h, f, args = None, extra_par=None):
    
    y_next = opt.root(lambda y_next: y_next - y - h * f(t + h/2, (y + y_next)/2, args), y).x

    return y_next

def RK4(t, y, h, f, args = None, extra_par=None):

    k1 = f(t, y, args)
    k2 = f(t + h/2, y + h/2 * k1, args)
    k3 = f(t + h/2, y + h/2 * k2, args)
    k4 = f(t + h, y + h * k3, args)

    y_next = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)

    return y_next

#rivedo se posso implementare meglio f_prev 
def AB2(t, y, h, f, args = None, f_prev=None):

    y_next = y + h/2 * (3 * f(t, y, args) - f_prev)

    return y_next

#miglioro implementazione f prev, e mi assicuro sia un array della dimensione giusta e ste cose qua
def AB4(t, y, h, f, args = None, f_prev=None):

    y_next = y + h/24 * (55 * f(t, y, args) - 59 * f_prev[0] + 37 * f_prev[1] - 9 * f_prev[2])

    return y_next

def AM1(t, y, h, f, args = None, extra_par=None ):

    y_next = opt.root(lambda y_next: y_next - y - (h/2) * (f(t, y, args) + f(t + h, y_next, args)), y).x

    return y_next

#rivedo molto implementazione
def BDF2(t, y, h, f, args = None, y_prev=None):

    y_next = opt.root(lambda y_next: y_next - (4/3) * y + (1/3) * y_prev - (2*h/3) * f(t + h, y_next, args), y).x

    return y_next

#ora lo sto implementando a caso, dovrei aggiungere controlli o eccezioni sulla forma di y e su quella di f
def leapfrog_KDK(t, y, h, f, args = None, extra_par=None):

    m = args[0] # massa
    if m is None:
        m = 1.0 # default

    dv_dq = args[1] # derivata di V rispetto a q

    v_half = y[1] - h/2 * dv_dq(t, y, args) # calcolo v a metà passo, e lo uso per calcolare q al passo successivo, e poi v al passo successivo.
    y_next_0 = y[0] + h * v_half/m
    y_next_1 = v_half - h/2 * f(t + h, [y_next_0, 0], args) # anche qua mi sembra implementata a caso

    y_next = np.array([y_next_0, y_next_1])

    return y_next

    
