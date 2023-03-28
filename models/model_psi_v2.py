import sympy as sp
from sympy import sin,cos
import numpy as np

'''
    Create the functions for simulation.
    These are real python functions to be passed to the integrator.

    Single Mode Rate Equations:
    dE/dt   = N*E + eta*cos(psi + DELTA*t)
    dpsi/dt = - b*N - eta*sin(psi + DELTA*t)/E
    dN/dt   = [P - N - P*(1 + 2*N)*(E**2)]/T

    parameters:
    P - Pumping rate
    DELTA - detuning factor
    b - Linewidth enhancement factor
    eta - Coupling coeff
    T - Ratio of carrier to photon lifetime
    t - time
'''
def setup(P, DELTA, b, eta, T):
    #Define variables and expressions
    E, psi, N, t = sp.var('E, psi, N, t')
    OMEGA = DELTA*np.sqrt(2*P/T)
    dE   = E*N + eta*cos(psi + OMEGA*t)
    dpsi = -b*N - eta*sin(psi + OMEGA*t)/E
    dN   = (P - N - P*(1 + 2*N)*(E**2))/T

    symbols = [t, E, psi, N]

    #Convert the expressions into lambda functions
    func1 = sp.lambdify(symbols,dE,"numpy")
    func2 = sp.lambdify(symbols,dpsi,"numpy")
    func3 = sp.lambdify(symbols,dN,"numpy")
    funcs=[func1,func2,func3]
    return funcs

def setup_ks(P, DELTA, b, eta, T):
    n = 1
    args = [P, DELTA, b, eta, T]
    for i, a in enumerate(args):
        try:
            n_ = len(a)
            n = n_
            break
        except TypeError:
            pass

    if n != 1:
        for i, a in enumerate(args):
            try:
                n_ = len(a)
                if n_ != n:
                    pass #should error
            except TypeError:
                args[i] = args[i]*np.ones(n)

    P, DELTA, b, eta, T = args
    #Define variables and expressions
    E, psi, N, t = sp.var('E, psi, N, t')
    OMEGA = DELTA*np.sqrt(2*P/T)
    dE   = E*N + eta*cos(psi + OMEGA*t)
    dpsi = -b*N - eta*sin(psi + OMEGA*t)/E
    dN   = (P - N - P*(1 + 2*N)*(E**2))/T

    symbols = [t, E, psi, N]

    #Convert the expressions into lambda functions
    func1 = sp.lambdify(symbols,dE,"numpy")
    func2 = sp.lambdify(symbols,dpsi,"numpy")
    func3 = sp.lambdify(symbols,dN,"numpy")
    funcs=[func1,func2,func3]
    return funcs

def setup_k(P, DELTA, b, eta, T):
    n = 1
    args = [P, DELTA, b, eta, T]
    for i, a in enumerate(args):
        try:
            n_ = len(a)
            n = n_
            break
        except TypeError:
            pass

    if n != 1:
        for i, a in enumerate(args):
            try:
                n_ = len(a)
                if n_ != n:
                    pass #should error
            except TypeError:
                args[i] = args[i]*np.ones(n)

    P, DELTA, b, eta, T = args
    #Define variables and expressions
    E, psi, N, t = sp.var('E, psi, N, t')
    OMEGA = DELTA*np.sqrt(2*P/T)
    dE   = lambda E, psi, N, t: E*N + eta*np.cos(psi + OMEGA*t)
    dpsi = lambda E, psi, N, t: -b*N - eta*np.sin(psi + OMEGA*t)/E
    dN   = lambda E, psi, N, t: (P - N - P*(1 + 2*N)*(E**2))/T

    funcs=[dE,dpsi,dN]
    return funcs

def get_names():
    return ['E', 'psi', 'N']

shortname = 'psi2'
longname = 'model_psi_v2'
