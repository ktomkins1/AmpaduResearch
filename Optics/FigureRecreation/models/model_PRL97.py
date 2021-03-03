import sympy as sp
from sympy import exp, I
import numpy as np

'''
    Create the functions for simulation.
    These are real python functions to be passed to the integrator.

    Single Mode Rate Equations:
    dE/dt   = (1 - j*b)*N*E + eta*exp(-j*DELTA*t)
    dN/dt   = [P - N - (1 + 2*N)*(abs(E)**2)]/T

    parameters:
    P - Pumping rate
    DELTA - detuning factor
    b - Linewidth enhancement factor
    eta - Coupling coeff
    T - Ratio of carrier to photon lifetime
    t - time
'''
def setup(P, DELTA, b, eta, T):
    #print('arguments are: {}, {}, {}, {}, {}'.format(P, DELTA, b, eta, T))
    #Define variables and expressions
    E, N, t = sp.var('E, N, t')
    dE   = (1 - I*b)*N*E + eta*exp(-I*DELTA*t)
    dN   = (P - N - (1 + 2*N)*(abs(E)**2))/T

    symbols = [t, E, N]

    #Convert the expressions into lambda functions
    func1 = sp.lambdify(symbols,dE,"numpy")
    func2 = sp.lambdify(symbols,dN,"numpy")
    funcs=[func1,func2]
    return funcs

def get_names():
    return ['E', 'N']
