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
def setup(P, DELTA, b, eta, T, get_exp=False):
    #Define variables and expressions
    E, psi, N, t = sp.var('E, psi, N, t')
    exp1 = E*N + eta*cos(psi + DELTA*t)
    exp2 = -b*N - eta*sin(psi + DELTA*t)/E
    exp3 = (P - N - P*(1 + 2*N)*(E**2))/T

    symbols = [t, E, psi, N]

    #Convert the expressions into lambda functions
    func1 = sp.lambdify(symbols,exp1,"numpy")
    func2 = sp.lambdify(symbols,exp2,"numpy")
    func3 = sp.lambdify(symbols,exp3,"numpy")
    funcs=[func1,func2,func3]
    if get_exp: funcs += [exp1,exp2,exp3]
    return funcs
    
def get_names():
    return ['E', 'psi', 'N']

if __name__ == '__main__':
    eta = np.linspace(0, 0.5, 100)
    funcs = setup(9.3, 0.1, 4.8, eta, 958.0, get_exp=True)
    print('our functions are:')
    for f in funcs[3:]:
        print(str(f))