import sympy as sp
from sympy import sin,cos

'''
    Create the functions for simulation.  
    These are real python functions to be passed to the integrator.
    
    Single Mode Rate Equations:
    dE/dt   = N*E + eta*cos(psi)
    dpsi/dt = DELTA - b*N - eta*sin(psi)/E
    dN/dt   = [P - N - P*(1 + 2*N)*(E**2)]/T
    
    parameters:
    P - Pumping rate
    DELTA - detuning factor
    b - Linewidth enhancement factor
    eta - Coupling coeff
    T - Ratio of carrier to photon lifetime
'''
def setup(P, DELTA, b, eta, T):
    #Define variables and expressions
    E, psi, N = sp.var('E, psi, N')
    exp1 = E*N + eta*cos(psi)
    exp2 = DELTA - b*N - eta*sin(psi)/E
    exp3 = (P - N - P*(1 + 2*N)*(E**2))/T

    symbols = [E, psi, N]

    #Convert the expressions into lambda functions
    func1 = sp.lambdify(symbols,exp1,"numpy")
    func2 = sp.lambdify(symbols,exp2,"numpy")
    func3 = sp.lambdify(symbols,exp3,"numpy")
    funcs=[func1,func2,func3]
    return funcs
    
def get_names():
    return ['E', 'psi', 'N']
