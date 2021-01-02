import numpy as np
import matplotlib

def sin(x):
    return np.sin(x)

def cos(x):
    return np.cos(x)

'''
    the following quantities describe the system:
    E_i         the E-field of the i-th block*                          
    dtX         the time derivative of the X (dtE_i would be of E_i)
    N_i         the carrier density of the i-th block*
    theta       the phase difference between blocks in a 2-block system
    alpha       the line-width enhancement factor    [typ 4]
    LAMBDA      the coupling constant*    [typ 10^-4 to 10^0]
    P_i         the excess pumping rate for the i-th block*
    DELTA       the cavity optical detuning
    T           the ratio of carrier lifetime to photon lifetime    [typ 10^3]
    t           the time of the system, normalized to photon lifetime
    tau_p       the photon lifetime    [typ 2x10^-3 ns]
    tau_c       the carrier lifetime    [typ 2 ns]
    
    * normalized
    
    (E_0, rho, theta_s) characterization of a stable phase-locked states:
            E_0 == E_1
            rho == E_2/E_1
    theta_s     the phase difference theta for a stable state numbered s
'''

'''
    the following functions govern the evolution of the system.
    (written similar to python expressions)
    
    (1) Phase-locked states of the system:
    dtE_1 = E_1*N_1 - LAMBDA*E_2*sin(theta)
    dtE_2 = E_2*N_2 + LAMBDA*E_1*sin(theta)
    dt_theta = DELTA - alpha*(N_2 - N_1) + LAMBDA*(E_1/E_2 - E_2/E_1)*cos(theta)
    T*dtN_1 = P_1 - N_1 - (1 + 2*N_1)*(E_1**2)
    T*dtN_2 = P_2 - N_2 - (1 + 2*N_2)*(E_2**2)
    
    (2) determination of DELTA based on characterization (E, r, t):
    DELTA = -1*alpha*LAMBDA*sin(theta_s)*(rho**-1 + rho) +\
            -1*LAMBDA*cos(theta_s)*(rho**-1 - rho)
    
    (3) determination of P1, P2:
    P_1 = E_0**2 + (1 + 2*(E_0**2))*LAMBDA*rho*sin(theta_s)
    P_2 = (rho**2)*(E_0**2) - (1 + 2*(rho**2)*(E_0**2))*\
                              LAMBDA*(rho**-1)*sin(theta_s)
    
    (4) steady state values of N_i
    N_1 = LAMBDA*rho*sin(theta_s)
    N_2 = -1*LAMBDA*(rho**-1)*sin(theta_s)
    
    (5) restriction if DELTA=0 for theta_s:
    theta_s = s*pi + arctan((1/alpha)*((rho**2 - 1)/(rho**2 + 1)))
    s = 0,1
    
    (6) restriction if P_1 = P_2 = P_0:
    E_0**2 = (LAMBDA*sin(theta_s)*(rho**2 + 1))/\
             (rho*((rho**2 - 1) - 4*LAMBDA*rho*sin(theta_s)))
    P_0 = E_0**2 + (1 + 2*(E_0**2)*LAMBDA*rho*sin(theta_s))
             
    
'''



if __name__ == '__main__':
    
    pass
    