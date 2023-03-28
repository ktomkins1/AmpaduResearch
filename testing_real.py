import sympy as sp
from sympy import sin,cos
import numpy as np
from numpy.lib.scimath import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode, solve_ivp


'''
    (1) Phase-locked states of the system:
    dtE_1 = E_1*N_1 - LAMBDA*E_2*sin(theta)
    dtE_2 = E_2*N_2 + LAMBDA*E_1*sin(theta)
    dt_theta = DELTA - alpha*(N_2 - N_1) + LAMBDA*(E_1/E_2 - E_2/E_1)*cos(theta)
    T*dtN_1 = P_1 - N_1 - (1 + 2*N_1)*(E_1**2)
    T*dtN_2 = P_2 - N_2 - (1 + 2*N_2)*(E_2**2)
'''

def setup(lmda,dlt,a,P1,P2,T):
    #Define variables and expressions
    E1, E2, th, N1, N2 = sp.var('E1, E2, th, N1, N2')
    exp1 = E1*N1 - lmda*E2*sin(th)
    exp2 = E2*N2 + lmda*E1*sin(th)
    exp3 = dlt - a*(N2 - N1) + lmda*(E1/E2 - E2/E1)*cos(th)
    exp4 = (P1 - N1 - (1 + 2*N1)*(E1**2))/T
    exp5 = (P2 - N2 - (1 + 2*N2)*(E2**2))/T

    symbols = [E1, E2, th, N1, N2]

    #Convert the expressions into lambda functions
    func1 = sp.lambdify(symbols,exp1,"numpy")
    func2 = sp.lambdify(symbols,exp2,"numpy")
    func3 = sp.lambdify(symbols,exp3,"numpy")
    func4 = sp.lambdify(symbols,exp4,"numpy")
    func5 = sp.lambdify(symbols,exp5,"numpy")
    funcs=[func1,func2,func3,func4,func5]
    return funcs

#Now define the step function(input to the scipy.integrate.ode function)
def int_step(t,y,funcs):
    result=[]
    for i in range(len(y)):
        result.append(funcs[i](*y))
    return result

#Here we define the integrator
def int_ode(y0,time,funcs):
    res=[]
    res.append(y0)
    r = ode(int_step)
    r.set_f_params(funcs)
    r.set_integrator('vode', method='bdf')
    r.set_initial_value(y0,time[0])
    for t in time[1:]:
        r.integrate(t)
        res.append(r.y)
    return res
    
def get_ivals_from_characteristics(E_0, rho, theta_s, alpha, LAMBDA):
    DELTA = -alpha*LAMBDA*np.sin(theta_s)*(rho**-1 + rho) - LAMBDA*np.cos(theta_s)*(rho**-1 - rho)
    
    P_1 = E_0**2 + (1 + 2*E_0**2)*LAMBDA*rho*np.sin(theta_s)
    R2 = rho**2
    P_2 = R2*(E_0**2) - (1 + 2*R2*(E_0**2))*LAMBDA*(rho**-1)*np.sin(theta_s)
    
    print("dlta:{0}, p1:{1}, p2:{2}".format(DELTA, P_1, P_2))
    return DELTA, P_1, P_2

def get_theta_s(s, rho, alpha=4):
    #in case of zero detuning
    return s*np.pi + np.arctan((1/alpha)*((rho**2 - 1)/(rho**2 + 1)))

def get_ref_amplitude(rho, LAMBDA, theta_s):
    #in case of symmetric pumping
    E_0sq = (LAMBDA*np.sin(theta_s)*(rho**2 + 1))/(rho*((rho**2 - 1) - 4*LAMBDA*rho*np.sin(theta_s)))
    return sqrt(E_0sq)
    
def get_sym_pr(E_0, LAMBDA, rho, theta_s):
    return E_0**2 + (1 + 2*E_0**2)*LAMBDA*rho*np.sin(theta_s)

def main(E_0, rho, theta_s, logLAMBDA, alpha=4, 
         llsim=0, ulsim=10000, llcyc=3000, ulcyc=4000,
         override_ivals=None, override_E_0=None, 
         override_theta_s=None):
    #constants!
    tau_p = 2*10**-3
    LAMBDA = 0#10**logLAMBDA
    
    if override_theta_s != None:
        theta_s = get_theta_s(override_theta_s, rho, alpha=alpha)
        print("theta_s:{0}".format(theta_s))
        
    if override_ivals == None:
        DELTA, P_1, P_2 = get_ivals_from_characteristics(E_0, rho, theta_s, alpha, LAMBDA)
    else:
        DELTA, P_1, P_2 = override_ivals
    
    
    if override_E_0 != None:
        E_0 = get_ref_amplitude(rho, LAMBDA, theta_s)
        P1 = get_sym_pr(E_0, LAMBDA, rho, theta_s)
        P2 = P1
    
    #Get the system of equations
    funcs = setup(LAMBDA, DELTA, alpha, P_1, P_2, 1000)
    
    #Define time list and  initial value
    init=[np.sqrt(0.5),np.sqrt(0.5),0,0,0]
    
    r, t = time_and_result_2(init, funcs, llsim, ulsim)
    time = t*tau_p
    
    names = ['E1','E2','theta','N1','N2']
    fig, axes = plt.subplots(5,1, sharex='col')
    fig.suptitle(','.join(map(str,init)))
    for k in range(len(r)):
        plot_data= r[k]
        ax = axes[k]
        ax.set_title(names[k])
        ax.set_ybound(lower=-1.0, upper=1.0)
        if k == 4: ax.set_xlabel("Time [ns]")
        ax.set_ylabel("Amplitude")
        ax.plot(time, plot_data)
        if k in [0,1]:
            ax.axvline(x=llcyc*tau_p, color='r')
            ax.axvline(x=ulcyc*tau_p, color='r')
    #plt.show()
    
    fig2, axes2 = plt.subplots(1)
    ax = axes2
    ax.plot(r[0][llcyc:ulcyc], r[1][llcyc:ulcyc])
    ax.set_title("Limit Cycle")
    ax.set_xlabel("E_1")
    ax.set_ylabel("E_2")
    plt.show()

def time_and_result_1(init, funcs, ll, ul):
    #Define time list and  initial value
    #time = [10**i for i in range(4)] +[0] + [-10**j for j in range(4)]
    #time = sorted(time)
    time_step = 0.0001
    time = np.arange(ll, ul,time_step)
    
    #Calculate the result
    r = int_ode(init,time,funcs)
    
    ret = [[lst[k] for lst in r] for k in range(5)]
    return ret, time

def time_and_result_2(init, funcs, ll, ul):
    #Calculate the result
    res_map = solve_ivp(int_step, (ll, ul), init, args=(funcs,), dense_output=True, max_step=1.0)
    return res_map['y'], res_map['t']
    

if __name__ == '__main__':
    main(0.5, 0.20, 1.89, -2.2, llcyc=6000, ulcyc=7000)
    main(0.5, 0.24, 2.75, -2.2)
    main(0.5, 0.05, 3.1415, -2.2, llcyc=6000, ulcyc=7000)
    main(0.5, 0.005, 3.1415, -2.2, llcyc=6000, ulcyc=7000)
    
    
