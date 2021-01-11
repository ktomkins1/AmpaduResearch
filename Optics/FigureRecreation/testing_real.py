import sympy as sp
from sympy import sin,cos
import numpy as np
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

def main():
    #Get the system of equations
    LAMBDA = 10**(-2.2)
    alpha = 4
    theta_s = 2.75
    rho = 0.24
    DELTA, P_1, P_2 = get_ivals_from_characteristics(0.5, rho, theta_s, alpha, LAMBDA)
    funcs = setup(LAMBDA, DELTA, alpha, P_1, P_2, 1000)
    
    #Define time list and  initial value
    init=[0.5,0.5,0,0,0]
    
    r, time = time_and_result_1(init, funcs, 0, 4000)

    names = ['E1','E2','theta','N1','N2']
    fig, axes = plt.subplots(5,1)
    fig.suptitle(','.join(map(str,init)))
    for k in range(len(r)):
        plot_data= r[k]
        ax = axes[k]
        ax.title.set_text(names[k])
        ax.plot(time, plot_data)
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
    res_map = solve_ivp(int_step, (ll, ul), init, args=(funcs,))
    return res_map['y'], res_map['t']
    

if __name__ == '__main__':
    main()
    
    
