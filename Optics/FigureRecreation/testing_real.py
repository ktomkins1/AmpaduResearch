import sympy as sp
from sympy import sin,cos
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode


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
    
def main_1():
    funcs = setup(10**(-3.5), 0, 4, 0.0002, 0.0002, 1000)
    
    #Define time list and  initial value
    #time = [10**i for i in range(4)] +[0] + [-10**j for j in range(4)]
    #time = sorted(time)
    time_step = 0.0001
    time = np.arange(0,200,time_step)
    init=[0.5,0.5,0,0.25,0.25]

    #Calculate the result
    r = int_ode(init,time,funcs)

    names = ['E1','E2','theta','N1','N2']
    fig, axes = plt.subplots(5,1)
    fig.suptitle(','.join(map(str,init)))
    for k in range(len(r[0])):
        plot_data= [lst[k] for lst in r]
        ax = axes[k]
        ax.title.set_text(names[k])
        ax.plot(time, plot_data)
    plt.show()

def main_2():
    pass

if __name__ == '__main__':
    main_1()
    
    
