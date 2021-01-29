import sympy as sp
from sympy import sin,cos
import numpy as np
from numpy.lib.scimath import sqrt
import matplotlib.pyplot as plt
from scipy.integrate import ode, solve_ivp
from scipy.signal import find_peaks


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

'''
    One single frame of integration
    
    based on previous values, calculate next values of the system
    
    parameters:
    t - the single time value
    y - the vectorized current state of the system [s1, s2, ... sn]
    funcs - the python functions corresponding to each diff eq F1..n
'''
def int_step(t,y,funcs):
    result=[]
    for i in range(len(y)):
        result.append(funcs[i](*y))
    return result

'''
    Perform simulation of the system using the Initial Value Problem approach
    
    Integrator uses Runge Kutta Method by default
    
    parameters:
    initial values - initial vector of the system
    funcs - the functions above (F1...n)
    ll - the initial time
    ul - the ending time
    step (default 1.0) the maximum time step for the integrator to use
    
    returns 2-tuple: 
    0 - list of traces created [s1(t), s2(t), ... sn(t)]
    1 - the associated time steps
'''
def simulate_functions(initial_values, funcs, ll, ul, step=1.0):
    #Calculate the result
    res_map = solve_ivp(int_step, [ll, ul], initial_values,
                        args=(funcs,), max_step=step, rtol=5e-9)
    return res_map['y'], res_map['t']

'''
    Perform frequency analysis on a single time series trace
    
    parameters:
    trace - the time series to perform on
    time_trace - the time steps associated with each point
    dt - the sampling period of the time series
'''
def freq_analysis_trace(trace, time_trace, dt):
    return np.abs(np.fft.fftshift(np.fft.fft(trace)))

'''
    For a single trace, establish a baseline and subtract it out
    
    This allows measuring the trace above 0 rather than a (moving) baseline
    
    For now, uses boxcar method (moving avg)
    
    parameters:
    trace - the time series
    width - the width of the boxcar. should be sufficient for the trace freqs
'''
def normalize_trace(trace, width):
    rolling_avg = np.convolve(np.ones(width)/width, trace, mode='same')
    return trace - rolling_avg

'''
    Create a plot for each of the associated traces in a list
    
    parameters:
    traces - list of n traces
    time_trace - time steps to plot against
    trace_names - list of names, one for each of n traces
    title - title of the plot
    Ebounds - for any traces appearing to be E-fields, the amplitude bounds
'''
def plot_n_traces(traces, time_trace, trace_names, title, Ebounds=[-1.0,1.5]):
    num_traces = len(traces)
    
    fig, axes = plt.subplots(num_traces,1,sharex='col')
    fig.suptitle(title)
    
    for k in range(num_traces):
        ax = axes[k]
        name = trace_names[k]
        if k == num_traces-1: 
            ax.set_xlabel("Normalized Time")
        ax.set_ylabel(name)
        ax.plot(time_trace, traces[k])
        if name[0] == 'E':
            # this is an E trace so use the bounds
            ax.set_ybound(lower=Ebounds[0], upper=Ebounds[1])

'''
    Plot one trace against another - useful for limit cycles
    
    parameters:
    trace1 - the x axis
    trace2 - the y axis
    name1 - the name associated with trace1
    name2 - the name associated with trace2
'''
def plot_limit_cycle(trace1, trace2, name1, name2):
    fig, axes = plt.subplots()
    fig.suptitle("Limit Cycle")
    axes.plot(trace1, trace2)
    axes.set_xlabel(name1)
    axes.set_ylabel(name2)

'''
    Plot the analysis of a single set of input parameters
    
    Simulates the system for a set type of initial conds and given input params
    
    parameters:
    E_0 - the references E-field
    P - the pumping rate
    DELTA - the detuning factor
    eta - the coupling coeff
    b (nom 4) - the linewidth enhancement factor
    T (nom 1000) - carrier : photon lifetime
    llsim - simulation starting point
    ulsim - simulation ending point
    llcyc - for creating a window into traces, the lower bound
    ulcyc - for creating a window into traces, the upper bound
'''
def plot_single_analysis(E_0, P, DELTA, eta, b=4, T=1000,
         llsim=0, ulsim=10000, llcyc=3000, ulcyc=5000):
    #constants!
    tau_p = 2*(10**-3) # tau_p is the photon lifetime
    
    #Get the system of equations
    funcs = setup(P, DELTA, b, eta, T)
    
    #Define initial values
    init=[np.sqrt(P),0,0]
    
    #simulate the system for bounds llsim and ulsim (lower and upper level)
    traces, time_trace = simulate_functions(init, funcs, llsim, ulsim)
    
    #normalize the electric field and perform a DFT on it
    norm_e = normalize_trace(traces[0], 100)
    fft_trace = freq_analysis_trace(norm_e, time_trace, 1.0)
    
    names = ['E','N','psi'] #names of each trace
    
    plot_n_traces(traces, time_trace, names)
    plot_limit_cycle(traces[0][llcyc:ulcyc], traces[1][llcyc:ulcyc],
                     names[0], names[1])
    
    #plot frequency analysis:
    f, ax = plt.subplots()
    f.suptitle("Frequency Spectrum")
    ax.plot(fft_trace)
    ax.set_xlabel("Freq (Hz?)")
    
    plt.show()

'''
    Sweep through eta values, recording the min and max points in the resulting
    waveforms, as well as the frequency spectrum
    
    parameters:
    e_min - starting eta value
    e_min - ending eta value
    e_step - the smallest change in eta
    P - the pumping rate
    DELTA - the detuning factor
    b (nom 4) - the linewidth enhancement factor
    T (nom 1000) - carrier : photon lifetime
    llsim - simulation starting point
    ulsim - simulation ending point
    llcyc - for creating a window into traces, the lower bound
    ulcyc - for creating a window into traces, the upper bound
    
    returns:
    e_values - the created axis (np array) of eta values
    bfdiag_points - the list of lists of peaks found for each eta
    freqs - the list (np array/matrix) of frequencies for each eta
'''
def eta_sweep(e_min, e_max, e_step, P, DELTA, b=4, T=1000, 
              reverse=False, llsim=0, ulsim=10000, sim_step=0.5,
              llcyc=5500, ulcyc=6000, continuation=False):
    #Define initial values
    init=[np.sqrt(P),0,0]
    
    #The Eta axis
    e_values = np.arange(e_min, e_max, e_step)
    if reverse: e_values = np.flip(e_values)
    e_size = e_values.size
    
    #The empty bfdiag:
    bfdiag_points = []
    
    #The empty frequency heatmap
    freqs = np.zeros((e_size, ulcyc-llcyc))
    
    f_out = None
    first_run = True
    for n in range(e_size):
        bfdiag_points.append([])
    for n, eta in enumerate(e_values):
        #if reverse: n = e_size - n - 1
        print('n is: {0}\r'.format(n), end='')
        #Get the system of equations
        funcs = setup(P, DELTA, b, eta, T)
        
        #Take a window into the relevant portion of the E-field trace
        y, t = simulate_functions(init, funcs, llsim, ulsim, step=sim_step)
        f_out = y[0][llcyc:ulcyc]
        if continuation: init = [i[-1] for i in y]
            #print('size of init: {0}\r'.format(len(init)))
        
        #Perform frequency analysis on that window
        #freqs[n] = freq_analysis_trace(f_out, times, 0)
        
        #Find all local minima and maxima and add them to the bifurcation diagrm
        bfdiag_points[n] += list(f_out[find_peaks(f_out)[0]][:4])
        bfdiag_points[n] += list(f_out[find_peaks(-f_out)[0]][:4])
    print('\nEta sweep complete')
    return e_values, bfdiag_points, freqs
    
def convert_bf_to_traces(bfdiag_points):
    #convert bfdiag scatter points to traces
    bfdiag_traces = []
    for n in range(len(bfdiag_points[0])):
        trace = []
        for m in range(len(bfdiag_points)):
            try:
                trace.append(bfdiag_points[m][n])
            except Exception as e:
                print(e)
                trace.append(0.0)
        bfdiag_traces.append(np.array(trace))
    return bfdiag_traces
    
'''
    Based on bifurcation analysis, plot all points
'''
def plot_bif_diag(e_values, bfdiag_points, override_fig=None, override_ax=None,
                  our_color='k'):
    if type(override_fig) is type(None) and type(override_ax) is type(None):
        fig, ax = plt.subplots()
        fig.suptitle("Bifurcation analysis of eta from {0} to {1}".format(
                        round(e_values[0], 4),
                        round(e_values[-1], 4)))
        ax.plot(e_values, np.ones_like(e_values), '--')
    else:
        fig = override_fig
        ax = override_ax
    
#    skipfrom = None
#    for n, eta in enumerate(e_values):
#        if eta >= 0.0093:
#            skipfrom=n
#            break
#    
#    for t in convert_bf_to_traces(bfdiag_points):
#        ax.plot(e_values[0:skipfrom], t[0:skipfrom], ':', c=our_color)
        
    for n, eta in enumerate(e_values):
        #if n < skipfrom: continue
        for pt in bfdiag_points[n]:
            try:
                ax.scatter(eta, pt, s=1, c=our_color)
            except Exception as e:
                print(e)
    ax.set_xlabel('Eta')
    ax.set_ylabel('Amplitude Extrema')
    return fig, ax

def plot_waterfall(e_values, freqs):
    fig, ax = plt.subplots()
    ax.set_title('Frequency Analisys')
    ax.pcolormesh(freqs.T, cmap='hot')
    ax.set_xticks(e_values)
    ax.set_xlabel('Eta')
    
def save_data(e_values, bfdiag_points, freqs):
    pass

if __name__ == '__main__':
    ev, bf, fr = eta_sweep(0.002, 0.0125, 0.000005, 0.375, 0, T=155)
    ev2, bf2, _ = eta_sweep(0.002, 0.0125, 0.000005, 0.375, 0, T=155,
                            continuation=True, our_color='b')
    
    save_data(ev, bf, fr)
    
    fig, ax = plot_bif_diag(ev, bf)
    plot_bif_diag(ev2, bf2, fig, ax)
    #plot_waterfall(ev, fr)
    plt.show()
    
    
