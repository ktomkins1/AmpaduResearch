import numpy as np
from numpy.lib.scimath import sqrt
from scipy.integrate import ode, solve_ivp
from maxima_determination import get_extrema

'''
    One single frame of integration
    
    based on previous values, calculate next values of the system
    
    parameters:
    t - the single time value
    y - the vectorized current state of the system [s1, s2, ... sn]
    funcs - the python functions corresponding to each diff eq F1..n
'''
def int_step(t,y,funcs):
    #return np.array([f(*y) for f in funcs])
    return np.array([f(*y) for f in funcs])

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
def simulate_functions(initial_values, funcs, ll, ul, step=1.0, tol=None):
    #Calculate the result
    res_map = solve_ivp(int_step, [ll, ul], initial_values,
                        args=(funcs,), max_step=step)#, rtol=tol)
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
def trace_single_analysis(setup, init, P, DELTA, eta, b=4, T=1000, sim_step=1.0,
                          llsim=0, ulsim=10000, llcyc=3000, ulcyc=5000):
    #Get the system of equations
    funcs = setup(P, DELTA, b, eta, T)
    
    #simulate the system for bounds llsim and ulsim (lower and upper level)
    traces, time_trace = simulate_functions(init, funcs, llsim, ulsim, step=sim_step)
    
    #normalize the electric field and perform a DFT on it
    norm_e = normalize_trace(traces[0], 100)
    fft_trace = freq_analysis_trace(norm_e, time_trace, 1.0)
    
    return time_trace, traces, fft_trace

#function eventually for getting a single trace
def trace(setup, config):
    pass

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
def eta_sweep(setup, e_min, e_max, e_step, P, DELTA, alpha=4, T=1000,       #TODO: take in axis, sweep other params
              reverse=False, llsim=0, ulsim=6000, sim_step=1.0,         #TODO: dynamic determination of start and end
              llcyc=5500, ulcyc=6000, continuation=True, ex_bias=0.001):
    #Define initial values
    init=[np.sqrt(P),0,0]
    
    #The Eta axis
    e_values = np.arange(e_min, e_max, e_step)
    if reverse: e_values = np.flip(e_values)
    e_size = e_values.size
    print('Performing {0} simulations...'.format(e_size))
    
    #The empty bfdiag:
    bfdiag_points = []
    
    #The empty frequency heatmap
    freqs = np.zeros((e_size, ulcyc-llcyc))
    
    f_out = None
    first_run = True
    for n in range(e_size):
        bfdiag_points.append([])
    for n, eta in enumerate(e_values):
        LAMBDA=eta
        print('n is: {0}\r'.format(n), end='') #TODO: time remaining calc
        #Get the system of equations
        funcs = setup(P, DELTA, alpha, LAMBDA, T)
        
        #Take a window into the relevant portion of the E-field trace
        y, t = simulate_functions(init, funcs, llsim, ulsim, step=sim_step)
        f_out = y[0][llcyc:ulcyc]
        if continuation: init = [i[-1] for i in y]
        
        #Perform frequency analysis on that window
        #freqs[n] = freq_analysis_trace(f_out, times, 0)
        
        #Find all local minima and maxima and add them to the bifurcation diagrm
        bfdiag_points[n] = get_extrema(f_out, ex_bias)
    print('\nEta sweep complete')
    return e_values, bfdiag_points, freqs
 
#def convert_bf_to_traces(bfdiag_points):  TODO: create a fn for more efficient plotting
#    #convert bfdiag scatter points to traces
#    bfdiag_traces = []
#    for n in range(len(bfdiag_points[0])):
#        trace = []
#        for m in range(len(bfdiag_points)):
#            try:
#                trace.append(bfdiag_points[m][n])
#            except Exception as e:
#                print(e)
#                trace.append(0.0)
#        bfdiag_traces.append(np.array(trace))
#    return bfdiag_traces

'''
    Perform a sweep of any value
    
    param setup: the callback to the setup which creates our model's fns
    param c:     the config dictionary
    param sweep_key: the key in the dictionary which is the value being swept
    param sweep_space: aruments to create the axis
    param axis_gen: the callback to function to create the axis
    
    returns: the axis created for sweeping
    returns: a list of lists of points for each value on the axis
    returns: a numpy array of frequencies which occurred at each point
'''
def general_sweep(setup, c, sweep_key, sweep_space, axis_gen=np.linspace):
    #Define initial values
    init=[c['E_0'],c['theta_0'],c['N_0']] #must be prepared
    
    #The Eta axis
    print('Sweep space is {}'.format(sweep_space))
    sweep_values = axis_gen(*sweep_space)
    if c['bf_reverse']: sweep_values = np.flip(sweep_values)
    sweep_size = sweep_values.size
    print('Performing {0} simulations...'.format(sweep_size))
    
    #The empty bfdiag:
    bfdiag_points = []
    
    #The empty frequency heatmap
    freqs = np.zeros((sweep_size, c['ulcyc']-c['llcyc']))
    
    f_out = None
    for n in range(sweep_size):
        bfdiag_points.append([])
    for n, val in enumerate(sweep_values):
        c[sweep_key] = val
        #print('n is: {0}\r'.format(n), end='') #TODO: time remaining calc
        #Get the system of equations
        funcs = setup(c['P'], c['DELTA'], c['alpha'], c['eta'], c['T'])
        
        #Take a window into the relevant portion of the E-field trace
        y, t = simulate_functions(init, funcs, c['llsim'], c['ulsim'], 
                                  step=c['sim_step'])
        f_out = y[0][c['llcyc']:c['ulcyc']]
        if c['bf_continuation']: init = [i[-1] for i in y]
        
        #Perform frequency analysis on that window
        #freqs[n] = freq_analysis_trace(f_out, times, 0)
        
        #Find all local minima and maxima and add them to the bifurcation diagrm
        bfdiag_points[n] = get_extrema(f_out, c['ex_bias'])
    print('Sweep of var {0} complete'.format(sweep_key))
    return sweep_values, bfdiag_points, freqs
    
