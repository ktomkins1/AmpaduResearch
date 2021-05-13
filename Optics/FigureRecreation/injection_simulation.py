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
    return [f(t,*y) for f in funcs]

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
    ftrace = np.abs(np.fft.fftshift(np.fft.fft(trace)))
    freqs = np.fft.fftshift(np.fft.fftfreq(len(time_trace), dt))
    return ftrace, freqs

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

''' function eventually for getting a single trace '''
def trace(setup, config):
    pass


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
        if n > 0 and n % 20 == 0:
            print('{}: n is {} out of {}'.format(c['bf_plot_id'],n,sweep_size))
        #Get the system of equations
        funcs = setup(c['P'], c['DELTA'], c['alpha'], c['eta'], c['T'])

        #Take a window into the relevant portion of the E-field trace
        y, t = simulate_functions(init, funcs, c['llsim'], c['ulsim'],
                                  step=c['sim_step'])
        f_out = y[0][c['llcyc']:c['ulcyc']]
        if c['bf_continuation']: init = [abs(i[-1]) for i in y] #y is trace list

        #Perform frequency analysis on that window
        #TODO: figure out dt
        freqs[n] = freq_analysis_trace(f_out, t[c['llcyc']:c['ulcyc']], 1.0)[0]

        #Find all local minima and maxima and add them to the bifurcation diagrm
        bfdiag_points[n] = get_extrema(f_out, c['ex_bias'])
    print('Sweep of var {0} complete'.format(sweep_key))
    return sweep_values, bfdiag_points, freqs

'''
    An alternative to np.linspace or arange for generating an axis

    param mag: the base for creating the axis is 1-10^-mag
'''
def get_exponential_axis(start, stop, num, mag=2):
    print('exp v1')
    r = 1 - np.power(10.0, -mag)
    ax = 1 - np.power(r, np.arange(num))
    ax *= (stop + start)
    ax -= start
    return ax
    
'''
    An alternative to np.linspace or arange for generating an axis
    
    param start_val: the value at 0%
    param stop_val: the value at 100%
    param start_pct: the percent to start at
    param stop_pct: the percent to stop at
    
'''
def get_pct_bounds_axis(start_val, stop_val, start_pct, stop_pct, num):
    p100 = stop_val - start_val
    lbound = p100*start_pct/100 + start_val
    rbound = p100*stop_pct/100 + start_val
    return np.linspace(lbound, rbound, num)
    
    
    
    
    
