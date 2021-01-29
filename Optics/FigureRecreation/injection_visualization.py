import matplotlib.pyplot as plt

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
def plot_single_analysis(P, DELTA, eta, b=4, T=1000,
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
    
    plot_n_traces(traces, time_trace, names, 'System for Eta of {0}'.format(eta))
    plot_limit_cycle(traces[0][llcyc:ulcyc], traces[1][llcyc:ulcyc],
                     names[0], names[1])
    
    #plot frequency analysis:
    f, ax = plt.subplots()
    f.suptitle("Frequency Spectrum")
    ax.plot(fft_trace)
    ax.set_xlabel("Freq (Hz?)")
    
    plt.show()
    
    
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
