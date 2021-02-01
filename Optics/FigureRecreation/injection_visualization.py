import matplotlib.pyplot as plt
import numpy as np

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
    return plt.show
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
    return plot.show
    
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
