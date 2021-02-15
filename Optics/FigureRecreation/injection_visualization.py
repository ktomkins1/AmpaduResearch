import matplotlib.pyplot as plt
import matplotlib.figure as figmod
import numpy as np
import os

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

    returns plot.show
'''
def plot_limit_cycle(trace1, trace2, name1, name2):
    fig, axes = plt.subplots()
    fig.suptitle("Limit Cycle")
    axes.plot(trace1, trace2)
    axes.set_xlabel(name1)
    axes.set_ylabel(name2)
    return plt.show

'''
    Based on bifurcation analysis, plot all points

    param e_values: The x axis, eta points to use
    param bfdiag_points: The points which are collected as results
    param override_fig: The figure object to use if not a new one
    param override_ax: The axis object to use if not a new one
    param our_color: the color of the points in the plot
    param get_fig: Whether to return the figure and axis objects for later use
    param get_plot: Whether to return the runtime showing function from plt

    returns fig, ax if get_fig
    returns plt.show if get_plot and also not get_fig
'''
def plot_bif_diag(values, bfdiag_points, value_name, config, save_loc,
                  override_fig=None, override_ax=None,
                  our_color='k', get_fig=False):
    if type(override_fig) is type(None) and type(override_ax) is type(None):
        fig = figmod.Figure(figsize=(15.0, 10.0), dpi=500)
        ax = fig.add_axes([0.05, 0.075, 0.9, 0.85])
        fig.suptitle("Bifurcation analysis of {0} from {1} to {2}\n".format(
                        value_name,
                        round(values[0], 4),
                        round(values[-1], 4)) + get_param_description(config))

    else:
        fig = override_fig
        ax = override_ax
        #TODO: since plotted after a different bfd, add a key to the figure

    for n, v in enumerate(values):
        #if n < skipfrom: continue
        for pt in bfdiag_points[n]:
            if config['bf_absv']: pt = np.abs(pt)
            try:
                ax.scatter(v, pt, s=1, c=our_color)
            except Exception as e:
                print(e)
                print('exception while plotting scatter point')
    ax.set_xlabel(value_name)
    ax.set_ylabel('Amplitude Extrema')

    #try to plot lines for the hopfs
    try:
        eta_FH = config['eta_FH']
        if min(values) < eta_FH < max(values):
            ax.axvline(eta_FH)
        eta_RH = config['eta_RH']
        if min(values) < eta_RH < max(values):
            ax.axvline(eta_RH)
    except Exception as e:
        print('hopf rev or fwd not present')


    if config['vis_save']:
        fstring = 'result_figure{}.png'.format(config['bf_plot_id'])
        fname = os.path.join(save_loc, fstring)
        fig.savefig(fname, dpi=500)
    if get_fig:
        return fig, ax
    if config['vis_show']:
        return plt.show

def get_fit(axis, bfdiag_points):
    #take the first and last points on average
    fpm = np.mean(bfdiag_points[0])
    lpm = np.mean(bfdiag_points[-1])

    #get slope
    m = (lpm - fpm)/(axis[-1] - axis[0])

    #fpy = m*fpx + b, b = fpy - m*fpx
    b = fpm - m*axis[0]

    #get points for the line to plot
    pts = lambda x: m*x + b
    return pts(axis)


def plot_waterfall(e_values, freqs):
    fig, ax = plt.subplots()
    ax.set_title('Frequency Analisys')
    ax.pcolormesh(freqs.T, cmap='hot')
    ax.set_xticks(e_values)
    ax.set_xlabel('Eta')
    return plt.show

def get_param_description(c):
    desc = ''
    desc += 'Pumping:' + str(c['P'])
    desc += ' Linewidth:' + str(c['alpha'])
    desc += ' T:' + str(c['T'])
    desc += ' Detuning:' + str(c['DELTA'])
    return desc

if __name__ == '__main__':
    #TODO: open a pickled data_set to plot
    pass
