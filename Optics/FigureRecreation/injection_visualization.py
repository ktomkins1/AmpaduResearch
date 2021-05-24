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


#TODO: Modularity...
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
def plot_bif_diag(results, value_name, config, save_loc, override_fig=None, override_ax=None,
                  our_color='k', get_fig=False):
    xaxis = results['axis']
    bfdiag_points = results['groups']
    bifurcations = []
    try:
        bifurcations = results['bifurcations']
    except KeyError:
        pass

    vsize, hsize = 9, 12

    if type(override_fig) is type(None) and type(override_ax) is type(None):
        try:
            vsize, hsize = config['vis_v'], config['vis_h']
        except KeyError:
            vsize, hsize = 9, 12
        fig = figmod.Figure(figsize=(hsize, vsize))
        ax = fig.add_axes([0.05, 0.075, 0.9, 0.85])
        fig.suptitle("Bifurcation analysis of {} from {} to {} (n={})\n".format(
                        value_name,
                        round(xaxis[0], 4),
                        round(xaxis[-1], 4),
                        len(xaxis)) + get_param_description(config))

    #TODO: find these optimally
    our_dpi = 500
    xpix, ypix = hsize*our_dpi, vsize*our_dpi


    zo = 0
    
    try:
        if config['vis_fimage']:
            #try to get frequency results out
            freqs_raw = np.abs(results['fr'])
            fr_axis = results['fx']
            hm = np.log10(freqs_raw.T + 0.001)

            #get the bounds to make it the right size
            xlen = len(xaxis)
            ylen = config['ulcyc'] - config['llcyc']
            if xlen < ylen: xlen = ylen

            #get extrema range
            bf_range, _, _ = convert_bf_to_array(bfdiag_points)
            #ymin = min(bf_range)
            #ymax = max(bf_range)
            ymin, ymax = ax.get_ybound()
            yspan = ymax-ymin
            ypt = yspan/(ymax*100)

            #sweep range
            xmin = min(xaxis)
            xmax = max(xaxis)
            xspan = xmax-xmin
            xpt = xspan/(xmax*100)

            #ax.set_xlim(xmin=xmin, xmax=xmax)
            #ax.set_ylim(ymin=ymin, ymax=ymax)

            #plot it as an image in the zo^th z plane
            zo += 1
            ax.imshow(hm, extent=(xmin-0.5*xpt, xmax-0.5*xpt, ymax-0.5*ypt, ymin-0.5*ypt), zorder=1)
            print('xlim: ' + str(ax.get_xlim()))
            print('ylim: ' + str(ax.get_ylim()))
            print('xbound: ' + str(ax.get_xbound()))
            print('ybound: ' + str(ax.get_ybound()))

            #set the axis titles for the frequency yaxis, the extrema y axis and the sweep axis
            #frequency
            transform = lambda x: np.round(x - ylen/2)
            ax_y2 = ax.secondary_yaxis('right', functions=(transform, transform))
            #extrema
            yticks = list(range(0, ylen, round(0.15*ylen)))
            ylist = np.round(np.linspace(ymin,ymax, ylen),2)
            ax.set_yticks(yticks)
            ax.set_yticklabels(ylist[yticks])
            #sweep
            xticks = list(range(0, xlen, round(0.2*xlen)))
            ax.set_xticks(xticks)
            ax.set_xticklabels(list(np.round(np.linspace(xaxis[0], xaxis[-1], len(xticks)),4)))
            ax.tick_params('x', labelrotation=-90.0)
            
    except KeyError as ke:
        print('did not show waterfall: ' + str(ke))
        
    zo += 1
    for n, v in enumerate(xaxis):
        for pt in bfdiag_points[n]:
            if config['bf_absv']: pt = np.abs(pt)
            ax.scatter(v, pt, zorder=2, s=1, c=our_color)#, marker=',')
    
    ax.set_xbound(lower=min(xaxis), upper=max(xaxis))
    vmode = 0
    try:
        vmode = config['vbounds']['mode']
    except KeyError:
        pass
    if vmode == 1:
        ax.set_ybound(lower=config['vbounds']['al'], 
                      upper=config['vbounds']['au'])
    if vmode == 2:
        #offset mode. offset has already been modified based on normalization
        try:
            c = config['vbounds']['o']
        except KeyError:
            c = 0
        d = config['vbounds']['+-']
        ax.set_ybound(lower=c-d, upper=c+d)
    
    #print('xlim: ' + str(ax.get_xlim()))
    #print('ylim: ' + str(ax.get_ylim()))
    #print('xbound: ' + str(ax.get_xbound()))
    #print('ybound: ' + str(ax.get_ybound()))

    ax.set_xlabel(value_name)
    ax.set_ylabel('Amplitude Extrema')

    try:
        if config['vis_showbfs']:
            for bx, by in bifurcations:
                ax.scatter(bx, by, color='r')
    except KeyError:
        pass

    #try to plot lines for the hopfs
    try:
        if value_name == 'eta' and config['vis_vlines']:
            eta_FH = float(config['eta_FH'])
            if min(xaxis) < eta_FH < max(xaxis):
                ax.axvline(eta_FH, color='r')
            eta_RH = float(config['eta_RH'])
            if min(xaxis) < eta_RH < max(xaxis):
                ax.axvline(eta_RH, color='r')
    except KeyError as ke:
        print('hopf rev or fwd not present: ' + str(ke))


    if config['vis_save']:
        all_pics_folder = os.path.join(config['root_dir'], 'all_pics')
        os.makedirs(all_pics_folder, exist_ok=True)
        local_fstring = 'result_figure{}.png'.format(config['bf_plot_id'])
        global_fstring = config['desc'] + '_' + str(config['bf_plot_id']) + '.png'
        for fname in [os.path.join(save_loc, local_fstring),
                      os.path.join(all_pics_folder, global_fstring)]:
            fig.savefig(fname, dpi=our_dpi)
    if get_fig:
        return fig, ax
    if config['vis_show']:
        return plt.show

def convert_bf_to_array(bf, rounding=3, xprop=4, yprop=3):

    #find largest and smallest points
    newarr = []
    for gr in bf:
        newarr += [round(i, rounding) for i in gr]

    mini = min(newarr)
    maxi = max(newarr)

    #determing minimum distance between points
    newarr = sorted(newarr)
    min_dist = np.inf
    for i in range(len(newarr)+1):
        if i >= len(newarr): break
        if i == 0: continue
        dist = newarr[i] - newarr[i-1]
        if 0 < dist < min_dist:
            min_dist = dist

    bf_range = np.arange(mini, maxi, min_dist)

    #find scaling factor(s) to elongate image to fit
    xlen = len(bf)
    ylen = len(bf_range)
    #should be 4:3
    newxlen = (xprop*ylen)/yprop
    newxlen = int(np.ceil(newxlen))
    xscale_factor = int(np.ceil(newxlen/xlen))

    our_image = np.zeros((len(bf_range), newxlen))
    for i, gr in enumerate(bf):
        for pt in gr:
            rinx = np.searchsorted(bf_range, pt)
            if rinx > 0:
                rinx -= 1
            for j in range(xscale_factor):
                dinx = i*xscale_factor + j
                if dinx >= newxlen: continue
                our_image[rinx, dinx] += 1.0

    return bf_range, our_image, newxlen

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


def plot_waterfall(axis, freqs, config, skey):
    fig, ax = plt.subplots()
    ax.set_title('Frequency Analisys')
    #ax.pcolormesh(np.log10(freqs.T), cmap='hot')
    ax.imshow(np.log10(freqs.T), cmap='hot')
    ax.set_xticks(axis)
    ax.set_xlabel(skey)
    save_loc = config['targetdir']
    all_pics_folder = os.path.join(config['root_dir'], 'freq_pics')
    os.makedirs(all_pics_folder, exist_ok=True)
    local_fstring = 'freq_figure{}.png'.format(config['bf_plot_id'])
    global_fstring = config['desc'] + '_' + str(config['bf_plot_id']) + '.png'
    for fname in [os.path.join(save_loc, local_fstring),
                  os.path.join(all_pics_folder, global_fstring)]:
        fig.savefig(fname)
    #return plt.show
    return

def get_param_description(c):
    desc = ''
    subdesc = ['P:',' alpha:',' T:',' Detuning:',' eta:']
    descvals = [c['P'],c['alpha'],c['T'],c['DELTA'],c['eta']]

    for i, sstr in enumerate(subdesc):
        if type(descvals[i]) != dict:
            desc += sstr + str(descvals[i])

    return desc

if __name__ == '__main__':
    #TODO: open a pickled data_set to plot
    pass
