import injection_simulation as sim
import injection_visualization as vis
import config_create as cc
import numpy as np
import pickle #for archiving
import json   #for configuration
#import importlib #for importing the model dynamically
import io
import os
import sys
import argparse
from matplotlib import pyplot as plt

#import model here
#TODO: dynamically import from config
import model_phys_review_1996 as model

def run_and_plot_sweep(): #for easy access based on a config
    pass

def run_and_compare_sweeps(): #for easy comparison on same plot
    pass

def save_single_sim_traces(config, time_axis, data):
    pass

def save_bifurcation_data(rdir, config, axis, data):
    our_obj = {'axis':axis, 'data':data}
    filepath = os.path.join(rdir,'bf_results_{}'.format(config['enc'])+'.pickle')
    with open(filepath, 'wb+') as fp:
        print('Saving data to disk....', end='')
        pickle.dump(our_obj, fp)
        print('Done')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', default='config.json',
                        help="The configuration file to load", dest='cfilename')

    args = parser.parse_args()

    #import config file
    '''
      Configs hold all of the parameters of the system.
      There can be parameters which are varied:
        one parameter can be given as a tuple: (start, end, number_of_steps)
        this is what is used as the sweep parameter for a bifurcation diagram
        
        one parameter can be given as a [list] of several values
        one bifurcation diagram or plot will be created for each.  these will
        then be plotted together (and plots of each alone will also be saved)
        
        if no list or range parameters are given, a plot of all of the traces
        will be created and stored.
    '''
    config = None
    with open(args.cfilename, 'r') as fp:
        config = json.load(fp)

    #ensure encoding and description are correct
    config['enc'] = cc.encode_config_hash(config)
    config['desc'] = cc.create_short_desc(config)

    #create folder under results
    results_dirname = os.path.join('results', config['desc'])
    os.makedirs(results_dirname, exist_ok=True)
    
    #find values to sweep and to display simulataneously
    bf_active = False #is this a single run or bifurcation sweep?
    comp_active = False #is there a value which we are comparing the keys
    skey, ckey = None, None
    bf_axis_arg = None
    comp_list = None
    for k in config.keys():
        tp = type(config[k])
        if tp == type({}):
            skey = k
            bf_active = True
            bf_axis_arg = config[k]['sweep']
        if tp == type([]):
            ckey = k
            comp_active = True
            comp_list = c[k]
    
    #call simulation(s)
    #TODO: This will need to be cleaned up
    results = []
    if comp_active:
        for cval in comp_list:
            config[ckey] = cval
            if bf_active:
                results.append(sim.general_sweep(model.setup, config, skey, bf_axis_arg))
            else:
                results.append(sim.single_sim(model.setup, config))
    elif bf_active:
        results.append(sim.general_sweep(model.setup, config, skey, bf_axis_arg))
    else:
        results.append(sim.single_sim(model.setup, config))

    fig, ax = plt.subplots()
    fig.suptitle('Variable sweep of {0}.  Other params: T={1}, alpha={2}, Zero Detuning'.format(skey, config['T'], config['alpha'])
    for n, color in enumerate(['k', 'r', 'b', 'g']):
        
    #call visualization(s)
    #show = vis.plot_bif_diag(vals, bf, 'Eta', config)
    #show()

    #save data
    if bf_active:
        for vals, bf, _ in results:
            save_bifurcation_data(results_dirname, config, vals, bf)

    #save configuration
    with open(results_dirname + '/config_' + config['enc'] + '.json', 'w+') as fp:
        json.dump(config, fp, indent=4)
