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

def dispatch(setup, config):
    cc.fix_config(config)
    if config['bf']:
        #run sweep
        skey = [s for s in config.keys() if type(config[s]) is dict][0]
        rn = np.linspace
        if config[skey].keys()[0] == 'arange':
            rn = np.arange
        
        #TODO: use calculations to determine the critical points of interest
        results = sim.general_sweep(setup, config, skey, 
                                    config[skey].values()[0], rn)
    else:
        results = sim.trace(setup, config)
    
    #save images of results
    plot_image(results, config)

def plot_image(results, config): #for easy access based on a config
    pass

def run_and_compare_sweeps(): #for easy comparison on same plot
    #this will need to be based on settings because many/most times bifurcation
    #   diagrams won't meaningfully fit on the same plot
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

def enumerate_configs(c):
    c['bf'] = False
    elist, ekey, clist = None, None, None
    
    #check each item in config c and see if c is a sweep
    #check for enumerator
    for k in c.keys():
        if type(c[k]) == type([]):
            elist = c[k]
            ekey = k
        elif type(c[k]) == type({}):
            c['bf'] = True
    
    if elist == None:
        clist = [c]
    else:
        for item in elist:
            d = c.copy()
            d[ekey] = item
            clist.append(d)
    return clist, elist, ekey

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

    #create folder under results
    results_dirname = os.path.join('results', config['desc'])
    os.makedirs(results_dirname, exist_ok=True)

    clist = enumerate_configs(config)
    
    #dispatch all threads
    threads = [Thread(target=dispatch, args=(model.setup, c) for c in clist]
    for t in threads: t.start()
    for t in threads: t.join()

