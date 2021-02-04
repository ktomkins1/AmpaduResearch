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
    filepath = os.path.join(rdir,'bf_results_{}'.format(config['desc'])+'.pickle')
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

    #if encoding is empty, create new encoding
    if 'enc' not in config.keys():
        config['enc'] = cc.encode_config_hash(config)
        config['desc'] = cc.create_short_desc(config)

    #create folder under results
    results_dirname = os.path.join('results', config['desc'])
    os.makedirs(results_dirname, exist_ok=True)
    
    #call simulation(s)
    eta_start = 0.0
    eta_end = 0.014
    d_eta = 0.0005
    vals, bf, fs = sim.eta_sweep(model.setup, eta_start, eta_end, d_eta, 
                                 config['P'], config['DELTA'], T=config['T'],
                                 llsim=config['llsim'], ulsim=config['ulsim'],
                                 llcyc=config['llcyc'], ulcyc=config['ulcyc'])

    #call visualization(s)
    #show = vis.plot_n_traces([ys[0], ys[2]], t, ["E1","N1"], "Electric Field for eta={0}, sim@6000,5700:5999".format(eta))
    show = vis.plot_bif_diag(vals, bf, 'Eta', config)
    show()

    #save data
    save_bifurcation_data(results_dirname, config, vals, bf)

    #save configuration
    with open(results_dirname + '/config_' + config['desc'] + '.json', 'w+') as fp:
        json.dump(config, fp, indent=4)
