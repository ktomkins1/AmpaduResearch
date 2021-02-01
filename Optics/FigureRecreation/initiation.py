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

def save_bifurcation_data(config, axis, data):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', default='config.json',
                        help="The configuration file to load", dest='cfilename')

    args = parser.parse_args()

    #import config file
    config = None
    with open(args.cfilename, 'r') as fp:
        config = json.load(fp)

    #if encoding is empty, create new encoding
    if 'enc' not in config.keys():
        config['enc'] = cc.encode_config_hash(config)
        config['desc'] = cc.create_short_desc(config)

    #create folder under results
    results_dirname = 'results/'+config['desc']
    os.makedirs(results_dirname, exist_ok=True)

    #call simulation(s)
    vals, bf, fs = sim.eta_sweep(model.setup, eta_start, eta_end, d_eta, P, DLT, T=155, llsim=0, ulsim=6000, llcyc=5500, ulcyc=5999)

    #call visualization(s)
    #show = vis.plot_n_traces([ys[0], ys[2]], t, ["E1","N1"], "Electric Field for eta={0}, sim@6000,5700:5999".format(eta))
    show = vis.plot_bif_diag(vals, bf, 'Eta', absv=True)
    show()

    #save data
    save_bifurcation_data(results_dirname, config, vals, bf)

    #save configuration
    with open(results_dirname + '/config_' + config['enc'] + '.json', 'w') as fp:
        json.dump(config, fp)
