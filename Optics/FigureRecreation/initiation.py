import injection_simulation as sim
import injection_visualization as vis
import numpy as np
import pickle #for archiving 
import json   #for configuration
import importlib #for importing the model dynamically
import io

#import model here
import model_antennas_2021 as model

def save_single_sim_traces(config, time_axis, data):
    pass

def save_bifurcation_data(config, axis, data):
    pass

if __name__ == '__main__':
    #import config file
    
    #define variables and parameters
    P = 0.375
    DLT = 0
    eta = 0.0105
    eta_start = 0.0
    eta_end = 0.1
    d_eta = 0.0001
    
    #call simulation(s)
    evals, bf, fs = sim.eta_sweep(model.setup, eta_start, eta_end, d_eta, P, DLT, b=2.3, T=155, llsim=0, ulsim=6000, llcyc=5500, ulcyc=5999, sim_step=1.0, continuation=False)
    
    #save data
    
    #call visualization(s)
    #show = vis.plot_n_traces([ys[0], ys[2]], t, ["E1","N1"], "Electric Field for eta={0}, sim@6000,5700:5999".format(eta))
    show = vis.plot_bif_diag(evals, bf)
    show()
