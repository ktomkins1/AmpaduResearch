import injection_simulation as sim
import injection_visualization as vis
import numpy as np
import pickle #for archiving 
import json   #for configuration

#import model here

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
    
    #call simulation(s)
    t, ys, fs = sim.trace_single_analysis([np.sqrt(P),0,0], P, DLT, eta, T=155, llsim=0, ulsim=6000, llcyc=5500, ulcyc=5999, sim_step=1.0)
    
    #save data
    
    #call visualization(s)
    show = vis.plot_n_traces([ys[0], ys[2]], t, ["E1","N1"], "Electric Field for eta={0}, sim@6000,5700:5999".format(eta))
    #show()
