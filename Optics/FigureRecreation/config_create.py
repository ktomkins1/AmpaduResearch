import io
import json
import pickle
import inspect

#change the values below before running
config = {
    #parameters of the system
    'E_0':'lambda P: np.sqrt(P)',  #initial e-field
    'theta_0':0.0,      #initial phase difference
    'N_0':0.0,          #initial carrier density
    'alpha':4.0,        #the line-width enhancement factor    [typ 4]
    'LAMBDA':'sweep',   #the coupling constant*    [typ 10^-4 to 10^0]
    'P':0.375,          #the excess pumping rate
    'DELTA':0.0,        #the cavity optical detuning
    'T':1000.0,         #the ratio of carrier lifetime to photon lifetime    [typ 10^3]
    'tau_p':0.002,      #the photon lifetime    [typ 2x10^-3 ns]
    'tau_c':2.0,        #the carrier lifetime    [typ 2 ns]

    #model to use
    'model':'model_phys_review_1996',
    'model_shortname':'PhyRev1996',

    #parameters for sweeping
    'bf_reverse':False,        #run in reverse for bifurcation sweeps
    'bf_continuation':True,    #TODO: check terminology

    #parameters of the simulator
    'llsim':0,              #the beginning time point?
    'ulsim':6000,           #end the simulator after this many points
    'sim_step':1.0,         #the size of a point?
    'llcyc':5500,           #begin a window to test the steady state behavior
    'ulcyc':6000,           #end window
    'ex_bias':0.001,        #the sensitivity of grouping extremas

    #parameters for plotting
    'bf_absv':False,           #plot the absolute value of the bf points
    'bf_fit_line':True         #plot line of best fit
}

def create_short_desc(c, sep='-'):
    desc = c['model_shortname']
    desc += sep + 'lw' + str(c['alpha'])
    desc += sep + 'T' + str(c['T'])
    desc += sep + 'DT' + str(c['DELTA'])
    desc += sep + 'cc' + str(c['LAMBDA'])
    r, cnt = c['bf_reverse'], c['bf_continuation']
    if r or cnt:
        desc += sep + r*'r' + cnt*'c'
    return desc

def encode_config_hash(c):
    pass

if __name__ == '__main__':
    #TODO: offer interactive through command line

    #serialize and create encoding
    enc = '0000000'
    #create description string
    desc_short = create_short_desc(config)
    #add string and encoding to object
    config['enc'] = enc
    config['desc'] = desc_short

    #save to json
    with open("config.json", 'w') as fp:
        json.dump(config, fp, indent=4)
