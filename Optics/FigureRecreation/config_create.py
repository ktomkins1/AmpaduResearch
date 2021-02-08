import io
import json
import pickle
import inspect
import lzma
from hashlib import blake2b

#change the values below before running
#one dict and one list are allowed.  
#These will determine values to sweep and compare, respectively
#dict option must be like {'linspace'|'arange':(start, end, n|d)}
#list option should be less than 8 options
#string option will be executed by exec and must be a statement which will be
#   prepended with config[<key name>]=
#string options should use c[<key>] if wanting to refer to variable <key>
#   for example, if refering to P in order to create E_0, use c[P]
config = {
    #parameters of the system
    'E_0':'np.sqrt(c[\'P\'])',  #initial e-field
    'theta_0':0.0,              #initial phase difference
    'N_0':0.0,                  #initial carrier density
    'alpha':4.8,                #the line-width enhancement factor    [typ 4]
    'eta':{'linspace':(0.0, 0.014, 10)},   #the coupling constant
    'P':[1.0, 2.0, 5.0, 9.0],   #the excess pumping rate
    'DELTA':0.0,                #the cavity optical detuning
    'T':958.0,                 #carrier lifetime to photon lifetime  [typ 10^3]
    'tau_p':0.002,              #the photon lifetime    [typ 2x10^-3 ns]
    'tau_c':2.0,                #the carrier lifetime    [typ 2 ns]

    #model to use
    'model':'model_phys_review_1996',
    'model_shortname':'PhyRev1996',

    #parameters for sweeping
    'bf_reverse':False,        #run in reverse for bifurcation sweeps
    'bf_continuation':True,    #use the final values as next initial values

    #parameters of the simulator
    'llsim':0,              #the beginning time point?
    'ulsim':6000,           #end the simulator after this many points
    'sim_step':1.0,         #the size of a point?
    'llcyc':5500,           #begin a window to test the steady state behavior
    'ulcyc':6000,           #end window
    'ex_bias':0.001,        #the sensitivity of grouping extremas

    #parameters for plotting
    'bf_absv':False,           #plot the absolute value of the bf points
    'bf_fit_line':True,         #plot line of best fit
    'vis_show':False,
    'vis_save':True
}

optional_params = ['desc', 'enc', 'root_dir']
required_params = ['E_0','theta_0','N_0','alpha','eta','P','DELTA','T',
                   'tau_p','tau_c','model','model_shortname','bf_reverse',
                   'bf_continuation','llsim','ulsim','sim_step','llcyc',
                   'ulcyc','ex_bias','bf_absv','bf_fit_line','vis_show']
known_str_params = ['desc', 'enc', 'root_dir', 'model', 'model_shortname']

def create_short_desc(c, sep='-'):
    desc = c['model_shortname']
    desc += sep + 'lw' + str(c['alpha'])
    desc += sep + 'T'  + str(c['T'])
    desc += sep + 'DT' + str(c['DELTA'])
    desc += sep + 'cc' + str(c['eta'])
    desc += sep + 'P'  + str(c['P'])
    r, cnt = c['bf_reverse'], c['bf_continuation']
    if r or cnt:
        desc += sep + r*'r' + cnt*'c'
    return desc.replace(' ', '').replace('.', '_')

def encode_config_hash(c):
    h = blake2b(digest_size=10)
    h.update(lzma.compress(pickle.dumps(c)))
    return h.hexdigest()
#    return bytes.hex(lzma.compress(pickle.dumps(c)))

def fix_config(c):
    c['desc'] = create_short_desc(c)
    c['enc'] = encode_config_hash(c)

def create_config(E_0='lambda P: np.sqrt(P)', theta_0=0.0, N_0=0.0, alpha=4.0,
                  eta='sweep', P=0.375, DELTA=0.0, T=1000.0,
                  tau_p=0.002,  tau_c=2.0, model='model_phys_review_1996',
                  model_shortname='PhyRev1996', bf_reverse=False,
                  bf_continuation=True, llsim=0, ulsim=6000, sim_step=1.0,
                  llcyc=5500, ulcyc=6000, ex_bias=0.001, bf_absv=False,
                  bf_fit_line=True):
    config = locals()
    fix_config(config)
    return config

if __name__ == '__main__':
    #TODO: offer interactive through command line

    #serialize and create encoding
    enc = encode_config_hash(config)[10:20] #hopefully different enough
    #create description string
    desc_short = create_short_desc(config)
    #add string and encoding to object
    config['enc'] = enc
    config['desc'] = desc_short

    #save to json
    with open("config.json", 'w') as fp:
        json.dump(config, fp, indent=4)

