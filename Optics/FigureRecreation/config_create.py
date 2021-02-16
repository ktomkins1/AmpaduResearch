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
    'T':958.0,                  #carrier lifetime to photon lifetime  [typ 10^3]
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
    'bf_absv':False,            #plot the absolute value of the bf points
    'bf_fit_line':False,        #plot line of best fit
    'vis_show':False,           #show the plot using plt runtime. DNU with MT
    'vis_save':True,            #save the plot as a picture
    'bf_plot_num': 1,           #how many plots to make for bf diagrams
    'bf_plot_id': 1             #which results of a multi-result plot is this?
}

optional_params = ['desc', 'enc', 'root_dir', 'gamma_r', 'omega_r',
                   'eta_FH', 'eta_RH', 'bf_plot_id']
required_params = ['E_0','theta_0','N_0','alpha','eta','P','DELTA','T',
                   'tau_p','tau_c','model','model_shortname','bf_reverse',
                   'bf_continuation','llsim','ulsim','sim_step','llcyc',
                   'ulcyc','ex_bias','bf_absv','bf_fit_line',
                   'vis_save','vis_show', 'bf_plot_num']
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
    for k in required_params:
        if k not in c.keys():
            print('Configuration Malformed. missing: {}'.format(k))
    c['enc'] = encode_config_hash(c)
    if c['bf_plot_num'] > 1 : return
    c['bf_plot_id'] = 0
    c['desc'] = create_short_desc(c)

def create_config(E_0='np.sqrt(c[\'P\'])', theta_0=0.0, N_0=0.0, alpha=4.8,
                  eta=0.01, P=1.0, DELTA=0.0, T=958.0,
                  tau_p=0.002,  tau_c=2.0, model='model_phys_review_1996',
                  model_shortname='PhyRev1996', bf_reverse=False,
                  bf_continuation=True, llsim=0, ulsim=6000, sim_step=1.0,
                  llcyc=5500, ulcyc=6000, ex_bias=0.001, bf_absv=False,
                  bf_fit_line=True, vis_save=True, vis_show=False,
                  bf_plot_num=1):
    config = locals()
    fix_config(config)
    return config

def save_config(c, name='config.json'):
    with open(name, 'w') as fp:
        json.dump(c, fp, indent=4)

if __name__ == '__main__':
    #TODO: offer interactive through command line

    fix_config(config)

    #save to json
    save_config(config)
