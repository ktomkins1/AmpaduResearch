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
    'P':9.3,   #the excess pumping rate
    'DELTA':0.0,                #the cavity optical detuning
    'T':958.0,                  #carrier lifetime to photon lifetime  [typ 10^3]
    'tau_p':0.002,              #the photon lifetime    [typ 2x10^-3 ns]
    'tau_c':2.0,                #the carrier lifetime    [typ 2 ns]

    #mode of the simulation
    'mode': 'bif',

    #model to use
    'model':'model_phys_review_1996',
    'model_shortname':'PhyRev1996',

    #parameters for sweeping
    'bf_reverse':False,        #run in reverse for bifurcation sweeps
    'bf_continuation':True,    #use the final values as next initial values
    'bf_cnb':0.1,

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
    'vis_type':'scatter',
    'bf_plot_num': 1,           #how many plots to make for bf diagrams
    'bf_plot_id': 1             #which results of a multi-result plot is this?
}

optional_params = ['desc', 'enc', 'root_dir', 'gamma_r', 'omega_r', 'mode',
                   'eta_FH', 'eta_RH', 'bf_plot_id', 'vis_type', 'bf_cnb']
required_params = ['E_0','theta_0','N_0','alpha','eta','P','DELTA','T',
                   'tau_p','tau_c','model','model_shortname','bf_reverse',
                   'bf_continuation','llsim','ulsim','sim_step','llcyc',
                   'ulcyc','ex_bias','bf_absv','bf_fit_line',
                   'vis_save','vis_show', 'bf_plot_num']
known_str_params = ['desc', 'enc', 'root_dir', 'model','mode',
                    'model_shortname', 'vis_type']
known_modes = ['single', 'bif', 'multi', 'stability']

def create_short_desc(c, sep='-'):
    desc = c['model_shortname']
    try:
        desc += c['mode'] + sep
    except KeyError:
        pass

    keys = ['alpha', 'T', 'DELTA', 'eta', 'P']
    for k in keys:
        if type(c[k]) == type({}):
            desc += k[0] + 'sweep'
        else:
            desc += k[0] + str(c[k])

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
    set_mode(c)
    c['enc'] = encode_config_hash(c)
    if c['bf_plot_num'] > 1 : return  #a (perhaps hackey) way of making many plots have same description
    c['bf_plot_id'] = 0
    c['desc'] = create_short_desc(c)

def set_mode(c):
    listcnt, dictcnt = 0,0
    for k in c.keys():
        if type(c[k]) == type({}): dictcnt += 1
        if type(c[k]) == type([]): listcnt += 1

    if dictcnt == 0: c['mode'] = 'single'
    if dictcnt == 1: c['mode'] = 'bif'
    if dictcnt > 1: c['mode'] = 'stability'
    if listcnt >= 1: c['mode'] = 'many_' + c['mode']


def create_config(E_0='np.sqrt(c[\'P\'])', theta_0=0.0, N_0=0.0, alpha=4.8,
                  eta=0.01, P=1.0, DELTA=0.0, T=958.0,
                  tau_p=0.002,  tau_c=2.0, model='model_phys_review_1996',
                  model_shortname='PhyRev1996', bf_reverse=False,
                  bf_continuation=True, llsim=0, ulsim=6000, sim_step=1.0,
                  llcyc=5500, ulcyc=6000, ex_bias=0.001, bf_absv=False,
                  bf_fit_line=True, vis_save=True, vis_show=False,
                  bf_plot_num=1, vis_type='scatter', bf_cnb=False):
    config = locals()
    fix_config(config)
    return config

def save_config(c, name='configs/config.json'):
    with open(name, 'w+') as fp:
        json.dump(c, fp, indent=4)

if __name__ == '__main__':
    #TODO: offer interactive through command line

    fix_config(config)

    #save to json
    save_config(config)
