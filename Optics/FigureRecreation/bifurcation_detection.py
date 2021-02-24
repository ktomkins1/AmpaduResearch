import numpy as np
import injection_simulation as sim
import initiation as init

def get_FRDPs(c):
    if init.detect_sweep_key(c) in ['alpha', 'P', 'T']: return None, None
    alpha, P, T = c['alpha'], c['P'], c['T']
    gamma_r = (1 + 2*P)/(2*T)
    omega_r = np.sqrt(complex((2*P/T) - gamma_r**2))
    c['gamma_r'] = gamma_r
    c['omega_r'] = omega_r
    return gamma_r, omega_r

def get_fwd_rev_hopf(c, gamma_r=None, omega_r=None):
    if init.detect_sweep_key(c) in ['alpha', 'P', 'T']: return None, None
    alpha, P, T = c['alpha'], c['P'], c['T']

    if gamma_r == None: gamma_r = c['gamma_r']
    if omega_r == None: omega_r = c['omega_r']

    eta_FH = 2*gamma_r*(np.sqrt(alpha**2 + 1)/(alpha**2 - 1))
    eta_RH = omega_r*np.sqrt(complex((alpha**2 - 1)/2))
    c['eta_FH'] = eta_FH
    c['eta_RH'] = eta_RH
    return eta_FH, eta_RH

def get_cnb_from_groups(results, sens=0.01):
    #cnb = character and bifurcations
    bifs = []
    axis = results['axis']
    ourts = Traces(sens, axis[1]-axis[0])
    for i, g in enumerate(results['groups']):
        # for each point in g, track it's trace and watch if there are bifurcations
        for pt in g: ourts.append_trace_near(axis[i], pt)
    results['character'], results['bifurcations'] = ourts.report(axis)

class Traces:
    def __init__(self, sensitivity, dx):
        self.s = sensitivity
        self.dx = dx
        self.latest_id = 0
        self.trace_list = []
        self.bifurcations = []

    def add_trace(self, xval, yval):
        new_trace = {
            'root':xval,
            'id':self.next_id(),
            'values':[yval,],
            'latest_x':xval
            }
        self.trace_list.append(new_trace)

    def append_trace_near(self, xval, yval):
        if self.get_active_num(xval) == 0:
            self.add_trace(xval, yval)
        else:
            for t in self.get_active(xval):
                #check if it is within y bounds
                if abs(t['values'][-1] - yval) < self.s:
                    #check if this trace has seen a point already
                    if t['latest_x'] == xval:
                        #this is a lower bifurcation
                        self.add_trace(xval, yval)
                        self.add_bif(xval, (yval + t['values'][-1])/2)
                    else:
                        #this point is part of the trace t
                        t['latest_x'] = xval
                        t['values'].append(yval)
                    return

    def next_id(self):
        nid = self.latest_id
        self.latest_id += 1
        return nid

    def get_active(self, xval):
        ret_list = []
        for t in self.trace_list:
            if (xval >= t['root']) and (xval - t['latest_x'] <= self.dx):
                ret_list.append(t)
        return ret_list

    def get_active_num(self, xval):
        return len(self.get_active(xval))

    def character(self, xaxis):
        return [self.get_active_num(x) for x in xaxis]

    def add_bif(self, xval, yval):
        self.bifurcations.append((xval, yval))

    def report(self, xaxis):
        return self.character(xaxis), self.bifurcations
