"""
Mean field model of microcircuit
"""

import numpy as np

def make_weights():
    """
    Creates network weights constrained by 
    mouse literature
    Arguments:
        (None): 

    Returns:
        params (dict): with keys 
            weights (5,5): Network weights, ordered as
            (soma, dend, pv, sst, vip)
            tau_e: glut time const
            tau_i: gaba time const

    """
    # 0: E, 1: D, 2: P, 3: S, 4: V
    wee = 0.42
    wed, wep = 1, 0.6*0.7
    wde, wds = 0.1*0.42, 0.25*1.96
    wpe, wpp, wps = 0.45 * 1, 0.5 * 1.5, 0.6 *1.3
    wse, wsv = 0.35*1, 0.7 * 0.25
    wve, wvs = 1, 0.7*0.25
    # Combine them
    Weights = np.zeros((5,5))
    Weights[0,0] = wee
    Weights[0,[1,2]] = wed, wep
    Weights[1,[0,3]] = wde, wds
    Weights[2, [0,2,3]] = wpe, wpp, wps
    Weights[3,[0, 4]] = wse, wsv
    Weights[4,[0, 3]]= wve, wvs
    Weights[:, 2:] *= -1 # make inhibition inhibitory
    params = {'weights': Weights, 'tau_e': 10., 'tau_i': 2.}
    return params


def find_rates(weights, inputs):
    """
    Creates network weights constrained by 
    mouse literature. Inverse of find_inputs
    Arguments:
        weights (n, n): recurrent connectivity
        inputs (n, ): external inputs 

    Returns:
        (5,): steady state rate for each cell/compartment
        (soma, dend, pv, sst, vip)

    """
    assert weights.shape[0] == weights.shape[1]
    n = weights.shape[0]
    return np.linalg.inv((np.eye(n) - weights))@inputs

def find_inputs(weights, target_rates):
    """
    Finds external inputs that will lead to 
    given steady state rates. Inverse of find_rates

    Arguments:
        weights (n, n): recurrent connectivity
        target_rates (n, ): rates we want to achieve

    Returns:
        (5,): external inputs to each cell/compartment
        (soma, dend, pv, sst, vip)

    """
    assert weights.shape[0] == weights.shape[1]
    assert weights.shape[0] == target_rates.shape[0]
    n = weights.shape[0]
    return (np.eye(n) - weights)@target_rates

def state_modulation(weights, x_base, x_mod):
    r0 = find_rates(weights, x_base)
    rm = find_rates(weights, x_base + x_mod)
    return (rm - r0) / (rm + r0)

# Simulation functions
def relu(x):
    """
    Transfer function
    """
    x[x < 0] = 0
    return x 

def soma_dynamics(params, re, rd, rp, Ie):
    wed, wep =  params['weights'][0,1], params['weights'][0,2]
    tau_e = params['tau_e']
    return (-re + wed*rd + wep*rp + Ie)/tau_e

def dend_dynamics(params, rd, re, rs, Id):
    wde, wds =  params['weights'][1,0], params['weights'][1,3]
    tau_e = params['tau_e']
    return (-rd + wde*re + wds*rs + Id)/tau_e

def pv_dynamics(params, rp, re, rs, Ip):
    wpe, wpp, wps = params['weights'][2,0], params['weights'][2,2], params['weights'][2,3]
    tau_i = params['tau_i']
    return (-rp + wpe*re + wpp*rp + wps*rs + Ip)/tau_i

def sst_dynamics(params, rs, re, rv, Is):
    wse, wsv = params['weights'][3,0], params['weights'][3,4]
    tau_i = params['tau_i']
    return (-rs + wse*re + wsv*rv + Is)/tau_i

def vip_dynamics(params, rv, re, rs, Iv):
    wve, wvs = params['weights'][4,0], params['weights'][4,3]
    tau_i = params['tau_i']
    return (-rv + wve*re + wvs*rs + Iv)/tau_i 

def update(params, r, inputs, dt = 0.1):
    """
    Simulates network dynamics one time step forward

    Arguments:
        params: dict with "weights" and "tau_e", "tau_i"
        r (n, ): rates
        inputs (n,): external inputs
        dt (float): time step

    Returns:
        r_new (n, ): updated rates

    """
    dredt = soma_dynamics(params, r[0], r[1], r[2], inputs[0])
    drddt = dend_dynamics(params, r[1], r[0], r[3], inputs[1])
    drpdt = pv_dynamics(params, r[2], r[0], r[3], inputs[2])
    drsdt = sst_dynamics(params, r[3], r[0], r[4], inputs[3])
    drvdt = vip_dynamics(params, r[4], r[0], r[3], inputs[4])
    r_new = r + dt * np.array([dredt, drddt, drpdt, drsdt, drvdt])
    r_new = relu(r_new)  # rectify
    return r_new

def simulate(params, time_steps, x_baseline, x_mod, stim_period, 
        dt = .1, noise_sd = 0):
    """
    Simulates network 

    Arguments:
        params: dict with "weights" and "tau_e", "tau_i"
        time_steps (int): in ms
        x_baseline (n, ): baseline inputs to each neuron/compartment
        x_mod (n, ): additional input
        stim_period (2, ): start 

    Returns:
        (time_steps, n): rates for each timestep and cell/compartment
    """
    assert params['weights'].shape[0] == params['weights'].shape[1]
    assert params['weights'].shape[0] == x_baseline.shape[0]
    assert x_baseline.shape[0] == x_mod.shape[0]
    assert stim_period[1] >= stim_period[0]
    assert stim_period[1] < time_steps
    n = params['weights'].shape[0]
    rates = np.zeros((time_steps, n), dtype=float)
    for t in range(time_steps-1):
        # Determine inputs
        if (t >= stim_period[0]) and (t <= stim_period[1]):
            inputs = x_baseline + x_mod
        else:
            inputs = x_baseline
        inputs += np.random.normal(size=(n,))*noise_sd * np.sqrt(dt)
        # simulate forward
        rates[t+1] = update(params, rates[t], inputs, dt)
    return rates


def simulate_variable_input(params,  
        x_baseline, x_mod, x_var, dt = .1):
    """
    Simulates network 

    Arguments:
        params: dict with "weights" and "tau_e", "tau_i"
        time_steps (int): in ms
        x_baseline (n, ): baseline inputs to each neuron/compartment
        x_mod (n, ): additional input
        x_var (timesteps, 2): time varying input to soma and dend
        stim_period (2, ): start 

    Returns:
        (time_steps, n): rates for each timestep and cell/compartment
    """
    assert params['weights'].shape[0] == params['weights'].shape[1]
    assert params['weights'].shape[0] == x_baseline.shape[0]
    assert x_baseline.shape[0] == x_mod.shape[0]
    n = params['weights'].shape[0]
    time_steps = x_var.shape[0]
    rates = np.zeros((time_steps, n), dtype=float)
    inputs = x_baseline + x_mod + x_var
    for t in range(time_steps-1):
        # simulate forward
        rates[t+1] = update(params, rates[t], inputs[t], dt)
    return rates

def simulate_both(params,  stim_period, 
        x_baseline, x_mod, x_var, dt = .1):
    """
    Simulates network 

    Arguments:
        params: dict with "weights" and "tau_e", "tau_i"
        time_steps (int): in ms
        x_baseline (n, ): baseline inputs to each neuron/compartment
        x_mod (n, ): additional input
        x_var (timesteps, 2): time varying input to soma and dend
        stim_period (2, ): start 

    Returns:
        (time_steps, n): rates for each timestep and cell/compartment
    """
    assert params['weights'].shape[0] == params['weights'].shape[1]
    assert params['weights'].shape[0] == x_baseline.shape[0]
    assert x_baseline.shape[0] == x_mod.shape[0]
    n = params['weights'].shape[0]
    time_steps = x_var.shape[0]
    rates = np.zeros((time_steps, n), dtype=float)
    inputs = x_baseline + x_mod + x_var
    for t in range(time_steps-1):
        if (t >= stim_period[0]) and (t <= stim_period[1]):# or (t >= stim_period[2]) and (t <= stim_period[3]):
            inputs = x_baseline + x_var[t] + x_mod
        else:
            inputs = x_baseline + x_var[t]
        # simulate forward
        rates[t+1] = update(params, rates[t], inputs, dt)
    return rates


