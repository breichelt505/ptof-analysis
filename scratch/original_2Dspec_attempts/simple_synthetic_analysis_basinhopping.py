# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 23:21:49 2023

@author: blr
"""

import astropy.units as u
import astropy.constants as c
import numpy as np
import numba
import time
import matplotlib.pyplot as plt
import source
import propagation
import utils
import detector
import scipy.stats as stats
import scipy.optimize as opt
import scipy.interpolate as interp

def output(N, mu_t, sig_t, mu_E, sig_E):
    dt = 0.01
    nE = 300
    m = 1.0
    x = 1.0
    E_min = max(mu_E-4*sig_E, 1e-4)
    
    t_arr, E_arr = propagation.get_t_and_E_arr(dt, nE, mu_t-5*sig_t, mu_t+5*sig_t, E_min, mu_E+5*sig_E, m, x)
    
    # Start with a simple double gaussian source
    d2YdtdE = propagation.d2YdtdE_func_vacprop(
        t_arr, E_arr, x, m, 
        source.d2YdtdE_0_func_doublegauss,
        np.r_[N, mu_t, sig_t, mu_E, sig_E]
        )
    
    # Assume simplest sensitivity-- constant
    dGdt = detector.dGdt_from_d2YdtdE(
        t_arr, E_arr, x, m, d2YdtdE, 
        detector.linear_sens,
        np.r_[0.0,1.0]
        )
    
    V = dGdt
    
    return t_arr, V

t_exp, V_exp = output(1.0, 0.0, 0.1, 0.0025, 0.0001)
sig = 0.01*V_exp.max() 
V_exp += np.random.randn(len(t_exp)) * sig
# plt.plot(t_exp, V_exp)

def rchisq(x):
    t_guess, V_guess = output(x[0],x[1],x[2], 0.0025, x[3])
    V_interp = interp.interp1d(
        t_guess.flatten(), V_guess.flatten(),bounds_error=False,
        fill_value = 0.0)(t_exp.flatten())
    rchisq = np.sum((V_interp-V_exp)**2/sig**2) / (len(t_exp)-5)
    return rchisq

sol = opt.minimize(rchisq, [1.0, 0.01, 0.11, 0.00, 0.0002],
                   bounds = [(0.1,10.0),(-1.0,1.0),(0.01,1.0),(0.001,0.005),(0.00005,0.0005)])