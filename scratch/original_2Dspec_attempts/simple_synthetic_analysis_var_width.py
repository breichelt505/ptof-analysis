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
    
    # V = interp.interp1d(t_arr, V)(t_in)
    
    return t_arr, V

def f(t_in, N, mu_t, sig_t, sig_E):
    t_sim, V_sim = output(N, mu_t, sig_t, 0.0025, sig_E)
    return interp.interp1d(t_sim.flatten(), V_sim.flatten(), 
                           bounds_error=False, fill_value=0.0)(t_in)

t_exp, V_exp = output(0.87, 0.0, 0.1, 0.0025, 0.0001)
V_exp += (np.random.randn(len(t_exp)) * V_exp.max() * 0.01)
plt.plot(t_exp, V_exp)
sigma = np.ones(len(t_exp)) * V_exp.max() * 0.01


lb = [1e-6,-np.inf,1e-6,1e-6]
ub = [np.inf]*4
p0 = [1.0,1.0,1.0,0.00012]

popt, pcov = opt.curve_fit(f, t_exp.flatten(), V_exp.flatten(), 
                           p0 = p0, #sigma = sigma,
                           #bounds = (lb,ub))
                           )
t_bestfit, V_bestfit = output(popt[0],popt[1],popt[2],0.0025,popt[3])
plt.plot(t_bestfit, V_bestfit)