# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 19:02:14 2023

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

Emin = 0.01
Emax = 0.09
dt = 0.01
nE = 300
m = 1.0
mu_t = 0.1
sig_t = 0.1
mu_E = 0.05
sig_E = 0.01
t0 = 0.1
N = 1.0

          
# plt.figure()
tic = time.time()
for x in np.linspace(0.0,30,2):        
    # plt.figure()
    t_arr, E_arr = propagation.get_t_and_E_arr(dt, nE, t0-0.5, t0+0.5, Emin, Emax, m, x)
    d2YdtdE = propagation.d2YdtdE_func_vacprop(
        t_arr, E_arr, x, m, 
        source.d2YdtdE_0_func_doublegauss,
        np.r_[N, mu_t, sig_t, mu_E, sig_E]
        )
    dGdt = detector.dGdt_from_d2YdtdE(
        t_arr, E_arr, x, m, d2YdtdE, 
        detector.quad_sens,
        np.r_[0.0,1.0,0.0]
        )
    # plt.pcolormesh(np.tile(t_arr, (1,E_arr.shape[1])), E_arr, d2YdtdE)
    # plt.colorbar()
    plt.plot(t_arr,dGdt/dGdt.max())
    
toc = time.time()    

print(toc-tic)