# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 09:27:10 2023

@author: blr
"""
import numpy as np

##############################################################################
########################### MODEL FUNCTIONS ##################################
##############################################################################

def d2YdtdE_0_func_doublegauss(t, E, N, mu_t, sig_t, mu_E, sig_E):
    d2YdtdE_0 = N/2/np.pi/sig_t/sig_E * np.exp(-0.5 * (((t-mu_t)/sig_t)**2 + ((E-mu_E)/sig_E)**2))
    return d2YdtdE_0

def dYdt_0_func_gauss(t, E, N, mu_t, sig_t):
    dYdt_0 = N/np.sqrt(2*np.pi)/sig_t * np.exp(-0.5 * ((t-mu_t)/sig_t)**2)
    return dYdt_0