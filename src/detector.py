# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 09:33:05 2023

@author: blr
"""
import numpy as np

##############################################################################
########################### MODEL FUNCTIONS ##################################
##############################################################################

def dGdt_from_d2YdtdE(t, E, x, m, d2YdtdE_arr, sens_func, params):
    sens_arr = sens_func(E, *params)
    dGdt_arr = simpson_int(d2YdtdE_arr*sens_arr, E)
    return dGdt_arr

def linear_sens(E, a, b):
    return E*a + b

# def V_from_dYdt()

##############################################################################
########################### HELPER FUNCTIONS #################################
##############################################################################

def simpson_int(arr, E):
    dE = (E[:,1]-E[:,0])
    integral = (3*dE/8)*(arr[:,0] +
                         3*np.sum(arr[:,1:-1:3], axis=1) +
                         3*np.sum(arr[:,2:-1:3], axis=1) +
                         2*np.sum(arr[:,3:-1:3], axis=1) +
                         arr[:,-1])
    return integral