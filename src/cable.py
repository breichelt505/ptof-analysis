# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 09:40:58 2023

@author: blr
"""
import numpy as np
import scipy.ndimage
import scipy.optimize as opt
import os
import pathlib
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
import astropy.constants as c

tstar = 1*u.ns
vstar = c.c


lstar = (tstar*vstar).si

##############################################################################
########################### MODEL FUNCTIONS ##################################
##############################################################################

def find_H_linear_phase_model(cable, length, f_arr):
    """
    Find transfer function of cable given a length, input time array, and type
    of cable. Assumes a perfectly linear model for phase (ie all component
    frequencies move at the same speed)

    Parameters
    ----------
    cable : str
        Cable name
    length : float
        Cable length in analysis units (ns*c = 0.2997 m)
    f_arr : ndarray of float
        Array of frequencies at which the transfer function is to be evaluated
        in analysis units (1/ns = GHz)

    Returns
    -------
    H : ndarray of float
        Transfer function of the cable-- the inverse fourier transform of this
        function provides the impulse response of the cable

    """
    file = pathlib.Path(f"{os.getcwd()}\\tables\\cable_response\\coax_info.xlsx")
    data = pd.read_excel(file,header=0,index_col=0)
    
    k1 = data["k1"][cable]
    k2 = data["k2"][cable]
    v = data["velocity factor"][cable]

    H = find_H_from_ks_linear_phase(k1, k2, v, length, f_arr)
    
    return H

def find_H_nonlinear_phase_model(cable, length, f_arr):
    file = pathlib.Path(f"{os.getcwd()}\\tables\\cable_response\\coax_info.xlsx")
    data = pd.read_excel(file,header=0,index_col=0)
    
    # k1 = data["k1"][cable]
    k2 = data["k2"][cable]
    v = data["velocity factor"][cable]
    ri = data["ri (in)"][cable]
    ro = data["ro (in)"][cable]
    roo = data["roo (in)"][cable]
    g_c = data["g_c (S*m)"][cable]
    
    eps_r_d = 1/v**2
    losstan_d = v*k2/2780
    
    Z0, gamma, H = find_H_from_cable_params(
        float(ri*u.imperial.inch/u.m),
        float(ro*u.imperial.inch/u.m),
        float(roo*u.imperial.inch/u.m), 
        g_c, 
        float(length*lstar/u.m), 
        eps_r_d, 
        losstan_d, 
        f_arr*1e9)
    
    return H


##############################################################################
########################### HELPER FUNCTIONS #################################
##############################################################################

def find_k1k2_from_s12(cable):
    file = pathlib.Path(f"{os.getcwd()}\\tables\\cable_response\\{cable}.csv")
    data = np.loadtxt(file,delimiter=",")
    fGHz = data[:,0]
    atten_dB100ft = data[:,1]
    i_sort = np.argsort(fGHz)
    fGHz = fGHz[i_sort]
    atten_dB100ft = atten_dB100ft[i_sort]
    
    def fit_func(f, k1, k2):
        return k1*np.sqrt(f) + k2*f
    
    popt, pcov = opt.curve_fit(fit_func, fGHz, atten_dB100ft)
    
    plt.scatter(fGHz,atten_dB100ft)
    plt.plot(np.linspace(0,50,1000),fit_func(np.linspace(0,50,1000), *popt))
    
    return popt

def find_H_from_ks_linear_phase(k1, k2, v, length, f_arr):
    lstar = c.c*u.ns
    atten = (k1*np.sqrt(f_arr) + k2*f_arr)*(length * lstar / (100*u.imperial.ft)).si.value
    H = 10**(-atten/20) * np.exp(-1j*2*np.pi*f_arr*length/v)
    
    return H

def find_H_from_cable_params(ri_si, ro_si, roo_si, g_c_si, x_si, eps_r_d, losstan_d, f_arr_GHz):
    eps0 = float(c.eps0/(u.F/u.m))
    mu0 = float(c.mu0.si/(u.H/u.m))
    omega = f_arr_GHz*2*np.pi
    
    Y = 2*np.pi*(eps_r_d*losstan_d*eps0*omega + 1j*omega*eps_r_d*eps0) / np.log(ro_si/ri_si)
    Rs = np.sqrt(mu0*omega / (2*g_c_si))
    Rc_hf = (Rs/(2*np.pi)) * (ro_si**-1 + ri_si**-1)
    Rc_lf = (g_c_si*np.pi)**-1 * (ri_si**-2 + (roo_si**2 - ro_si**2)**-1)
    R = np.maximum(Rc_hf, Rc_lf)
    
    L_hf = (mu0/(2*np.pi))*np.log(ro_si/ri_si) + (Rs/(2*np.pi*omega)) * (ro_si**-1 + ri_si**-1)
    Zapprox = R + 1j*omega*L_hf

    Zapprox[np.isnan(Zapprox)] = Zapprox[(~np.isnan(Zapprox)).nonzero()[0][0]]
    
    gamma = np.sqrt(Zapprox*Y)
    H = np.exp(-gamma*x_si)
    Z0 = np.sqrt(Zapprox/Y)
    
    return Z0, gamma, H
