# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 09:40:34 2023

@author: blr
"""
import astropy.units as u
import astropy.constants as c
import numpy as np
import numba
import time
import matplotlib.pyplot as plt
import scipy.interpolate as interp

# Unless otherwise specified, all units are assumed to be defined by
# c = 1, 1ns = 1, m_neutron = 1, and 1V = 1. These are referred to hereafter as 
# analysis units


##############################################################################
########################### MODEL FUNCTIONS ##################################
##############################################################################

def d2YdtdE_func_vacprop(t, E, x, m, d2YdtdE_0_func, params):
    """
    Describes propagation of massive particles described by a joint distribution
    over time and energy as they move through a vacuum.

    Parameters
    ----------
    t : TYPE
        DESCRIPTION.
    E : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.
    m : TYPE
        DESCRIPTION.
    d2YdtdE_0_func : function [f8(params), f8(params)]
        A function of type f(t, E, m, *params) that describes
    params : TYPE
        DESCRIPTION.

    Returns
    -------
    d2YdtdE_arr : TYPE
        DESCRIPTION.

    """
    tret = t - tof_func(E, m, x)
    d2YdtdE_arr = d2YdtdE_0_func(tret, E, *params)
    return d2YdtdE_arr

def dYdt_IRF_vacprop_gauss(t, x, m, sig_E, mu_E):
    E = E_func(t, m, x)
    dEdtof = dEdtof_func(E, m, x)
    IRF = dEdtof * np.exp(-0.5*((E-mu_E)/sig_E)**2) / np.sqrt(2*np.pi*sig_E**2)
    return IRF

def dYdt_IRF_vacprop_arbdYdE(t, x, m, E_arr, dYdE_arr):
    E_tof = E_func(t, m, x)
    dEdtof = dEdtof_func(E_tof, m, x)
    dYdE_tof = interp.interp1d(E_arr, dYdE_arr,bounds_error=False,fill_value=0.0)(E_tof)
    IRF = -dEdtof * dYdE_tof
    return IRF

def dYdt_func_vacprop_gaussdYdE(t, E, x, m, dYdt):
    """
    Describes propagation of massive particles described by a product of a distribution
    over time and a gaussian distribution over energy

    Parameters
    ----------
    t : TYPE
        DESCRIPTION.
    E : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.
    m : TYPE
        DESCRIPTION.
    d2YdtdE_0_func : function [f8(params), f8(params)]
        A function of type f(t, E, m, *params) that describes
    params : TYPE
        DESCRIPTION.

    Returns
    -------
    d2YdtdE_arr : TYPE
        DESCRIPTION.

    """
    tret = t - tof_func(E, m, x)
    d2YdtdE_arr = d2YdtdE_0_func(tret, E, *params)
    return d2YdtdE_arr



##############################################################################
########################### HELPER FUNCTIONS #################################
##############################################################################

# @numba.njit()
def E_func(tof, m, x):
    """
    Energy of a particle given its time of flight, mass, and distance traveled
    
    Parameters
    ----------
    tof : f8 value or array
        Particle time of flight in analysis units
    m : f8 value or array
        Particle mass in analysis units
    x : f8 value or array
        Particle distance traveled in analysis units

    Returns
    -------
    E : TYPE
        DESCRIPTION.

    """
    E = m * ((1 - (x/tof)**2)**(-1/2) - 1)
    return E

@numba.njit()
def v_func(E, m):
    if m != 0:
        gamma = E/m + 1
        v = np.sqrt(1-gamma**-2)
    return v

# @numba.njit()
def dEdtof_func(E, m, x):
    gamma = E/(m) + 1
    dEdt = -m *  (v_func(E,m)*gamma)**3 / x       
    return dEdt

# @numba.njit()
def tof_func(E, m, x):
    t = x/v_func(E,m)
    return t

# @numba.njit()
def get_t_and_E_arr(dt, nE, t0_min, t0_max, E0_min, E0_max, m, x):
    tmin = tof_func(E0_max, m, x)+t0_min
    tmax = tof_func(E0_min, m, x)+t0_max
    nt = int(((tmax-tmin)/dt))
    t_arr = np.linspace(tmin, tmax, nt)[:,np.newaxis]
    if tmin <= t0_max or np.any(t_arr - (t0_max +t0_min)/2 <= 0):
        E_arr = np.tile(np.linspace(E0_min, E0_max, nE)[np.newaxis, :], (nt,1))
        
    else:
        Emid = E_func(t_arr - (t0_max +t0_min)/2, m, x)
        delta_E =  (t0_max - t0_min) * -dEdtof_func(Emid, m, x)
        delta_E = np.minimum(delta_E, E0_max-E0_min)
        # Emid = np.maximum(delta_E/2+1*u.keV, Emid)
        E_arr = np.tile(np.linspace(-0.5,0.5,nE)[np.newaxis,:], (nt,1)) * delta_E + Emid
    
    return t_arr, E_arr

