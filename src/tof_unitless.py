# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 09:09:06 2023

@author: blr
"""
# All units are assumed to be such that m_neutron=1, c=1, e=1, and 1 ns = 1

import astropy.units as u
import astropy.constants as c
import numpy as np
import scipy.integrate as sciint
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as interp
# import numba
import time

def E(tof, m, x):
    E = m * ((1 - (x/tof)**2)**(-1/2) - 1)
    return E

def v(E, m):
    gamma = E/m + 1
    v = np.sqrt(1-gamma**-2)
    return v

def dEdtof(E, m, x):
    gamma = E/(m) + 1
    dEdt = -m *  (v(E,m)*gamma)**3 / x       
    return dEdt

def tof(E, m, x):
    t = x/v(E,m)
    return t.to(u.ns)

def d2YdtdE(t, E, x, m, d2YdtdE_0_func, params):
    
    tret = t - tof(E, m, x)
    
    d2YdtdE_arr = d2YdtdE_0_func(tret, E, params)
    return d2YdtdE_arr

def get_t_and_E_arr(dt, nE, t0_min, t0_max, E0_min, E0_max, m, x):
    tmin = tof(E0_max, m, x)+t0_min
    tmax = tof(E0_min, m, x)+t0_max
    nt = int(((tmax-tmin)/dt))
    t_arr = np.linspace(tmin, tmax, nt)[:,np.newaxis]
    if tmin <= t0_max or np.any(t_arr - (t0_max +t0_min)/2 <= 0*u.s):
        E_arr = np.tile(np.linspace(E0_min, E0_max, nE)[np.newaxis, :], (nt,1))
        
    else:
        Emid = E(t_arr - (t0_max +t0_min)/2, m, x)
        delta_E =  (t0_max - t0_min) * -dEdtof(Emid, m, x)
        delta_E = np.minimum(delta_E, E0_max-E0_min)
        Emid = np.maximum(delta_E/2+1*u.keV, Emid)
        E_arr = np.tile(np.linspace(-0.5,0.5,nE)[np.newaxis,:], (nt,1)) * delta_E + Emid
    
    return t_arr, E_arr
    

def simpson_int(arr, E):
    dE = (E[:,1]-E[:,0])
    integral = (3*dE/8)*(arr[:,0] +
                         3*np.sum(arr[:,1:-1:3], axis=1) +
                         3*np.sum(arr[:,2:-1:3], axis=1) +
                         2*np.sum(arr[:,3:-1:3], axis=1) +
                         arr[:,-1])
    return integral

def dGdt(t, E, x, m, d2YdtdE_0_func, sens_func, params):
    d2YdtdE_arr = d2YdtdE(t, E, x, m, d2YdtdE_0_func, params)
    sens_arr = sens_func(E)
    dGdt_arr = simpson_int(d2YdtdE_arr*sens_arr, E)
    return dGdt_arr



Emin = 1.5*u.MeV
Emax = 3.5*u.MeV
dt = 1*u.ps
nE = 10
m = c.m_p
t0 = 10.0*u.ns
t_model = np.linspace(0,0.3, 10)*u.ns
Ti_model = np.r_[5.0, 5.0, 6.0, 7.0, 8.0, 9.0, 3.0, 15.0, 2.0, 5.0]*1.0*u.keV
dYdt_model = np.r_[0.0, 0.1, 0.3, 0.7, 1.0, 2.0, 3.0, 1.2, 0.2, 0.0]
params = {"omega0" : 82.542*u.keV**0.5,
          "E0" : (2.4495+0.01)*u.MeV,
          "alpha_d_omega" : {"1" : 1.7013e-3, "2" : 0.16888, "3" : 0.49, "4" : 7.9469e-4},
          "alpha_d_eth" : {"1" : 4.69515, "2" : -0.040729, "3" : 0.47, "4" : 0.81844},
          "t0" : t0,
          "t_model" : t_model,
          "Ti_model" : Ti_model,
          "dYdt_model" : dYdt_model
          }
          

tic = time.time()
for x in np.linspace(7,20,5)*u.cm:        
    # plt.figure()
    t_arr, E_arr = get_t_and_E_arr(dt, nE, t0+0.0*u.ns, t0+0.3*u.ns, Emin, Emax, m, x)
    # plt.pcolor(np.tile(t_arr, (1,E_arr.shape[1])).si.value, E_arr.si.value, d2YdtdE(t_arr, E_arr, x, m, d2YdtdE_0_ballabio_gaussian, params).si.value)
    # plt.colorbar()
    
    dGdt_arr  = dGdt(t_arr, E_arr, x, m, d2YdtdE_0_ballabio_gaussian, sens_protons, params)
    # plt.figure()
    # plt.plot(t_arr.flatten()-t_arr.flatten()[dGdt_arr==dGdt_arr.max()], dGdt_arr/dGdt_arr.max())
    # plt.plot(t_arr.flatten(), dGdt_arr/dGdt_arr.max())
toc = time.time()    

print(toc-tic)

