# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 09:27:10 2023

@author: blr
"""
import numpy as np
import scipy

##############################################################################
########################### MODEL FUNCTIONS ##################################
##############################################################################

def d2YdtdE_0_func_doublegauss(t, E, N, mu_t, sig_t, mu_E, sig_E):
    d2YdtdE_0 = N/2/np.pi/sig_t/sig_E * np.exp(-0.5 * (((t-mu_t)/sig_t)**2 + ((E-mu_E)/sig_E)**2))
    return d2YdtdE_0

def dYdt_0_func_gauss(t, E, N, mu_t, sig_t):
    dYdt_0 = N/np.sqrt(2*np.pi)/sig_t * np.exp(-0.5 * ((t-mu_t)/sig_t)**2)
    return dYdt_0

def ballabio_spectrum(reaction, Ti, Y, num_energies):
"""return a fusion product energy spectrum from plasma w/ ion temperature Ti [keV] and yield Y"""
    if reaction == "DDn":
        a1E = 4.69515
        a2E = -0.040729
        a3E = 0.47
        a4E = 0.81844
        a1w = 1.7013 * 1e-3
        a2w = 0.16888
        a3w = 0.49
        a4w = 7.9460 * 1e-4
        E0 = 2.4495 # MeV
        w0 = 82.542 # kev^(1/2)
    
        # interpolation formulas for deviation from standard results
        deltaE = (a1E/(1 + a2E * (Ti**a3E))) * (Ti**(2./3)) + a4E*Ti # keV
        Emean = E0 + deltaE*1e-3 #MeV
        deltaw = (a1w/(1 + a2w * (Ti**a3w))) * (Ti**(2./3)) + a4w*Ti
        FWHM = w0 * (1 + deltaw) * np.sqrt(Ti)
        Estdev  = FWHM / (2 * np.sqrt(2 * np.log(2))) * 1e-3 # MeV
    
    elif reaction == "DTn":
        a1E = 5.30509
        a2E = 2.4736 * 1e-3
        a3E = 1.84
        a4E = 1.3818
        a1w = 5.1068 * 1e-4
        a2w = 7.6223 * 1e-3
        a3w = 1.78
        a4w = 8.7691 * 1e-5
        E0 = 14.021 #Mev
        w0 = 177.259 #kev^(1/2)
        
        # interpolation formulas for deviation from standard results
        deltaE = (a1E/(1 + a2E * (Ti**a3E))) * (Ti**(2./3)) + a4E*Ti # keV
        Emean = E0 + deltaE*1e-3 #MeV
        deltaw = (a1w/(1 + a2w * (Ti**a3w))) * (Ti**(2./3)) + a4w*Ti
        FWHM = w0 * (1 + deltaw) * np.sqrt(Ti)
        Estdev  = FWHM / (2 * np.sqrt(2 * np.log(2))) * 1e-3 # MeV
    
    # full formula for D3He not given by Ballabio
    # so more approximate version is used, which should still be accurate to < 5% for Ti < 40keV
    elif reaction == "D3Hep": 
        E0 = 14.630 # Mev
        w0 = 180.985 # kev^(1/2)
        deltaE = (9 * Ti **(2./3)+Ti)*1e-3 # MeV
        Estdev = ((w0*(Ti)^(1/2))/(2*np.sqrt(2*np.log(2))))*1e-3 # MeV
        Emean = E0 + deltaE # MeV

    energies = np.linspace(Emean - 5*Estdev, Emean + 5*Estdev, num_energies)
    spectrum = scipy.stats.norm.pdf(energies, Emean, Estdev) * Y
    return energies, spectrum

##############################################################################
########################### HELPER FUNCTIONS #################################
##############################################################################