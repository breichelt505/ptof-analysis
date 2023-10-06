# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 23:41:53 2023

@author: blr
"""

import numpy as np
import matplotlib.pyplot as plt
import src.cable
import astropy.units as u
import astropy.constants as c
import os
import pathlib
import pandas as pd
import scipy.ndimage

# file_path = r"C:\Users\blr\Dropbox (MIT)\MIT\Research\Projects\ptof_analysis\tables\detector_response\megaray_05-22\SC1000A_avg_xrays_nofilter.txt"
file_path = r"C:\Users\blr\Dropbox (MIT)\MIT\Research\Projects\ptof_analysis\tables\detector_response\megaray_05-22\SC1000A_avg_electrons_on_cap.txt"
megaray_IRF = np.loadtxt(file_path)

ARC_IRF = pd.read_excel(r"C:\Users\blr\Dropbox (MIT)\PTOF_GUI_Ben\Codes\Detector configs\IRFs\90-078 IRFs\SC1000 PTOF -500V IRF N221115-timing N221115-IRF-SC1000A.xlsx",header=0,sheet_name="IRF")

# t = np.linspace(-1000,1000,100000)
# t = t_avg*1e9
# n = 1000000
# tmax = 1000
# t_arr = np.linspace(0.0, tmax, n)

t_arr = megaray_IRF[:-1,0]*1e9
V = megaray_IRF[:-1,1]

lstar = c.c*u.ns

n = len(t_arr)
dt = t_arr[1]-t_arr[0]

fGHz = np.fft.rfftfreq(n,dt)

H_snout = src.cable.find_H_nonlinear_phase_model("RG402",(3100*u.mm/lstar).si.value,fGHz)
H_DLP = src.cable.find_H_nonlinear_phase_model("UFB293C",(150*u.imperial.inch/lstar).si.value,fGHz)
H_track = src.cable.find_H_nonlinear_phase_model("UFB293C",(25*u.imperial.ft/lstar).si.value,fGHz)
H_infra = src.cable.find_H_nonlinear_phase_model("LMR600",(35*u.m/lstar).si.value,fGHz)
H_rack1 = src.cable.find_H_nonlinear_phase_model("LL290HF",(120*u.imperial.inch/lstar).si.value,fGHz)
H_rack2 = src.cable.find_H_nonlinear_phase_model("LL290HF",(72*u.imperial.inch/lstar).si.value,fGHz)

H_pre_bias = H_snout*H_DLP*H_track*H_infra
H_tot = H_snout*H_DLP*H_track*H_infra*H_rack1*H_rack2

# plt.plot(t_arr,np.fft.irfft(H_pre_bias))
# plt.plot(t_arr,np.fft.irfft(H_tot))
t_arr2 = np.array(ARC_IRF["T"])
V_sig2 = np.array(ARC_IRF["V"])

plt.figure()
V_sig = np.fft.irfft(H_tot*np.fft.rfft(V))
V_sig2_plot = V_sig2[(V_sig2>=np.nanmax(V_sig2*0.5)).nonzero()[0][0]]
plt.plot(t_arr-t_arr[(V_sig*np.nanmax(V_sig2)/np.nanmax(V_sig)>=V_sig2_plot).nonzero()[0][0]],V_sig/V_sig.max())

plt.plot(t_arr2-t_arr2[(V_sig2>=np.nanmax(V_sig2*0.5)).nonzero()[0][0]],V_sig2/np.nanmax(V_sig2))
#%%

plt.figure()
dt = t_arr2[1]-t_arr2[0]
V_broadened = scipy.ndimage.gaussian_filter1d(V_sig2/np.nanmax(V_sig2),sigma=0.25/2.355/dt)
t_shift = t_arr2[(V_broadened>=np.nanmax(V_broadened*0.5)).nonzero()[0][0]]
plt.plot(t_arr2-t_shift,V_broadened/np.nanmax(V_broadened))