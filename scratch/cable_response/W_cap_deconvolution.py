# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 13:26:57 2023

@author: blr
"""

import numpy as np
import matplotlib.pyplot as plt
import src.cable
import astropy.units as u
import astropy.constants as c
import os
import pathlib
import skimage
import pandas as pd
import scipy.interpolate as interp
#%%

file_path = r"C:\Users\blr\Dropbox (MIT)\MIT\Research\Projects\ptof_analysis\tables\detector_response\megaray_05-22\SC1000A_avg_xrays_nofilter.txt"
psf = np.loadtxt(file_path)[15990:16100,1]
t = np.loadtxt(file_path)[15990:16100,0]*1e9

file_path = r"C:\Users\blr\Dropbox (MIT)\MIT\Research\Projects\ptof_analysis\tables\detector_response\megaray_05-22\SC1000A_avg_electrons_on_cap.txt"
full = np.loadtxt(file_path)[15990:16100,1]

deconv = skimage.restoration.richardson_lucy(full, psf, clip=False,num_iter=18)

plt.figure()
plt.plot(t,psf,"r:")
plt.plot(t,full,"k--")
plt.plot(t,deconv,"b")
reconstructed = np.convolve(psf,deconv,mode="same")
plt.plot(t,reconstructed/reconstructed.max(),"c--")

plt.figure()
file_path = r"C:\Users\blr\Dropbox (MIT)\MIT\Research\Projects\ptof_analysis\tables\detector_response\megaray_05-22\1000F_avg_xrays_nofilter.txt"
V_1000F = np.loadtxt(file_path)[15850:16900,1]
t_1000F = np.loadtxt(file_path)[15850:16900,0]*1e9
reconstructed = np.convolve(deconv,V_1000F,mode="full")
plt.plot(deconv/deconv.max())
plt.plot(V_1000F/V_1000F.max())

plt.figure()
dt = t[1]-t[0]
t_recon = np.linspace(0,dt*(len(reconstructed)-1),len(reconstructed))
plt.plot(t_recon-0.541,reconstructed/reconstructed.max())
plt.plot(t_1000F,V_1000F)

#%% Deconvolve for david
file_path = r"C:\Users\blr\Dropbox (MIT)\MIT\Research\Projects\ptof_analysis\tables\detector_response\megaray_05-22\SC1000A_avg_electrons_on_cap.txt"
megaray_IRF = np.loadtxt(file_path)
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


plt.figure()
V_sig = np.fft.irfft(H_tot*np.fft.rfft(V))
# V_sig[(t_arr>35)]=0.0
# V_sig[(t_arr<16)]=0.0
# plt.plot(t_arr,V_sig)

other_IRF = pd.read_excel("C:\\Users\\blr\\Downloads\\waveform (19).xlsx",
                          header=3)

t_david_data = other_IRF["Time (s)"]*1e9
V_david_data = other_IRF["Signal (V)"]-0.0005
# V_david_data = V_david_data[(t_david_data<=89)*(t_david_data>=79.0)]
V_david_data[np.logical_or(t_david_data>89,t_david_data<85.0)] = 0.0
# t_david_data = t_david_data[(t_david_data<=89)*(t_david_data>=79.0)]
plt.plot(t_david_data,V_david_data)

V_sig -= 0.004 - 0.0002
V_IRF_interp = interp.interp1d(t_arr, V_sig,bounds_error=False,fill_value=0.0)(t_david_data)
V_IRF_interp[t_david_data>26.3]=0.0
plt.plot(t_david_data,V_IRF_interp)

V_IRF_interp = V_IRF_interp[(t_david_data<=30)*(t_david_data>=10)]
V_david_data = V_david_data[(t_david_data<=95)*(t_david_data>=75)]
t_david_data = t_david_data[(t_david_data<=30)*(t_david_data>=10)]

#%%


deconv = skimage.restoration.richardson_lucy(V_david_data, V_IRF_interp, clip=False,num_iter=500)
plt.scatter(t_david_data,V_david_data/V_david_data.max())
# plt.plot(t_david_data,V_IRF_interp/V_IRF_interp.max())
# plt.plot(t_david_data,deconv/deconv.max(),"b",lw=5)
# recon = np.convolve(deconv,V_IRF_interp,mode="same")
plt.plot(t_david_data,recon/recon.max(),"r--")

plt.xlabel("time (ns) [arbitrary offset]")
plt.ylabel("signal (arb)")

