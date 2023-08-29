# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 14:17:40 2023

@author: blr
"""

from matplotlib.widgets import TextBox
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pathlib
import scipy.io
import copy
import matplotlib.pylab as pl
import scipy.ndimage
import h5py
import astropy.units as u
import astropy.constants as c
import propagation
import scipy.stats
import scipy.optimize as opt




def get_wf_data(f):
    i_interp = f["DATA"]["SHOT"]["CHANNEL_01"]["SCALED"]["INTERPOLATED_FLAGS"][:]
    i_sat = f["DATA"]["SHOT"]["CHANNEL_01"]["SCALED"]["SATURATED_FLAGS"][:]
    t = f["DATA"]["SHOT"]["CHANNEL_01"]["SCALED"]["X_AXIS"][:]
    V = f["DATA"]["SHOT"]["CHANNEL_01"]["SCALED"]["DATA"][:]
    V[i_interp==1] = V.max()
    V[i_sat==1] = V.max()
    return t,V
        
def find_waveforms(expname, DIM="90-078"):
    shotnum = expname.split(" ")[0]
    
    fname =  "C:\\Users\\blr\\Dropbox (MIT)\\PTOF_GUI_Ben\\Analysis Files\\"
    fname += expname
    fname += f"\\TD_TC0{DIM}_PCD_SCOPE-01-DB_{shotnum}.h5"
    try:
        with h5py.File(fname) as f:
            t_tek, V_tek = get_wf_data(f)
    except:
        t_tek, V_tek = [np.nan, np.nan]
    
    fname =  "C:\\Users\\blr\\Dropbox (MIT)\\PTOF_GUI_Ben\\Analysis Files\\"
    fname += expname
    fname += f"\\TD_TC0{DIM}_PCD_CRT-SCOPE-01-DB_{shotnum}.h5"
    try:
        with h5py.File(fname) as f:
            t_FTD, V_FTD = get_wf_data(f)
    except:
        t_FTD, V_FTD = [np.nan, np.nan]
        
    return t_tek, V_tek, t_FTD, V_FTD

def shift_t_rising(t,V,thresh=0.1):
    i_shift = (V>V.max()*thresh).nonzero()[0][0]
    dt = (t[i_shift]-t[i_shift-1])/(V[i_shift]-V[i_shift-1]) * (V[i_shift]-thresh*V.max())
    t -= t[i_shift] - dt
    return t  

def shift_t(t,V,thresh=0.25,remove_neg=True):
    t_pk = t[np.argmax(V)]
    if remove_neg:
        i_count = np.arange(max(0,np.argmax(V)-400,((V<0)*(t<t_pk)).nonzero()[0][-1]),min(np.argmax(V)+300,len(V),((V<0)*(t>t_pk)).nonzero()[0][0]))
    else:
        i_count = np.arange(max(0,np.argmax(V)-300),min(np.argmax(V)+300,len(V)))
    i_shift = (V[i_count]>V.max()*thresh).nonzero()[0][-1]
    t -= t[i_count][i_shift]
    return t[i_count], V[i_count] 


fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.2)

colors = {"1000I" : "r",
          "SC1000A" : "k",
          "1000G" : "mediumseagreen",
          "1000H" : "c",
          "1000K" : "m"}

#%%

for d in pathlib.Path(r"C:\Users\blr\Dropbox (MIT)\PTOF_GUI_Ben\Analysis Files").iterdir():
    expname = d.stem
    if "_ARC_" not in expname and "AuDisk" not in expname:
        continue
    for f in d.iterdir():
        if "Traveler" in f.stem:
            data = pd.read_excel(f)
            detector = data["Unnamed: 4"][12]
            print(detector, expname)
    if detector.upper() not in ["1000I","SC1000A","1000G","1000H","1000K"]:
        continue
    t_tek, V_tek, t_FTD, V_FTD = find_waveforms(expname)
    plot_FTD = "n"
    plot_tek = "n"
    
    try:
        t_tek= shift_t_rising(t_tek, V_tek,0.5)
        ax.plot(t_tek, V_tek/V_tek.max(), colors[detector.upper()],lw=4)
        plt.xlim(-2e-9,10e-9)
        plt.pause(0.2)
        plot_tek = input("Do you want this tek line (y/n) : ")
        ax.lines[-1].remove()
    except:
        pass
    
    try:
        t_FTD= shift_t_rising(t_FTD, V_FTD,0.5)
        ax.plot(t_FTD, V_FTD/V_FTD.max(),colors[detector.upper()],lw=4)
        plt.xlim(-2e-9,10e-9)
        plt.pause(0.2)
        plot_FTD = input("Do you want this line (y/n) : ")
        ax.lines[-1].remove()
            
    except:
        pass
    
    if plot_tek == "y" and plot_FTD == "y":
        ax.plot(t_tek, V_tek/V_tek.max(), colors[detector.upper()],lw=4)
        ax.plot(t_FTD[(V_FTD>=V_FTD.max()*0.8).nonzero()[0][-1]:]+0.3e-9, V_FTD[(V_FTD>=V_FTD.max()*0.8).nonzero()[0][-1]:]/V_tek.max()/10**1.5,colors[detector.upper()],lw=4)
    elif plot_tek == "y" and plot_FTD == "n":
        ax.plot(t_tek, V_tek/V_tek.max(), colors[detector.upper()],lw=4)
    elif plot_tek == "n" and plot_FTD == "y":
        ax.plot(t_FTD, V_FTD/V_FTD.max(), colors[detector.upper()],lw=4)
#%%
for d in pathlib.Path(r"C:\Users\blr\Dropbox (MIT)\PTOF_GUI_Ben\Analysis Files").iterdir():
    expname = d.stem
    if "_ARC_" not in expname and "AuDisk" not in expname:
        continue
    for f in d.iterdir():
        if "Traveler" in f.stem:
            data = pd.read_excel(f)
            detector = data["Unnamed: 4"][12]
            print(detector, expname)
    if detector.upper() not in ["1000I","SC1000A","1000G","1000H","1000K"]:
        continue
    t_tek, V_tek, t_FTD, V_FTD = find_waveforms(expname,"90-124")
    plot_FTD = "n"
    plot_tek = "n"
    
    try:
        t_tek= shift_t_rising(t_tek, V_tek,0.5)
        ax.plot(t_tek, V_tek/V_tek.max(), colors[detector.upper()],lw=4,linestyle="--")
        plt.xlim(-2e-9,10e-9)
        plt.pause(0.2)
        plot_tek = input("Do you want this tek line (y/n) : ")
        ax.lines[-1].remove()
    except:
        pass
    
    try:
        t_FTD= shift_t_rising(t_FTD, V_FTD,0.5)
        ax.plot(t_FTD, V_FTD/V_FTD.max(),colors[detector.upper()],lw=4,linestyle="--")
        plt.xlim(-2e-9,10e-9)
        plt.pause(0.2)
        plot_FTD = input("Do you want this line (y/n) : ")
        ax.lines[-1].remove()
            
    except:
        pass
    
    if plot_tek == "y" and plot_FTD == "y":
        ax.plot(t_tek, V_tek/V_tek.max(), colors[detector.upper()],lw=4)
        ax.plot(t_FTD[(V_FTD>=V_FTD.max()*0.8).nonzero()[0][-1]:]+0.3e-9, V_FTD[(V_FTD>=V_FTD.max()*0.8).nonzero()[0][-1]:]/V_tek.max()/10**1.5,colors[detector.upper()],lw=4)
    elif plot_tek == "y" and plot_FTD == "n":
        ax.plot(t_tek, V_tek/V_tek.max(), colors[detector.upper()],lw=4)
    elif plot_tek == "n" and plot_FTD == "y":
        ax.plot(t_FTD, V_FTD/V_FTD.max(), colors[detector.upper()],lw=4)
    
    
    # plt.draw()
    
