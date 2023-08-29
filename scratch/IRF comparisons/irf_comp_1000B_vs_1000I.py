# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 12:28:37 2023

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
        
def find_waveforms(expname):
    shotnum = expname.split(" ")[0]
    
    fname =  "C:\\Users\\blr\\Dropbox (MIT)\\PTOF_GUI_Ben\\Analysis Files\\"
    fname += expname
    fname += f"\\TD_TC090-078_PCD_SCOPE-01-DB_{shotnum}.h5"
    try:
        with h5py.File(fname) as f:
            t_tek, V_tek = get_wf_data(f)
    except:
        t_tek, V_tek = [np.nan, np.nan]
    
    fname =  "C:\\Users\\blr\\Dropbox (MIT)\\PTOF_GUI_Ben\\Analysis Files\\"
    fname += expname
    fname += f"\\TD_TC090-078_PCD_CRT-SCOPE-01-DB_{shotnum}.h5"
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

for d in pathlib.Path(r"C:\Users\blr\Dropbox (MIT)\PTOF_GUI_Ben\Analysis Files").iterdir():
    expname = d.stem
    if "_ARC_" not in expname and "AuDisk" not in expname:
        continue
    for f in d.iterdir():
        if "Traveler" in f.stem:
            data = pd.read_excel(f)
            detector = data["Unnamed: 4"][12]
            print(detector, expname)
    if detector.upper() not in ["1000I"]:
        continue
    t_tek, V_tek, t_FTD, V_FTD = find_waveforms(expname)
    try:
        t_tek= shift_t_rising(t_tek, V_tek,0.5)
        ax.plot(t_tek, V_tek/V_tek.max(), "r")
        plt.xlim(-2e-9,10e-9)
        plt.pause(0.2)
        keep_line = input("Do you want this line (y/n) : ")
        if keep_line == "n":
            ax.lines[-1].remove()
        # plt.xlabel("time (ns)")
        # plt.ylabel("normalized signal")
    except:
        pass
    
    try:
        t_FTD= shift_t_rising(t_FTD, V_FTD/V_FTD.max(),0.5)
        ax.plot(t_FTD, V_FTD,"b")
        plt.xlim(-2e-9,10e-9)
        plt.pause(0.2)
        keep_line = input("Do you want this line (y/n) : ")
        if keep_line == "n":
            ax.lines[-1].remove()
    except:
        pass
    # plt.draw()
    


expname = "N220406-002-999 Fa_LSR_ARC_Point_S26a"
t_tek, V_tek, t_FTD, V_FTD = find_waveforms(expname)
l, = ax.plot(t_tek, V_tek, lw=2)


class ShiftScale:
    t = t_tek
    y = V_tek
    t_shift = 0.0
    y_scale = 1.0

    def shift(self, shiftval):
        self.t_shift = float(shiftval)
        l.set_xdata(self.t+self.t_shift)
        plt.draw()

    def scale(self, scaleval):
        self.y_scale = float(scaleval)
        l.set_ydata(self.y*self.y_scale)
        plt.draw()

callback = ShiftScale()
axshift = fig.add_axes([0.3, 0.03, 0.2, 0.06])
axscale = fig.add_axes([0.7, 0.03, 0.2, 0.06])
tscale = TextBox(axscale,"scale",initial="1.0",label_pad=0.05)
tscale.on_submit(callback.scale)
tshift = TextBox(axshift,"shift", initial="0.0",label_pad=0.05)
tshift.on_submit(callback.shift)

plt.show()