# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 14:01:35 2023

@author: blr
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pathlib
import scipy.io
import copy
import matplotlib.pylab as pl
import scipy.ndimage
import skimage

def get_waveform(file):
    other_IRF = pd.read_excel(file, header=3)
    
    t = other_IRF["Time (s)"]*1e9
    V = other_IRF["Signal (V)"]
    
    return t, V
 
plt.figure()
waveforms = [24,25,26,28]
exp =["N230611-001-999","N230708-001-999","N230723-001-999","N231001-003-999"]
GRH_BW = np.r_[270,183,147,np.nan]
colors = ["mediumseagreen","blue","red","black"]
distances = np.r_[811,420,np.nan,np.nan]
PTOF_20to80 = np.r_[275,230,170,225]

#%%
for i,wf in enumerate(waveforms):
    try:
        other_IRF = pd.read_excel(f"C:\\Users\\blr\\Downloads\\waveform ({wf}).xlsx",
                                  header=3)
        
        # plt.figure()
        t = np.array(other_IRF["Time (s)"])*1e9
        V = np.array(other_IRF["Signal (V)"])
        
        plt.plot(t-t[(V>=V.max()*0.5).nonzero()[0][0]],V/V.max(),color=colors[i],lw=5)
    except:
        pass

