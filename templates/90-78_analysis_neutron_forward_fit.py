#%%
import src.file_utilities as futil
import src.waveform_utilities as wfutil
import matplotlib.pyplot as plt
import numpy as np
import skimage
import astropy.units as u
import src.units
import astropy.constants as c
import src.deconvolution
import src.propagation
import src.source
import scipy.optimize as opt
import pandas as pd
import src.file_download

#%%
%reload_ext autoreload
%autoreload 2
%matplotlib qt

#%% Important values
shotnum = "N240104-001-999"
expdir = r"C:\Users\blr\Dropbox (MIT)\MIT\Research\Projects\ptof_analysis\scratch\manual_analysis_tests\N240104-001-999"
t_fidu_cal = -382.563
cable_file_path = r".\tables\cable_response\N231017-001-999_deconvolved_cable_chain_IRF.txt"
neut_sensitivity_path = r"tables\detector_response\sensitivities\Neutron Sensitivity.xlsx"
t_delay_cal = 210.692
dist = src.units.unitful_to_analysis((417.2)*u.mm)
T_DD = 3.03

#%% Read in files and perform stitching
shifts = np.r_[0.0, 0.145, 0.145+0.162, 0.145+0.162+0.155]
scales = np.r_[2, 4, 8, 16]
file = f"C:\\Users\\blr\\Dropbox (MIT)\\MIT\\Research\\Projects\\ptof_analysis\\scratch\\manual_analysis_tests\\{shotnum}\\TD_TC090-078_PCD_SCOPE-02-DB_{shotnum}.h5"
ts, Vs, sat_flags, vclip_lbs, vclip_ubs = futil.read_tek_file(file)
Vs_corr = [wfutil.subtract_const_bg(ts[i],Vs[i],-np.inf,-1000) for i in range(len(ts))]
Vs_corr = [np.interp(ts[i]+shifts[i], ts[i], Vs_corr[i])*scales[i] for i in range(len(ts))]
V_stitch = src.waveform_utilities.stitch_waveforms_twosided(Vs_corr,vclip_lbs*scales,vclip_ubs*scales, 0.9, 0.95)

#%% Get experiment data
data = futil.get_setup_info(expdir, "090-078")
print(data)
print("\n\n\n")
print(f"{shotnum} {src.file_download.get_exp_ID(shotnum,expdir)}")

# %% Find scales and shifts of fidus
shifts = []
scales = []
shift_and_scale_fns = []
for i, wf in enumerate(Vs):
    popt, pcov, shift_and_scale_fn = wfutil.shift_and_scale(ts[0],Vs[0],np.inf,ts[i],Vs[i],-560,-460)
    shifts.append(popt[0])
    scales.append(popt[1])

for i in range(3):
    print(f"shift and scale ch {i+1} to ch {i} is {shifts[i]-shifts[i+1]}, {scales[i]/scales[i+1]}")
    
print("\n")
shifts = []
scales = []
for i, wf in enumerate(Vs):
    popt, pcov, gauss = wfutil.fit_gaussian_fidu(np.array(ts[i]),np.array(Vs[i]),-560,-460, 0.2)
    shifts.append(popt[2])
    scales.append(popt[0])

for i in range(3):
    print(f"gauss fit ch {i+1} to ch {i} is {shifts[i]-shifts[i+1]}, {scales[i+1]/scales[i]}")


#%% Generate analysis signal by correcting fidu timing
plt.rcParams.update({"font.size" : 16})
t_shift = shifts[0]-t_fidu_cal
V_analysis = np.copy(V_stitch)
t_analysis = ts[0] - t_shift

plt.figure(figsize=(10,5))
plt.plot(t_analysis,V_analysis,"k-",lw=4)
for V in Vs_corr:
    plt.plot(t_analysis,V,"r",lw=2,alpha=0.3)
plt.xlim(200,250)
plt.title(f"{shotnum} Analysis Signal")
plt.xlabel("Time (ns)")
plt.ylabel("Signal (V)")
plt.tight_layout()
plt.savefig(expdir + r"/analysis_signal.pdf")
plt.savefig(expdir + r"/analysis_signal.png")

# %% Find cable chain IRF shape
cable_chain_IRF_matrix = np.loadtxt(cable_file_path)
t_cable_chain_IRF = cable_chain_IRF_matrix[:,0]
V_cable_chain_IRF = cable_chain_IRF_matrix[:,1]

# add up cable chain component delays
t_delay = 0.0
t_delay += 14.838 # snout cable -- 14017032	
t_delay += 15.49575 # PCD cable -- 13254732-200
t_delay += 30.624 # cable track cable -- rough guess
t_delay += 137.350 # infrastructure cable -- rough guess
t_delay += 0.0 # bias T -- ask eddie exact SN/TDR
t_delay += 12.21 # rack cable -- ask eddie exact SN
t_delay += 0.0 + 0.354 # attenuation -- mystery 3db atten + 20 dB barth 13390504
t_delay += 0.0 # clipper -- mystery, ask eddie exact SN
t_delay += 0.171 # splitter -- ask eddie exact SN
t_delay += 0.0 # cable -- ask eddie exact SN
t_delay += 0.0 # misc

t_cable_chain_IRF += t_delay - t_delay_cal
plt.figure()
plt.plot(t_cable_chain_IRF, V_cable_chain_IRF)


# %% Do DDn deconv
file_path = r".\tables\detector_response\1000I_N211201_deconvolved_IRF.txt"
V_detector_IRF = np.loadtxt(file_path)[:,1]
t_detector_IRF = np.loadtxt(file_path)[:,0]
trise = wfutil.find_rise_time(t_detector_IRF, V_detector_IRF, 0.1)
V_detector_IRF = np.interp(np.arange(-50,50,0.04)+trise, t_detector_IRF, V_detector_IRF,0,0)
t_detector_IRF = np.arange(-50,50,0.04)

# There is not a 3 cm W filter on this shot
# V_filter_IRF = np.zeros_like(t_detector_IRF)
# V_filter_IRF[1250] = 0.64/0.04
# t_fall = 0.071
# V_filter_IRF += 0.30*(t_detector_IRF>=0.0)*np.exp(-t_detector_IRF/t_fall)/t_fall
# t_fall = 0.300
# V_filter_IRF += 0.06*(t_detector_IRF>=0.0)*np.exp(-t_detector_IRF/t_fall)/t_fall

t_TOF_IRF = np.arange(-50,50,0.04)+10*dist
E, spec = src.source.ballabio_spectrum("DDn", T_DD, 1.0, 100)
E = src.units.unitful_to_analysis(E*u.MeV)
V_TOF_IRF = src.propagation.dYdt_IRF_vacprop_arbdYdE(t_TOF_IRF,dist,1.0,E,spec)
V_TOF_IRF[t_TOF_IRF<1.5*dist]=0.0

sens_data = pd.read_excel(neut_sensitivity_path)
E_sens = src.units.unitful_to_analysis(np.array(sens_data["E(MeV)"])*u.MeV)
sens = np.array(sens_data["Sens (b-MeV)"])
E_TOF = src.propagation.E_func(t_TOF_IRF,1.0,dist)
E_TOF[t_TOF_IRF<dist*1.5] = 0.0
sens_TOF = np.interp(E_TOF, E_sens, sens)
V_TOF_IRF *= sens_TOF


t_total_IRF = src.deconvolution.time_basis_conv_samemode(t_TOF_IRF,t_cable_chain_IRF)
V_total_IRF = np.convolve(V_TOF_IRF,V_cable_chain_IRF,"same")

t_total_IRF = src.deconvolution.time_basis_conv_samemode(t_total_IRF,t_detector_IRF)
V_total_IRF = np.convolve(V_total_IRF,V_detector_IRF,"same")

# t_total_IRF = src.deconvolution.time_basis_conv_samemode(t_total_IRF,t_detector_IRF)
# V_total_IRF = np.convolve(V_total_IRF,V_filter_IRF,"same")

# plt.figure()
# plt.plot(t_cable_chain_IRF, V_cable_chain_IRF)
# plt.plot(t_total_IRF, V_total_IRF)
# plt.plot(t_detector_IRF, V_filter_IRF)

# V_analysis[(t_analysis<232.0)] = 0.0

mask1 = (t_analysis>=233.0)*(t_analysis<=253.0)
mask2 = (t_total_IRF>=220.0)*(t_total_IRF<=240.0)
t_fit_mins = [229, 242]
t_fit_maxs = [234, 245]
V_bg_subtract = wfutil.subtract_exp_bg(t_analysis,V_analysis,t_fit_mins,t_fit_maxs)

def gauss(t, A, mu, sig):
    return A*np.exp(-(t-mu)**2/2/sig**2)/sig

def fit_fn(t, A, mu, sig):
    V_fit = np.convolve(V_total_IRF[mask2],gauss(t,A,mu,sig),"same")
    return V_fit

hist_t = src.deconvolution.time_basis_deconv_samemode(t_analysis[mask1],t_total_IRF[mask2])
# hist_V = 

popt, pcov = opt.curve_fit(fit_fn, hist_t, V_bg_subtract[mask1], p0=[1.0,5.2,0.1],sigma=0.06*np.ones(len(hist_t)),absolute_sigma=True)

plt.figure(figsize=(10,5))
plt.plot(t_analysis,V_bg_subtract,"r", lw=3)
V_fit = fit_fn(hist_t,*popt)
plt.plot(t_analysis[mask1], V_fit,"b", lw=3)
print(popt, np.sqrt(np.diag(pcov)))

plt.ylim(-0.1,0.6)
plt.xlim(234,244)
plt.title(f"{shotnum} DDn Signal and Fit")
plt.xlabel("Time (ns)")
plt.ylabel("Signal (V)")
plt.tight_layout()
plt.savefig(expdir + r"/DDn_fit.pdf")
plt.savefig(expdir + r"/DDn_fit.png")

# %%
laser_data = np.loadtxt(r"scratch\manual_analysis_tests\N240104-001-999\N240104-001-999_Total_Power.tsv",skiprows=1)
t_laser = laser_data[:,0]
P_laser = laser_data[:,2]

plt.figure(figsize=(10,5))
plt.plot(t_laser, P_laser,"k", lw=5, alpha=0.7)

# plt.ylim(-0.1,0.6)
plt.xlim(0, 8)
plt.title(f"{shotnum} Laser Power and Burn History")
plt.xlabel("Time (ns)")
plt.ylabel("Laser Power (TW)",color="k")
plt.twinx()
dYdt =  gauss(t_laser, *popt)
plt.ylabel("Burn History (normalized)",color="b")
plt.plot(t_laser,dYdt/dYdt.max(), "b--", lw=5, alpha=0.7)
plt.tick_params(axis="y", color="b", labelcolor="b")
plt.tight_layout()
plt.savefig(expdir + r"/laser_power.pdf")
plt.savefig(expdir + r"/laser_power.png")

# %%
