import numpy as np
import scipy.optimize as opt
import scipy.stats
import scipy.interpolate as interp
import numba

def subtract_const_bg(t, V, t_fit_min=-np.inf, t_fit_max=-1000):
    mask = (t>t_fit_min)*(t<t_fit_max)
    return V - np.mean(V[mask])

def fit_gaussian_fidu(t, V, t_fit_min=-np.inf, t_fit_max=-300, amp_thresh = 0.2):
    mask = (t>t_fit_min)*(t<t_fit_max)
    t_max_0 = (t[mask])[np.argmax(V[mask])]
    V_max_0 = V[mask].max()
    sigma_0 = 0.1

    i_lb = (V > V_max_0*amp_thresh).nonzero()[0][0]
    i_ub = (V > V_max_0*amp_thresh).nonzero()[0][-1]

    t_fit = t[i_lb:i_ub]
    V_fit = V[i_lb:i_ub]

    p0 = [V_max_0, sigma_0, t_max_0]

    def gaussian(t, A, sigma, mu):
        return A*np.exp(-(t-mu)**2/2/sigma**2)
    
    popt, pcov = opt.curve_fit(gaussian, t_fit, V_fit, p0)

    return popt, pcov, gaussian

def fit_exponnorm_fidu(t, V, t_fit_min=-np.inf, t_fit_max=-300, amp_thresh = 0.2):
    mask = (t>t_fit_min)*(t<t_fit_max)
    t_max_0 = (t[mask])[np.argmax(V[mask])]
    V_max_0 = V[mask].max()
    sigma_0 = 0.1

    i_lb = (V > V_max_0*amp_thresh).nonzero()[0][0]
    i_ub = (V > V_max_0*amp_thresh).nonzero()[0][-1]

    t_fit = t[i_lb:i_ub]
    V_fit = V[i_lb:i_ub]

    p0 = [V_max_0, sigma_0, t_max_0, sigma_0]

    def exponnorm(t, A, sigma, mu, lamb):
        fit = scipy.stats.exponnorm.pdf(t, K=1/sigma/lamb, loc=mu, scale=sigma)
        return A*fit/fit.max()
    
    popt, pcov = opt.curve_fit(exponnorm, t_fit, V_fit, p0)

    return popt, pcov, exponnorm


def shift_and_scale(t_ref, V_ref, V_ref_clip, t_scale, V_scale, t_fit_min, t_fit_max):
    mask = (t_ref>t_fit_min)*(t_ref<t_fit_max)

    def shift_scale_fn(t, shift, scale):
        return np.interp(t+shift, t_scale, V_scale)*scale
    
    sigma = np.ones_like(t_ref)
    sigma[V_ref>V_ref_clip*0.95] = 1000

    popt, pcov = opt.curve_fit(shift_scale_fn, t_ref[mask], V_ref[mask], [0.0,1.0], sigma=sigma[mask])

    return popt, pcov, shift_scale_fn

@numba.njit()
def switching_function(V, V_clips, start_thresh = 0.4, end_thresh = 0.8):
    weights = np.zeros_like(V_clips)
    clips_from_0 = np.zeros(len(V_clips)+1)
    clips_from_0[1:] = V_clips

    V_starts = np.diff(clips_from_0)*start_thresh + clips_from_0[:-1]
    V_ends = np.diff(clips_from_0)*end_thresh + clips_from_0[:-1]


    passed_switch_starts = (V>=V_starts).nonzero()[0]
    passed_switch_ends = (V>=V_ends).nonzero()[0]
    
    if len(passed_switch_starts) == 0 and len(passed_switch_ends)==0:
        weights[0] = 1.0
    elif len(passed_switch_ends) == 0:
        dV = V_ends[0] - V_starts[0]
        weights[1] = (V - V_starts[0]) / dV
        weights[0] = 1 - weights[1]
    elif len(passed_switch_starts) == len(V_clips) or len(passed_switch_ends)==len(V_clips):
        i = passed_switch_starts[-1]
        weights[-1] = 1.0
    elif len(passed_switch_starts) == len(passed_switch_ends):
        i = passed_switch_starts[-1]
        weights[i+1] = 1.0
    elif len(passed_switch_starts) != len(passed_switch_ends):
        i = passed_switch_starts[-1]
        dV = V_ends[i] - V_starts[i]
        weights[i+1] = (V - V_starts[i]) / dV
        weights[i] = 1 - weights[i+1]

    return weights

@numba.njit()
def stitch_waveforms(Vs, V_clips):
    V_stitch = np.zeros_like(Vs[0])
    for i in range(len(Vs[0])): 
        Vis = np.array([V[i] for V in Vs])
        V_max = max(Vis)
        V_stitch[i] = np.sum(switching_function(V_max, V_clips)*Vis)
    return V_stitch

