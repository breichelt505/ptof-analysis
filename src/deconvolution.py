import numpy as np



def time_basis_conv_fullmode(t1, t2):
    if len(t1) != len(t2):
        raise ValueError("len(t1) must equal len(t2)")
    n_full = len(t1) + len(t2) - 1
    tf_full = np.linspace(t1[0]+t2[0], t1[-1]+t2[-1], n_full)
    return tf_full

def time_basis_conv_samemode(t1, t2):
    if len(t1) != len(t2):
        raise ValueError("len(t1) must equal len(t2)")
    tf_full = time_basis_conv_fullmode(t1, t2)
    n_full = len(t1) + len(t2) - 1
    n_tf = max([len(t1), len(t2)])
    tf_same = tf_full[int(np.floor(n_full/4)):int(np.floor(n_full/4))+n_tf]
    return tf_same


def time_basis_deconv_samemode(tf, t1):
    if len(t1) != len(tf):
        raise ValueError("len(t1) must equal len(t2)")
    n2 = len(tf)
    dt = -(t1[0]+t1[-1])/2 
    if len(tf) % 2 == 0:
        dt += (t1[1]-t1[0])/2
    t2_same = np.linspace(tf[0]+dt, tf[-1]+dt, n2)    
    return t2_same


#%%

# %%
