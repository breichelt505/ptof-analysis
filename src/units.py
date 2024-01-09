#%%
import astropy.units as u
import astropy.constants as c


tstar = 1*u.ns
vstar = 1*c.c
mstar = 1*c.m_p
Vstar = 1*u.V

lstar = (tstar*vstar).si
Estar = (mstar*vstar**2).si
qstar = Estar / lstar
Istar = qstar / tstar


# %%

def unitful_to_analysis(qty):
    if type(qty) not in [u.quantity.Quantity, c.codata2018.CODATA2018]:
        return qty
    else:
        powers = qty.si.unit.decompose().powers
        bases = qty.si.unit.decompose().bases

        scale = 1.0
        for i, base in enumerate(bases):
            if base == u.m:
                scale *= (u.m/lstar).si.value**powers[i]
            if base == u.A:
                scale *= (u.A/Istar).si.value**powers[i]
            if base == u.kg:
                scale *= (u.kg/mstar).si.value**powers[i]
            if base == u.s:
                scale *= (u.s/tstar).si.value**powers[i]
        
        return qty.si.value * scale



# %%
