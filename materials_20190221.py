import meep as mp
import numpy as np
import math
import cmath
import sys

######################################
# function to calculate permittivity
######################################

def symmetric_matrix_from_2vectors(v1=[1,1,1],v2=[0,0,0]):
    a,b,c=v1
    u,v,w=v2
    return np.complex128(np.matrix([[a,u,v],[u,b,w],[v,w,c]]))

# reduce matrix to a scalar if it is a constant matrix
def reduce_matrix_to_scalar(M):
    a=M.diagonal().mean()
    if (M.shape[0]==M.shape[1]) and np.allclose(M, a*np.eye(M.shape[0])):
        return a
    else:
        return None

def get_Permittivty_from_Medium(freq,medium):
    sigma=symmetric_matrix_from_2vectors(medium.D_conductivity_diag)
    prod0=np.identity(3)+1j*sigma/(2*math.pi*freq)
    
    epsi0=symmetric_matrix_from_2vectors(medium.epsilon_diag,
                                         medium.epsilon_offdiag)
    for susc in medium.E_susceptibilities:
        f0=susc.frequency
        gamma=susc.gamma
        sigma=symmetric_matrix_from_2vectors(susc.sigma_diag,
                                             susc.sigma_offdiag)
        if susc.__class__.__name__=='LorentzianSusceptibility':
            epsi0+=sigma*f0*f0/(f0*f0-freq*freq-1j*freq*gamma)
        elif susc.__class__.__name__=='DrudeSusceptibility':
            epsi0+=sigma*f0*f0/(-freq*freq-1j*freq*gamma)
        else:
            return None

    return prod0.dot(epsi0)


def get_refractive_index(freq,medium):
	epsi_t=get_Permittivty_from_Medium(freq,medium)
	epsi=reduce_matrix_to_scalar(epsi_t)
	return cmath.sqrt(epsi)

########################################################
########################################################
########################################################
########################################################
# material library
########################################################
########################################################
########################################################
########################################################

um_scale = 1.0
eV_um_scale = um_scale/1.23984193
metal_range = mp.FreqRange(min=um_scale/12.4, max=um_scale/0.2)
############################
# Indium 295K
############################
metal_range = mp.FreqRange(min=um_scale/10, max=um_scale/0.1)
In295_plasma_frq = 7.462*um_scale
In295_frq0 = 1e-10
In295_gam0 = 0.147*um_scale
In295_sig0 = In295_plasma_frq**2/In295_frq0**2

In295_susc = [mp.DrudeSusceptibility(frequency=In295_frq0, 
                                  gamma=In295_gam0, sigma=In295_sig0),
       ]

In295 = mp.Medium(epsilon=3.363, E_susceptibilities=In295_susc, 
               valid_freq_range=metal_range)


############################
# Indium 4.2K
############################
metal_range = mp.FreqRange(min=um_scale/10, max=um_scale/0.1)
In4p2_plasma_frq = 6.8353*um_scale
In4p2_frq0 = 1e-10
In4p2_gam0 = 0.0471*um_scale
In4p2_sig0 = In4p2_plasma_frq**2/In4p2_frq0**2

In4p2_susc = [mp.DrudeSusceptibility(frequency=In4p2_frq0, 
                                  gamma=In4p2_gam0, sigma=In4p2_sig0),
       ]

In4p2 = mp.Medium(epsilon=1.000, E_susceptibilities=In4p2_susc, 
               valid_freq_range=metal_range)

############################
# palladium (Pd)
############################
Pd_plasma_frq = 9.72*eV_um_scale
Pd_f0 = 0.330
Pd_frq0 = 1e-10
Pd_gam0 = 0.008*eV_um_scale
Pd_sig0 = Pd_f0*Pd_plasma_frq**2/Pd_frq0**2
Pd_f1 = 0.649
Pd_frq1 = 0.336*eV_um_scale      # 3.690 um
Pd_gam1 = 2.950*eV_um_scale
Pd_sig1 = Pd_f1*Pd_plasma_frq**2/Pd_frq1**2
Pd_f2 = 0.121
Pd_frq2 = 0.501*eV_um_scale      # 2.475 um
Pd_gam2 = 0.555*eV_um_scale
Pd_sig2 = Pd_f2*Pd_plasma_frq**2/Pd_frq2**2
Pd_f3 = 0.638
Pd_frq3 = 1.659*eV_um_scale      # 0.747 um
Pd_gam3 = 4.621*eV_um_scale
Pd_sig3 = Pd_f3*Pd_plasma_frq**2/Pd_frq3**2
Pd_f4 = 0.453
Pd_frq4 = 5.715*eV_um_scale      # 0.217 um
Pd_gam4 = 3.236*eV_um_scale
Pd_sig4 = Pd_f4*Pd_plasma_frq**2/Pd_frq4**2

Pd_susc = [mp.DrudeSusceptibility(frequency=Pd_frq0, 
                                  gamma=Pd_gam0, sigma=Pd_sig0),]
Pd = mp.Medium(epsilon=1.0, E_susceptibilities=Pd_susc, 
               valid_freq_range=metal_range)
