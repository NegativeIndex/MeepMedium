"""A meep simulation file. We can calculate the reflection from a
medium based on the time-domain simulaiton. Then compare the result
with the analytic results from *save_refidx*.

"""

import meep as mp
import argparse
import math
import cmath
import sys,os
import numpy as np
import datetime

sys.path.append("/Users/wdai11/python-study")
import  MeepMedium.dwt_materials as mym
sys.path.insert(0,'/Users/wdai11/function')
import my_output as my
import meep.materials as mat

#########################################
# common parameters
#########################################
class common:
    pols=["Ez"]
    geoms=["Empty","Pda"]
    fcens=np.arange(0.25,0.28,0.02)
    dfs=np.ones(fcens.size)*0.02
  
    pol='Ez'
    geom='Pda'
    fcen=0.25
    df=0.02
    theta=0
    r=1

    if geom == 'Empty':
        sig='{}_{}_F{:0.3f}Th{:02.1f}'.format(geom,pol,fcen,theta)
    else :
        sig='{}{:0.2f}_{}_F{:0.3f}Th{:02.1f}'.format(geom,r,pol,fcen,theta)


#########################################
# simulation function
#########################################
def simulation_fun():
    resolution=200
    geom=common.geom
    sig=common.sig

    # source parameters
    fcen=common.fcen
    df=common.df
    df2=df*1.4             # source frequency width
    nfreq=101              # number of frequency bins
 
    pol=common.pol         # polarization
    src_cmpt=eval('mp.{}'.format(pol))   

    theta=common.theta     # incident angle 
    theta_r=math.radians(theta)

    # material
    n_GaSb=3.8   
    mat_src=mp.Medium(index=n_GaSb)
    if geom=="Pda":
        mat_metal=mym.lossy_Pd(common.r)
    elif geom in ("Pd","Empty"):
        mat_metal=mat.Pd
    else:
        sys.exit("Wrong materials")
   
    # geometry
    # from bottom to top
    sx=0.1
    sxx=0.1 # periodic along x direction

    d_back=0.3
    d_metal=0.7
    d_src=1
    dpml=2
    sy=d_back+d_metal+d_src
    syy=sy+2*dpml

    cell_size = mp.Vector3(sxx,syy,0)
    boundary_layers = [mp.PML(thickness=dpml,direction=mp.Y)]

    if geom=="Empty":
        geometry=[mp.Block(size=mp.Vector3(sxx,syy,mp.inf),
                           center=mp.Vector3(0,0,0),
                           material=mat_src)]
    else:        
        tts=[d_back,d_metal]
        ccs=my.center_from_thickness(tts,-sy/2)
    
        geometry= [mp.Block(size=mp.Vector3(sxx,syy,mp.inf),
                             center=mp.Vector3(0,0,0),
                             material=mat_src),
                   mp.Block(size=mp.Vector3(sxx,tts[1],mp.inf),
                            center=mp.Vector3(0,ccs[1],0),
                            material=mat_metal),]
   

    # bloch k vector
    k=mp.Vector3(math.sin(theta_r),math.cos(theta_r),0).scale(
        fcen*mym.get_refractive_index(fcen,mat_src).real)
    amp_func=lambda x: cmath.exp(1j*2*math.pi*k.dot(x))

    # oblique source
    y0=-sy/2+d_back+d_metal # metal/GaSb interface
    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df2),
                         component=src_cmpt,
                         center=mp.Vector3(0,y0+d_src*0.8),
                         size=mp.Vector3(sxx,0),
                         amp_func=amp_func
                     )]

    # setup simulations
    sim = mp.Simulation(cell_size=cell_size,
                        geometry=geometry,
                        boundary_layers=boundary_layers,
                        sources=sources,
                        k_point=k,
                        force_complex_fields=True,
                        filename_prefix=sig,
                        resolution=resolution)

    # flux 
    y_frs=[y0+d_src*0.4,]
    frs=[mp.FluxRegion(center=mp.Vector3(0,y,0),
                       size=mp.Vector3(sxx,0,0),
                       weight=-1) for y in y_frs]
    trans=[sim.add_flux(fcen,df,nfreq,fr) for fr in frs]


    ####
    # define step function to display fluxes
    my_display_count=0
    step_time_flush=10
    step_time_print=200
    step_time_terminate=50
    def my_display_fluxes(sim):
        nonlocal my_display_count
        print('='*40)
        print('='*40)
  
        freqs = mp.get_flux_freqs(trans[0])
        fluxes=[mp.get_fluxes(tran) for tran in trans]
        data=np.column_stack((freqs,*fluxes,))
        my.matrix_output(None,data,"{:10.3e}","flux")
      
        my_display_count+=1
        print('='*40)
        print('No. {} display at t={}'.format(
            my_display_count,step_time_print))
        my.my_flush_step(sim)

    # run simulations
    pt_field=mp.Vector3(0.1*sxx,y0+d_src*0.4,0) # monitor point
    sim.run(# mp.at_beginning(mp.output_epsilon),
        mp.at_every(step_time_flush, my.my_flush_step),
        mp.at_every(step_time_print, my_display_fluxes),
        until_after_sources=mp.stop_when_fields_decayed(
            step_time_terminate, src_cmpt, pt_field, 1e-9))
    sys.stdout.flush()

    # collect data
    freqs = mp.get_flux_freqs(trans[0])
    fluxes=[mp.get_fluxes(tran) for tran in trans]
    data=np.column_stack((freqs,*fluxes,))


    # output
    fname=sig+".dat"
    my.matrix_output(fname,data,"{:10.3e}","flux")

######################### 
# main function
##########################
if __name__=='__main__':
    time0=datetime.datetime.now()
    print('-'*50)
    print(common.sig)
    simulation_fun()

    time1=datetime.datetime.now()
    dtime=time1-time0
    print("Simulation used {:0.2f} seconds".format(
        dtime.total_seconds()))
    sys.stdout.flush()

