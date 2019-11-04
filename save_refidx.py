#!/usr/bin/env python
"""Calcluate and save the permittivity and refractive index directly
from Meep mdium class.

"""
import meep as mp
import numpy as np
import math
import cmath
import sys

sys.path.append("/Users/wdai11/python-study")
import  MeepMedium.dwt_materials as mym
import  MeepMedium.dwt_output as my
import meep.materials as mat


def RefIdx_from_medium(medium,name):
    """Calculate the refractive index and reflection from the medium;
    save the result to a file"""
    wls=np.arange(0.5,10,0.5)
    nref=[mym.get_refractive_index(1/wl,medium) for wl in wls]
    n0=3.8
    ref=[abs((n-n0)/(n+n0))**2 for n in nref]

    nref=np.array(nref)
    ref=np.array(ref)
    fname='Ref_{}.dat'.format(name)
    data=np.column_stack((wls,nref.real,nref.imag,ref))
    print('{} Refectlion'.format(name))
    my.matrix_output(fname,data,"{:10.5e}",sig=None)
    print('-'*30)

#########################
def main(args):
    RefIdx_from_medium(mat.Al,'Al')
    RefIdx_from_medium(mym.In4p2,'In')
    RefIdx_from_medium(mym.lossy_Pd(0.5),'Pd')

#########################
# main function
if __name__=='__main__':
    main(sys.argv)
