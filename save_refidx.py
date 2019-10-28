#!/Users/wdai11/anaconda3/bin/python3

import meep as mp
import numpy as np
import math
import cmath
import sys

sys.path.insert(0,'/Users/wdai11/function')
import my_output as my
import materials_20190221 as mym
sys.path.insert(0,'/Users/wdai11/meep-master/meep-master/python/examples')
import materials_library as mat


def RefIdx_from_medium(medium,name):
    wls=np.arange(0.5,10,0.02)
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


RefIdx_from_medium(mat.Al,'Al')
RefIdx_from_medium(mat.Pd,'Pd')
RefIdx_from_medium(mat.Au,'Au')
RefIdx_from_medium(mat.Ag,'Ag')

# Ref_from_medium(mym.In295,'In295')
# Ref_from_medium(mym.In4p2,'In4p2')
# Ref_from_medium(mat.Al,'Al', n0=1.5)
# Ref_from_medium(mat.Au,'Al', n0=1.5)
