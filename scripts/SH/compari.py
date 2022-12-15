
from dmipy.core import modeling_framework
from dmipy.core.acquisition_scheme import acquisition_scheme_from_bvalues
from os.path import join
import numpy as np
from dipy.io.image import load_nifti
from dmipy.signal_models import gaussian_models
from dmipy.core.modeling_framework import MultiCompartmentModel
from dmipy.core import modeling_framework
import nibabel as nib
import os
import sys 

try:
   dmri= sys.argv[1]
   bval= sys.argv[2]
   bvec= sys.argv[3]
   output= sys.argv[4]
except ValueError:
   sys.exit(1)
   

class DwiInfo:
    def __init__(self,dmri,bval,bvec,output,delta = 0.01293,Delta = 0.03666):
        self.dmri= dmri
        self.bval= bval
        self.bvec = bvec
        self.delta = delta
        self.Delta = Delta
        self.output = output
    
    def getScheme(self):
        acquisition_path = modeling_framework.GRADIENT_TABLES_PATH
        bvalues = np.loadtxt(self.bval)
        bvalues_SI = bvalues * 1e6  # SI units as s/m^2
        gradient_directions = np.loadtxt(self.bvec)
        gradient_directions[2] *= -1
        return (acquisition_scheme_from_bvalues(bvalues_SI, gradient_directions.T, self.delta, self.Delta))
        
    def loadNifty(self):
        return (load_nifti(self.dmri))
        
def isotropicFit(dmri,bval,bvec,output):
    info = DwiInfo(dmri,bval,bvec,output)
    scheme = info.getScheme()
    scheme.print_acquisition_info
    dwi, affine = info.loadNifty()

    ###########################################
    if (os.path.exists(output) == False):
       os.makedirs(output) 
       print("The", output, " is created!")

    ###########################################
    ball = gaussian_models.G1Ball()

    my_mod = modeling_framework.MultiCompartmentModel([ball])
    my_mod.set_parameter_optimization_bounds('G1Ball_1_lambda_iso', [0.1e-9, 3e-9])


    mymod_fit = my_mod.fit(scheme,dwi, mask=dwi[..., 0]>0)
    iso = mymod_fit.fitted_parameters['G1Ball_1_lambda_iso']

    img = nib.Nifti1Image(iso, affine)
    nib.save(img,output+'/iso.nii.gz')
    
    print ('the isotopic fitting is done!')
    
isotropicFit(dmri,bval,bvec,output)

