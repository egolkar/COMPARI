
from dmipy.core import modeling_framework
from dmipy.core.acquisition_scheme import acquisition_scheme_from_bvalues
from os.path import join
import numpy as np
from dipy.io.image import load_nifti
from dmipy.core.acquisition_scheme import gtab_dmipy2dipy
from dipy.reconst.csdeconv import auto_response_ssst
from dmipy.core import modeling_framework
import nibabel as nib
import argparse
import os
import sys 
import pickle

class DwiInfo:
    def __init__(self,args):
        self.dmri= args.dMRI_filename
        self.mask = args.mask
        self.bval= args.bvals
        self.bvec = args.bvecs
        self.delta = args.delta
        self.separation = args.separation
    
    def loadDWI(self):
        return (load_nifti(self.dmri))
    
    def getMask(self):
        return (load_nifti(self.mask))
        
    def getScheme(self):
            acquisition_path = modeling_framework.GRADIENT_TABLES_PATH
            bvalues = np.loadtxt(self.bval)
            bvalues_SI = bvalues * 1e6  # SI units as s/m^2
            gradient_directions = np.loadtxt(self.bvec)
            if (self.delta == None or self.separation == None):
                return acquisition_scheme_from_bvalues(bvalues_SI, gradient_directions.T)
            return (acquisition_scheme_from_bvalues(bvalues_SI, gradient_directions.T, self.delta, self.separation))

    def responseFunction(self):
            scheme = self.getScheme()
            dwi, affine = self.loadDWI()           
            gtab = gtab_dmipy2dipy(scheme)
            response, ratio = auto_response_ssst(gtab, dwi, roi_center=None, roi_radii=10, fa_thr=0.7)
            lambdas = response[0]
            return lambdas[0] * 1e-6

def argumentParser():
    p = argparse.ArgumentParser()
    p.add_argument('-i','--data',action='store',metavar='data',dest='dMRI_filename',
                    type=str, required=True, 
                    help='Input dMRI data file (Nifti format).'
                    )
                    
    p.add_argument('-m','--mask',action='store',metavar='mask', dest='mask',
                    type=str, required=True, 
                    help='Brain mask file (Nifti format).'
                    )
    p.add_argument('-b','--bvals', action='store', metavar='bvals', dest='bvals',
                    type=str, required=True, default=None,
                    help='B-values given in s/mm^2 (.bval file).'
                    )
    p.add_argument('-r','--bvecs', action='store', metavar='bvecs', dest='bvecs',
                    type=str, required=True, default=None,
                    help='Gradient directions (.bvec file).'
                    )
    p.add_argument('-d','--delta', action='store', metavar='delta', dest='delta',
                    type=float, required=False, default=None,
                    help='Pulse duration time in seconds.'
                    )
    p.add_argument('-s','--Delta', action='store', metavar='Delta', dest='separation',
                    type=float, required=False, default=None,
                    help='Pulse separation time in seconds.'
                    )
    p.add_argument('-o', '--output', action='store', metavar='output', dest='output',
                    type=str, required=True,
                    help='Output folder.'
                    )
    return p   
       
def main():
    parser = argumentParser()
    args = parser.parse_args()
    info = DwiInfo(args)	
   
    if (os.path.exists(args.output) == False):
        os.makedirs(args.output) 
        print("The", args.output, " is created!")
       
    if (os.path.exists(args.output + '//' + 'DWIs') == False):
        os.makedirs(args.output+ '//' + 'DWIs') 
    if (os.path.exists(args.output + '//' + 'Masks') == False):
        os.makedirs(args.output+ '//' + 'Masks') 

    dwi, affine = info.loadDWI()
    mask, affine = info.getMask()
    for i in range (0,mask.shape[1]):
        if np.any(mask[:,i:i+1,:]):
            img = nib.Nifti1Image(dwi[:,i:i+1,:,:], affine)
            nib.save(img,args.output+ '//' + 'DWIs'+r'/DWI'+str(i)+'.nii.gz')
            img2 = nib.Nifti1Image(mask[:,i:i+1,:], affine)
            nib.save(img2,args.output+ '//' + 'Masks'+r'/Mask'+str(i)+'.nii.gz')
       
    with open(args.output+'/pickle.pkl', 'wb') as f: 
        pickle.dump(info.responseFunction(), f)   
        
    if (os.path.exists(args.output + '//' + 'results') == False):
        os.makedirs(args.output+ '//' + 'results') 
   
if __name__ == '__main__':
    main()