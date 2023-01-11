
from dmipy.core import modeling_framework
from dmipy.core.acquisition_scheme import acquisition_scheme_from_bvalues
from os.path import join
import numpy as np
from dipy.io.image import load_nifti
from dmipy.core.acquisition_scheme import gtab_dmipy2dipy
from dipy.reconst.csdeconv import auto_response_ssst
from dmipy.signal_models import gaussian_models,cylinder_models
from dmipy.core.modeling_framework import MultiCompartmentModel
from dmipy.core import modeling_framework
import nibabel as nib
import argparse
import os
import sys 

from dmipy.distributions import distribute_models
import nibabel as nib
from dmipy.distributions.distribute_models import SD2BinghamDistributed
import pickle
from dmipy.distributions.distribute_models import BundleModel
from dipy.data import get_sphere
from dipy.reconst import shm
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

DESCRIPTION = '''
COMPARI = Crossing fiber Modeling in the Peritumoral Area using clinically feasible dMRI
'''

class DwiInfo:
    def __init__(self,args):
        #delta = 0.01293,Delta = 0.03666
        self.dmri= args.dMRI_filename
        self.mask = args.mask
        self.bval= args.bvals
        self.bvec = args.bvecs
        self.delta = args.delta
        self.separation = args.separation
        self.lmax = args.lmax
        self.xflip = args.xflip
        self.yflip = args.yflip
        self.zflip = args.zflip
        self.output = args.output
        
        #COMPARI fixed parameters 
        self.freeWaterDiffusivity = 3.0e-9
        self.minDiffusivity = 0.1e-9
        
    
    def getScheme(self):
        acquisition_path = modeling_framework.GRADIENT_TABLES_PATH
        bvalues = np.loadtxt(self.bval)
        bvalues_SI = bvalues * 1e6  # SI units as s/m^2
        gradient_directions = np.loadtxt(self.bvec)
        gradient_directions[0] *= float(self.xflip) # flip the X coordinates, based on your data
        gradient_directions[1] *= float(self.yflip) # flip the Y coordinates, based on your data
        gradient_directions[2] *= float(self.zflip) # flip the Z coordinates, based on your data
        return (acquisition_scheme_from_bvalues(bvalues_SI, gradient_directions.T, self.delta, self.separation))
        
    def loadNifty(self):
        return (load_nifti(self.dmri))
        
    def getMask(self):
        return (load_nifti(self.mask))
    
    def getFreeWaterDiffusivity(self):
        return self.freeWaterDiffusivity
        
    def getMinDiffusivity(self):
        return self.minDiffusivity
        
    def getlmax(self):
        return int(self.lmax)
        
    def isotropicFit(self):
        scheme = self.getScheme()
        scheme.print_acquisition_info
        dwi, affine = self.loadNifty()
        dmask, aff  = self.getMask()
        
        ball = gaussian_models.G1Ball()

        my_mod = modeling_framework.MultiCompartmentModel([ball])
        my_mod.set_parameter_optimization_bounds('G1Ball_1_lambda_iso', [self.minDiffusivity, self.freeWaterDiffusivity])

        mymod_fit = my_mod.fit(scheme, dwi, mask = dmask.any()>0)
        iso = mymod_fit.fitted_parameters['G1Ball_1_lambda_iso']
      
        print ('The isotopic ball was fitted!')
        return iso
        
    def COMPARI(self):
        scheme = self.getScheme()
        dwi, affine = self.loadNifty()  
        dmask, aff  = self.getMask()
        iso = self.isotropicFit()
        sphere = get_sphere(name='symmetric724')
        
        gtab = gtab_dmipy2dipy(scheme)
        response, ratio = auto_response_ssst(gtab, dwi, roi_center=None, roi_radii=10, fa_thr=0.7)
        lambdas = response[0]
        lambda_par = lambdas[0] * 1e-6 #Lambda_parallel   
        print('The lambda_par is :',lambda_par)
        
        edemaIndex = np.where((iso >lambda_par) & (iso<2.7e-9))
        nonEdemaIndex = np.where(iso <= lambda_par  )
        
        ###########################################################
        dwiShape = dwi.shape
        partial_volume_0 = np.zeros(dwiShape[0:3], dtype=np.float)
        partial_volume_1 = np.zeros(dwiShape[0:3], dtype=np.float)
        vf_FW = np.zeros(dwiShape[0:3], dtype=np.float)
        vf_HW = np.zeros(dwiShape[0:3], dtype=np.float)
        vf_stick = np.zeros(dwiShape[0:3], dtype=np.float)
        vf_zepp = np.zeros(dwiShape[0:3], dtype=np.float)
        fod = np.zeros((dwiShape[0], dwiShape[1], dwiShape[2], int((self.lmax+1)*(self.lmax+2)/2) ), dtype=np.float)   
        BundleModel_1_partial_volume_0 =  np.zeros(dwiShape[0:3], dtype=np.float)
        BundleModel_2_partial_volume_0 =  np.zeros(dwiShape[0:3], dtype=np.float)
        AI = np.zeros(dwiShape[0:3], dtype=np.float)
        ###########################################################        
        ball_csf = gaussian_models.G1Ball()
        ball_hindred= gaussian_models.G1Ball()
        stick = cylinder_models.C1Stick()
        zeppelin = gaussian_models.G2Zeppelin()
        ################## Bundle_aniso ###################
        Bundle_aniso = distribute_models.BundleModel([stick,zeppelin])
        Bundle_aniso.set_equal_parameter('C1Stick_1_lambda_par', 'G2Zeppelin_1_lambda_par') 
        Bundle_aniso.set_tortuous_parameter('G2Zeppelin_1_lambda_perp','C1Stick_1_lambda_par','partial_volume_0')
        Bundle_aniso.set_fixed_parameter('C1Stick_1_lambda_par',lambda_par)

        
        if (len(edemaIndex[0])>0):
            print('fitting to Edema voxels ' )
            listIndex = list(zip(edemaIndex[0],edemaIndex[1],edemaIndex[2]))
            ################## Bundle_iso ###################  
            for i in listIndex: 
                bundle_iso = distribute_models.BundleModel([ball_csf,ball_hindred])            
                bundle_iso.set_fixed_parameter('G1Ball_1_lambda_iso', self.freeWaterDiffusivity )
                bundle_iso.set_fixed_parameter('G1Ball_2_lambda_iso', float(iso[i]) )
                my_mod = modeling_framework.MultiCompartmentModel([bundle_iso,Bundle_aniso])
                mymod_fit = my_mod.fit(scheme,dwi[i], mask = dmask.any()>0)  
                #vf_edema[i, slc, j] = 1
                partial_volume_0[i] = mymod_fit.fitted_parameters['partial_volume_0']
                partial_volume_1[i] = mymod_fit.fitted_parameters['partial_volume_1']
                BundleModel_1_partial_volume_0[i] = mymod_fit.fitted_parameters['BundleModel_1_partial_volume_0']
                BundleModel_2_partial_volume_0[i] = mymod_fit.fitted_parameters['BundleModel_2_partial_volume_0']
                
                SH_mod = modeling_framework.MultiCompartmentSphericalHarmonicsModel(models=[bundle_iso, Bundle_aniso], sh_order=self.lmax)           
                SH_mod.set_fixed_parameter('partial_volume_0', mymod_fit.fitted_parameters['partial_volume_0']) 
                SH_mod.set_fixed_parameter('partial_volume_1', mymod_fit.fitted_parameters['partial_volume_1']) 
                SH_mod.set_fixed_parameter('BundleModel_1_partial_volume_0', mymod_fit.fitted_parameters['BundleModel_1_partial_volume_0'])
                SH_mod.set_fixed_parameter('BundleModel_2_partial_volume_0', mymod_fit.fitted_parameters['BundleModel_2_partial_volume_0'])
                sh_mod_fit = SH_mod.fit(scheme,dwi[i], mask = dmask.any()>0,solver='csd')
                
                vf_FW[i]= BundleModel_1_partial_volume_0[i] * partial_volume_0[i]
                vf_HW[i]= (1-BundleModel_1_partial_volume_0[i]) * partial_volume_0[i]
                
                vf_stick[i]= BundleModel_2_partial_volume_0[i] * partial_volume_1[i]
                vf_zepp[i]= (1-BundleModel_2_partial_volume_0[i]) * partial_volume_1[i]
            
                fod[i,:] = sh_mod_fit.fod_sh()[:]
                AI[i] = sh_mod_fit.anisotropy_index()
        
            
        ############### fitting to nonEdema voxels ####################
        if (len(nonEdemaIndex[0])>0):
            print('fitting to nonEdema voxels ') 
            my_mod = modeling_framework.MultiCompartmentModel([ball_csf,Bundle_aniso])
            my_mod.set_fixed_parameter('G1Ball_1_lambda_iso',self.freeWaterDiffusivity)
            
            mymod_fit = my_mod.fit(scheme,dwi[nonEdemaIndex], mask = dmask.any()>0)  
            partial_volume_0[nonEdemaIndex] = mymod_fit.fitted_parameters['partial_volume_0']
            partial_volume_1[nonEdemaIndex] = mymod_fit.fitted_parameters['partial_volume_1']
            #aniso 
            BundleModel_1_partial_volume_aniso = mymod_fit.fitted_parameters['BundleModel_1_partial_volume_0']                
            lambda_perp = (1 - BundleModel_1_partial_volume_aniso) * lambda_par

            SH_mod = modeling_framework.MultiCompartmentSphericalHarmonicsModel(models=[ball_csf, Bundle_aniso], sh_order= self.lmax)
            SH_mod.set_fixed_parameter('G1Ball_1_lambda_iso',self.freeWaterDiffusivity)
            SH_mod.set_fixed_parameter('partial_volume_0', mymod_fit.fitted_parameters['partial_volume_0']) 
            SH_mod.set_fixed_parameter('partial_volume_1', mymod_fit.fitted_parameters['partial_volume_1']) 
            SH_mod.set_fixed_parameter('BundleModel_1_partial_volume_0', BundleModel_1_partial_volume_aniso)
            sh_mod_fit = SH_mod.fit(scheme,dwi[nonEdemaIndex], mask = dmask.any()>0,solver='csd')

            BundleModel_1_partial_volume_0[nonEdemaIndex] = mymod_fit.fitted_parameters['BundleModel_1_partial_volume_0']
            BundleModel_2_partial_volume_0[nonEdemaIndex] =  (1 - mymod_fit.fitted_parameters['BundleModel_1_partial_volume_0']) * lambda_par  
            vf_FW[nonEdemaIndex] = mymod_fit.fitted_parameters['partial_volume_0']


            vf_stick[nonEdemaIndex]= BundleModel_1_partial_volume_0[nonEdemaIndex] * partial_volume_1[nonEdemaIndex]
            vf_zepp[nonEdemaIndex]= (1-BundleModel_1_partial_volume_0[nonEdemaIndex]) * partial_volume_1[nonEdemaIndex]
            
            listIndex = list(zip(nonEdemaIndex[0],nonEdemaIndex[1],nonEdemaIndex[2]))
            index =0
            for i in listIndex:
                fod[i[0],i[1],i[2],:] = sh_mod_fit.fod_sh()[index,:]
                index +=1
            
            AI[nonEdemaIndex] = sh_mod_fit.anisotropy_index()
        
        sf = shm.sh_to_sf(fod,sphere,sh_order=self.lmax,basis_type = 'tournier07')
        fods = shm.sf_to_sh(sf,sphere,sh_order=self.lmax,basis_type = 'tournier07')
        ############################################################
        
        
        return partial_volume_0,partial_volume_1,vf_FW,vf_HW,vf_stick,vf_zepp,fods,AI,aff 

        
        

def argumentParser():
    p = argparse.ArgumentParser(description=DESCRIPTION)
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
                    type=float, required=True, default=None,
                    help='Pulse duration time in seconds.'
                    )
    p.add_argument('-s','--Delta', action='store', metavar='Delta', dest='separation',
                    type=float, required=True, default=None,
                    help='Pulse separation time in seconds.'
                    )
    p.add_argument('-l','--lmax', action='store', metavar='lmax', dest='lmax',
                    type=int, required=True, default=8,
                    help='Harmonic order. Defalut is 8.'
                    )      
    p.add_argument('-fx','--xflip',action='store',metavar='flip',dest='xflip',
                    type=str, required=False, default=1, 
                    help='flip gradient directions through X axis. Defaluit is 1 otherwise use -1'
                    )
    p.add_argument('-fy','--yflip',action='store',metavar='flip',dest='yflip',
                    type=str, required=False, default=1, 
                    help='flip gradient directions through Y axis. Defaluit is 1 otherwise use -1'
                    )
    p.add_argument('-fz','--zflip',action='store',metavar='flip',dest='zflip',
                    type=str, required=False, default=1, 
                    help='flip gradient directions through Z axis. Defaluit is 1 otherwise use -1'
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
    partial_volume_0,partial_volume_1,vf_FW,vf_HW,vf_stick,vf_zepp,fods,AI,aff = info.COMPARI()
    
    img= nib.Nifti1Image(partial_volume_0, aff)
    nib.save(img,args.output+'/vf_bundle_iso.nii.gz' )

    img= nib.Nifti1Image(partial_volume_1, aff)
    nib.save(img,args.output+'/vf_bundle_aniso.nii.gz' )
    
    img= nib.Nifti1Image(vf_FW, aff)
    nib.save(img,args.output+'/vf_FW.nii.gz' )
    
    img= nib.Nifti1Image(vf_HW, aff)
    nib.save(img,args.output+'/vf_HW.nii.gz' )

    img= nib.Nifti1Image(vf_stick, aff)
    nib.save(img,args.output+'/vf_intra.nii.gz' )

    img= nib.Nifti1Image(vf_zepp, aff)
    nib.save(img,args.output+'/vf_extra.nii.gz' )


    img= nib.Nifti1Image(fods, aff)
    nib.save(img,args.output+'/fods.nii.gz' )

    img= nib.Nifti1Image(AI, aff)
    nib.save(img,args.output+'/AI.nii.gz' )
    
if __name__ == '__main__':
    main()
    