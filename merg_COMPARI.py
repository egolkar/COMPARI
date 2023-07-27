# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 21:47:51 2023

@author: GolkarE
"""

import nibabel as nib
import numpy as np
import os

from dipy.data import get_sphere
from dipy.reconst import shm
from dipy.direction import peak_directions
import sys 


vol= str(sys.argv[1])

#vol = '300-800'
#ext='/noflip'

#ext= '/'+str(sys.argv[1])

path = str(sys.argv[2]) 
# '/cbica/home/golkare/Data/Multishell3-062'

mask = str(sys.argv[3]) # mask file

lmax = int(sys.argv[4])
sphere = get_sphere(name='symmetric724')




nii = nib.load(mask)
img = nii.get_fdata()
ii= np.where(img==0)

vf_bundle_aniso = np.zeros(img.shape, dtype=np.float)
vf_bundle_iso = np.zeros(img.shape, dtype=np.float)
vf_FW = np.zeros(img.shape, dtype=np.float)
vf_HW = np.zeros(img.shape, dtype=np.float)
vf_intra= np.zeros(img.shape, dtype=np.float)
vf_extra= np.zeros(img.shape, dtype=np.float)
fods = np.zeros((img.shape[0], img.shape[1], img.shape[2], int((lmax+1)*(lmax+2)/2) ), dtype=np.float)
AI= np.zeros(img.shape, dtype=np.float)


#print('img.shape[1]',img.shape[1])

for i in range(1,img.shape[1],1):
    
    if (os.path.exists(path+'/'+vol + '/results/'+str(i)+'/vf_bundle_aniso.nii.gz') == False):
        continue
    print('vol:',vol,' slice i:',i)
    nii = nib.load(path+'/'+vol + '/results/' +str(i)+'/vf_bundle_aniso.nii.gz')   
    bundle_aniso = nii.get_fdata()
    vf_bundle_aniso[:,i:i+1,:]=bundle_aniso
    
    nii = nib.load(path+'/'+vol + '/results/' +str(i)+'/vf_bundle_iso.nii.gz')   
    bundle_iso = nii.get_fdata()
    vf_bundle_iso[:,i:i+1,:]=bundle_iso
    
    nii = nib.load(path+'/'+vol + '/results/' +str(i)+'/vf_FW.nii.gz')   
    FW = nii.get_fdata()
    vf_FW[:,i:i+1,:]=FW
    
    nii = nib.load(path+'/'+vol + '/results/' +str(i)+'/vf_HW.nii.gz')   
    HW = nii.get_fdata()
    vf_HW[:,i:i+1,:]=HW
    
    nii = nib.load(path+'/'+vol + '/results/' +str(i)+'/vf_intra.nii.gz')   
    intra = nii.get_fdata()
    vf_intra[:,i:i+1,:]=intra
    
    nii = nib.load(path+'/'+vol + '/results/' +str(i)+'/vf_extra.nii.gz')   
    extra = nii.get_fdata()
    vf_extra[:,i:i+1,:]=extra
    
    nii = nib.load(path+'/'+vol + '/results/' +str(i)+'/fods.nii.gz')   
    fod = nii.get_fdata()
    fods[:,i:i+1,:,:]=fod
    
    nii = nib.load(path+'/'+vol + '/results/' +str(i)+'/AI.nii.gz')   
    ai = nii.get_fdata()
    AI[:,i:i+1,:]=ai

vf_bundle_aniso[ii] = 0     
img_vf_bundle_aniso = nib.Nifti1Image(vf_bundle_aniso, nii.affine)
nib.save(img_vf_bundle_aniso,path+'/'+vol+'/results' +'/vf_bundle_aniso.nii.gz')

vf_bundle_iso[ii]=0
img_vf_bundle_iso = nib.Nifti1Image(vf_bundle_iso, nii.affine)
nib.save(img_vf_bundle_iso,path+'/'+vol+'/results/vf_bundle_iso.nii.gz')

vf_FW[ii]=0
img_vf_FW = nib.Nifti1Image(vf_FW, nii.affine)
nib.save(img_vf_FW,path+'/'+vol+'/results' +'/vf_FW.nii.gz')

vf_HW[ii]=0
img_vf_HW = nib.Nifti1Image(vf_HW, nii.affine)
nib.save(img_vf_HW,path+'/'+vol+'/results' +'/vf_HW.nii.gz')

vf_intra[ii]=0
img_vf_intra = nib.Nifti1Image(vf_intra, nii.affine)
nib.save(img_vf_intra,path+'/'+vol+'/results' +'/vf_intra.nii.gz')

vf_extra[ii]=0
img_vf_extra = nib.Nifti1Image(vf_extra, nii.affine)
nib.save(img_vf_extra,path+'/'+vol+'/results' +'/vf_extra.nii.gz')

AI[ii]=0
img_AI = nib.Nifti1Image(AI, nii.affine)
nib.save(img_AI,path+'/'+vol+'/results' +'/AI.nii.gz')


freeWater = np.where(vf_bundle_iso>.9)
listIndex = list(zip(freeWater[0],freeWater[1],freeWater[2]))
index =0
for i in listIndex:
    fods[i[0],i[1],i[2],:] = np.zeros((1, 1, 1, int((lmax+1)*(lmax+2)/2) ), dtype=np.float)
    index +=1


wmfod= np.zeros(fods.shape)
sf = shm.sh_to_sf(fods,sphere,sh_order=lmax,basis_type = 'tournier07')
for i in range(1,img.shape[0],1):
    print ( 'vol: ', vol, ' i:',i )
    for j in range(1,img.shape[1],1):
        for k in range(1,img.shape[2],1):               
                ps = peak_directions(sf[i,j,k,:], sphere, relative_peak_threshold=.25)
                if (len(ps[1]) != 0 and ps[1][0]<3.5):
                    wmfod[i,j,k,:] = fods[i,j,k,:]
                if (img[i,j,k] == 0 ):
                    wmfod[i,j,k,:] = np.zeros((1, 1, 1, int((lmax+1)*(lmax+2)/2) ), dtype=np.float)
          


img_wmfod = nib.Nifti1Image(wmfod, nii.affine) 
nib.save(img_wmfod,path+'/'+vol+'/results' +'/fods_lmax'+str(lmax)+'.nii.gz')
     








