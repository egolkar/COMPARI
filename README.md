# COMPARI

This Python package contains COMPARI, a tool for Crossing Fiber Modeling in the Peritumoral Area using clinically feasible dMRI 

COMPARI uses NIfTI as its data format for MR and brain mask images, using the "nii.gz" file extension. It uses the FSL convention for "bval" and "bvec" text files.

# Dependencies
+ Python >= 2.7
+ NumPy >= 1.10.4
+ Scipy >= 0.17.1
+ Nibabel >= 2.0.1
+ Dipy >= 0.10.1
+ [Dmipy](https://github.com/AthenaEPI/dmipy)
# Installation instruction

+ local install (user-wise, does not require root priviledge)

          run "python setup.py install --prefix=~/local"
          add "~/local/lib/python2.7/site-packages/" to your PYTHONPAT
          add "~/local/bin" to your PATH

+ system install (requires root privilege)
          run "sudo python setup.py install"
        
+ Input arguments

          -h, --help    Show this help message and exit
          -i  --data    Input dMRI data file (Nifti format).
          -m  --mask    Brain mask file (Nifti format).
          -b  --bvals   B-values given in s/mm^2 (.bval file).
          -r  --bvecs   Gradient directions (.bvec file).
          -d  --delta   Pulse duration time in seconds.
          -s  --Delta   Pulse separation time in seconds.
          -l  --lmax    Harmonic order. Defalut is 8.
          -o  --output  Output folder 
  
+ Optional arguments 

          -fx --xflip   Flip gradient directions through X axis. Defaluit is 1 otherwise use -1
          -fy --yflip   Flip gradient directions through Y axis. Defaluit is 1 otherwise use -1
          -fz --zflip   Flip gradient directions through Z axis. Defaluit is 1otherwise use -1





