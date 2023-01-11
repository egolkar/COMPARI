from distutils.core import setup
from distutils.extension import Extension
import numpy as np
import os


DESCRIPTION = 'A package to model crossing fiber in the peritumoral refion from dMRI' 

setup(name='compari',
    description=DESCRIPTION,
    version='0.1.dev0',
    author='Ehsan Golkar',
    author_email='egolkar@gmail.com',
    scripts = [ os.path.join('scripts', 'compari.py'),
                ],
    )
