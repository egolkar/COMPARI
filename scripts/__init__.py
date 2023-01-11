from compari import main
from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'A package to model peritumoal context '


setup(
    name="convrsn",
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    author="Ehsan Golkar",
    author_email="egolkar@gmail.com",
    license='MIT',
    packages=find_packages(),
    install_requires=[],
    keywords='conversion',

)
