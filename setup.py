from numpy.distutils.core import setup

import sys, os

version = '1.0'
from numpy.distutils.core import Extension

ext1 = Extension ( 'rtmodel_ad_trans2', \
    sources = ['semidiscrete/rtmodel_ad_trans2.pyf', \
    "semidiscrete/dataSpec_P5.f90",  "semidiscrete/leaf.f90", \
    "semidiscrete/msnadimbrf_ad.f90", "semidiscrete/nadimmod_ad.f90", \
    "semidiscrete/rtmodel_ad_trans.f90", "semidiscrete/soil.f90", \
    "semidiscrete/tav_abs.f90", "semidiscrete/nadimbrf_ad.f",  \
    "semidiscrete/nadimtools_ad.f" ] )

if __name__ == "__main__":
    from numpy.distutils.core import setup
    sys.argv.extend ( ["config_fc", "--fcompiler=gnu95", 
#        "--f90exec=/usr/bin/gfortran44", 
#        "--f77exec=/usr/bin/gfortran44", 
        "--f77flags='-ffixed-form -ffixed-line-length-none -fdefault-real-8 -O3 -march=native'"] )
    #"--f77flags='-ffree-form -ffixed-line-length-none -fdefault-real-8 '"])
    #"--f90flags='-ffree-form -ffixed-line-length-none -fdefault-real-8 '" ])
    
    DISTNAME = 'semidiscrete1'
    DESCRIPTION = 'SemiDiscrete python wrappers'
    LONG_DESCRIPTION = open('README.txt').read()
    MAINTAINER = 'Jose Gomez-Dans/NCEO & University College London'
    MAINTAINER_EMAIL = "j.gomez-dans@ucl.ac.uk"
    URL = 'http://github.com/jgomezdans/semidiscrete'
    LICENSE = 'Undecided'
    VERSION = "1.0.2"
    DOWNLOAD_URL="https://github.com/jgomezdans/semidiscrete/zipball/master"
    
    setup ( name=DISTNAME,
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    version=VERSION,
    download_url=DOWNLOAD_URL,
    long_description=LONG_DESCRIPTION,
    packages=['semidiscrete1'],
    package_dir={'semidiscrete1': '.'},
    ext_modules = [ext1,],
    ext_package="semidiscrete1",
    classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved',
        'Programming Language :: Fortran',
        'Programming Language :: Python',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS'
        ]
    )
