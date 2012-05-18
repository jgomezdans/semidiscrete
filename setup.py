from numpy.distutils.core import setup

import sys, os

version = '1.0'

def configuration ( parent_package='', top_path=None ):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('semidiscrete',parent_package,top_path, )
    
    config.add_extension ('rtmodel_ad_trans2', 
    sources = ['semidiscrete/rtmodel_ad_trans2.pyf', \
        "semidiscrete/dataSpec_P5.f90",  "semidiscrete/leaf.f90", \
        "semidiscrete/msnadimbrf_ad.f90", "semidiscrete/nadimmod_ad.f90", \
        "semidiscrete/rtmodel_ad_trans.f90", "semidiscrete/soil.f90", \
        "semidiscrete/tav_abs.f90", "semidiscrete/nadimbrf_ad.f",  \
        "semidiscrete/nadimtools_ad.f" ] )
    return config

if __name__ == "__main__":    
    # This is a hack, as F77/F90 flags can be passed in the config object
    # but only on recent versions of numpy/f2py
    sys.argv.extend ( ["config_fc", "--fcompiler=gnu95", 
            #"--f90exec=/usr/bin/gfortran44", 
            "--f77flags='-ffixed-form -ffixed-line-length-none'" ])

    DISTNAME = 'semidiscrete'
    DESCRIPTION = 'semidiscrete python wrappers'
    LONG_DESCRIPTION = open('README.rst').read()
    MAINTAINER = 'Jose Gomez-Dans/NCEO & University College London'
    MAINTAINER_EMAIL = "j.gomez-dans@ucl.ac.uk"
    URL = 'http://github.com/jgomezdans/semidiscrete'
    LICENSE = 'Undecided'
    VERSION = "1.0.1"
    DOWNLOAD_URL="https://github.com/jgomezdans/semidiscrete/zipball/master"
    setup ( configuration = configuration,
        name=DISTNAME,
        maintainer=MAINTAINER,
        include_package_data=True,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        long_description=LONG_DESCRIPTION,
        zip_safe=False, # the package can run out of an .egg file
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
            ], packages=["semidiscrete"]
    )
       
    