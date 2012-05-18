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
            "--f90exec=/usr/bin/gfortran44", 
            "--f77flags='-ffixed-form -ffixed-line-length-none'" ])

    name         = "semidiscrete"  # name of the generated python extension (.so)
    description  = "semidiscrete python wrappers"
    author       = "J Gomez-Dans/NCEO & University College London"
    author_email = "j.gomez-dans@ucl.ac.uk"
    url="http://fapar.jrc.ec.europa.eu/WWW/Data/Pages/FAPAR_Software/FAPAR_Software_RTModels_1-2Discrete.php"
    
    setup( name=name,\
        description=description, \
        author=author, \
        author_email = author_email, \
        url=url,
        configuration = configuration, version="1.0",\
        packages=["semidiscrete"])
    