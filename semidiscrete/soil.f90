subroutine soil(wavelength,srefl,s1,s2,s3,s4)
      use dataspec_p5
      implicit none
      integer :: notok(nl)
      double precision   :: wavelength(nl)
      double precision   :: srefl(nl),s1,s2,s3,s4

      srefl=s1*phis1+s2*phis2+s3*phis3+s4*phis4
end

subroutine soil_ad(wavelength,drefl_ds1,drefl_ds2,drefl_ds3,drefl_ds4,s1,s2,s3,s4)
      use dataspec_p5
      implicit none
      double precision   :: wavelength(nl)
      double precision   :: srefl(nl),s1,s2,s3,s4
      double precision   :: drefl_ds1(nl),drefl_ds2(nl),drefl_ds3(nl),drefl_ds4(nl)

      drefl_ds1 = phis1
      drefl_ds2 = phis2
      drefl_ds3 = phis3
      drefl_ds4 = phis4

end

