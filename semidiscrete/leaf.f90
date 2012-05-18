subroutine printleaf_ad(wavelength,drdkab,drdkar,drdkw,drdm,dtdkab,dtdkar,dtdkw,dtdm,drdn,dtdn)
      use dataspec_p5
      double precision   :: dr1(nl),dt1(nl),drn(nl),dtn(nl),wavelength(nl)
      double precision   :: drdkab(nl), drdkar(nl), drdkw(nl), drdm(nl)
      double precision   :: dtdkab(nl), dtdkar(nl), dtdkw(nl), dtdm(nl),drdn(nl),dtdn(nl)
      do i=1,nl
        print*,wavelength(i),drdkab(i), drdkar(i), drdkw(i), drdm(i),drdn(i),dtdkab(i), dtdkar(i), dtdkw(i), dtdm(i),dtdn(i)
      enddo

end

subroutine printleaf(wavelength,refl, trans)
      use dataspec_p5
      double precision   :: kk(nlmax),k(nl)
      double precision   :: wavelength(nl)
      double precision   :: refl(nl)
      double precision   :: trans(nl)

      common /kleaf/ kk

      k = kk(1:nl)
      do i=1,nl
        print*,wavelength(i),refl(i),trans(i),k(i)
      enddo
end

subroutine leaf(wavelength,refl, trans, kab,kar,kw,km,n)
      use dataspec_p5
      double precision   :: wavelength(nl)
      double precision   :: refl(nl)
      double precision   :: trans(nl)
      double precision   :: kk(nlmax),k1(nl),k2(nl),k3(nl)
      double precision   :: transmission(nl)
      double precision   :: rek(nlmax),rik(nlmax)
      double precision   :: re(nl),ri(nl),k(nl)
      double precision   :: tn(nl),w(nl),delta(nl)
      double precision   :: ptheta(nl)
      double precision   :: kab,kar,kw,km,n
      common /leaf_params/ rek,rik
      common /kleaf/ kk

	k = kk(1:nl)
	re = rek(1:nl)
	ri = rik(1:nl)
! this is an approximation to the transmission, theta
      ptheta = 0.354824 + 0.0315842*(nrefrac-1.35)
      wavelength = lambda
      k1 = 1.194494
      k2 =  -0.0343164638713112* (nrefrac-1.35)
      k3 =  - .0110496402235069*((nrefrac-1.35)**2)
      k = (acab*kab + acar*kar + acw*kw + acm*km)*(k1+k2+k3)
      transmission = (1 - ptheta)*exp(-k)/(1 - ptheta*exp(-k))
! end of this approximation
      tn = (n*ri)*(transmission**(2))
      w = (1 - ri) * transmission / (1 - ri * transmission)
      delta = (tn) / (transmission**n + tn)
      refl = re + (1 - re) *w * delta
      trans = (1 - re) *w *(1 - delta )  
end

subroutine leaf_ad_anc(kab,kar,kw,km,n)
      use dataspec_p5
      double precision   :: rek(nlmax),rik(nlmax),kk(nlmax),e1(nl),n1(nl),e2(nl),e3(nl),e4(nl),e5(nl),e6(nl),e7(nl)
      double precision   :: mult(nl),k1(nl),k2(nl),k3(nl),dr1k(nlmax),dt1k(nlmax),drnk(nlmax),dtnk(nlmax)
      double precision   :: kab,kar,kw,km,n
      double precision   :: ptheta(nl),p2(nl),p3(nl),re(nl),ri(nl),k(nl)
      common /leaf_params/ rek,rik
      common /leaf_advars/ dr1k,dt1k,drnk,dtnk

	re = rek(1:nl)
	ri = rik(1:nl)
	k = kk(1:nl)

      ptheta = 0.354824 + 0.0315842*(nrefrac-1.35)
      k1 = 1.194494
      k2 =  -0.0343164638713112* (nrefrac-1.35)
      k3 =  - .0110496402235069*((nrefrac-1.35)**2)
      mult = (k1+k2+k3)
      k = (acab*kab + acar*kar + acw*kw + acm*km)*(k1+k2+k3)
      e1 = exp(k)
      p2 = (1-ptheta)*(1-ptheta)
      p3 = -p2*(1-ptheta)
      n1 = ((-1+ptheta)/(ptheta-e1))**n
      e2 = (e1 + ptheta*(-1 + ri) - ri)**2
      e3 = (e1 - ptheta)**2
      e4 = (e3*n1+n*p2*ri)**2
      e5 = (e1 - ptheta)**3
      e6 = (1 + (n*p2*ri)/ (e3* n1))**2
      e7 = log((-1 + ptheta)/ (-e1 + ptheta))
      dr1k(1:nl) = mult * (e1*n*p3*(re-1)*(ri-1)*ri*(n*p2*ri-(e1-ptheta)*n1*(e1*(n-3)-(n-3)*ptheta+(n-2)*(ptheta-1)*ri)))/(e2*e4)
      dt1k(1:nl) = mult * (e1*(ptheta-1)*(re-1)*(ri-1)*(1+(n*p2*ri*(e1*(n-1)+ptheta-n*ptheta+(n-2)*(ptheta-1)*ri))/(e5*n1)))/(e2*e6)
      drnk(1:nl) =   (e3* p3*n1* (-1 + re)*(-1 + ri)*ri*(-1 + n*e7))/((e1 + ptheta*(-1 + ri) - ri)* e4)
      dtnk(1:nl) =  - ((e3* p3*n1* (-1 + re)*(-1 + ri)*ri*(-1 + n*e7))/ ((e1 + ptheta*(-1 + ri) - ri)*e4))
end

subroutine leaf_ad(wavelength,drdkab,drdkar,drdkw,drdm,dtdkab,dtdkar,dtdkw,dtdm,drdn,dtdn,kab,kar,kw,km,n)
      use dataspec_p5
      implicit none
      double precision   :: dr1(nlmax),dt1(nl),drn(nl),dtn(nl),wavelength(nl)
      double precision   :: drdkab(nl), drdkar(nl), drdkw(nl), drdm(nl)
      double precision   :: dtdkab(nl), dtdkar(nl), dtdkw(nl), dtdm(nl),drdn(nl),dtdn(nl)
      double precision   :: kab,kar,kw,km,n
      double precision   :: refl(nl),trans(nl),refl1(nl),trans1(nl)
      double precision   :: dr1k(nlmax),dt1k(nlmax),drnk(nlmax),dtnk(nlmax)
      common /leaf_advars/ dr1k,dt1k,drnk,dtnk

      wavelength = lambda
      call leaf_ad_anc(kab,kar,kw,km,n)
      drdkab = acab * dr1k(1:nl)
      drdkar = acar * dr1k(1:nl)
      drdkw = acw * dr1k(1:nl)
      drdm = acm * dr1k(1:nl)
      dtdkab  = acab * dt1k(1:nl)
      dtdkar = acar * dt1k(1:nl)
      dtdkw = acw * dt1k(1:nl)
      dtdm = acm * dt1k(1:nl)
      drdn = drnk(1:nl)
      dtdn = dtnk(1:nl)
end


subroutine leafprep(nbands_to_use,bands_to_use)
      use dataspec_p5
      implicit none

      integer nbands_to_use
      integer bands_to_use(nbands_to_use)

      double precision   :: rek(nlmax),rik(nlmax)
      double precision   :: t1(nbands_to_use),t2(nbands_to_use)
      double precision    :: ang1
      double precision    :: ang2
      common /leaf_params/ rek,rik

      nl = nbands_to_use
      call subsetspectra(nbands_to_use,bands_to_use)
      ang1=90.0
      ang2=40.0
      call tav_abs(t2,ang2,nrefrac)
      call tav_abs(t1,ang1,nrefrac)
      rek(1:nl) = 1.0 - t2
      rik(1:nl) = 1.0 - t1/(nrefrac*nrefrac)
       
end
