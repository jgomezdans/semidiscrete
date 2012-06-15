subroutine msnadimbrf(nbands,nl,nv,lad,kab,kar,kw,km,n,s1,s2,s3,s4,theta_i,phi_i,theta_v,phi_v,xhc,xlai,rpl,brf,bandpass)
! multispectral version of nadimbrf
!      use dataspec_p5
      use  mo_nad
      implicit none
! passed terms
      integer :: lad,nv,nl
      integer :: nbands
      double precision :: kab,kar,kw,km,n,s1,s2,s3,s4
      double precision :: theta_v(nv),phi_v(nv),theta_i,phi_i !variable dimension
      double precision :: xlai,xhc,rpl
      double precision :: brf(nv,nbands) !changing to variable dimension
      double precision :: bandpass(nbands,nl)

! locals
      double precision :: xrl(nl),xtl(nl),xrs(nl)
      double precision :: wavelength(nl)
      double precision :: lrefl(nl)
      double precision :: ltrans(nl)
      double precision :: srefl(nl)
      integer :: j,thiswavelength
      integer :: k,l0,l1
      integer :: wavesamples(nl)
      double precision :: thisbandpass(nl)
      integer :: thisnbands
      double precision :: allbrf(nv,nl),somebrf(nl),result(nl),field(nl)
! initialise locals to 0
      wavesamples = 0
      xrl=0.
      xtl=0.
      xrs=0.
      wavelength=0.
      lrefl=0.
      ltrans=0.
      srefl=0.
      thiswavelength=0
      l0=0
      l1=0

! we have an int array of size [nl] called wavesamples
! this contains 1 if we need to simulate that wavelength and 0 otherwise
! The (normalised) bandpass functions are bandpass[nbands,nl]
! which are weightings 

! so first, sort out wavesamples as a mask
      do j=1,nbands
	thisbandpass(:) = bandpass(j,:)
        where(thisbandpass .gt. 0) wavesamples = 1
      enddo
      thisnbands = sum(wavesamples)
! when we come to call the rt model, we only need to process
! those wavebands where wavesamples == 1
! start the code
! The output, brf is of dimension [nv,nbands]
! but we have to use thisnbands
! for calculating the brf of dimensions allbrf[nv,thisnbands]

! the leaf & soil models work with all wavebands
      call leaf(wavelength,lrefl,ltrans, kab,kar,kw,km,n)
      call soil(wavelength,srefl,s1,s2,s3,s4)

! now pack up xrs, xrl,xtl
      xrs(1:thisnbands) = pack(srefl,wavesamples /= 0)
      xrl(1:thisnbands) = pack(lrefl,wavesamples /= 0)
      xtl(1:thisnbands) = pack(ltrans,wavesamples /= 0)

      nw = thisnbands
      call nadimbrf(theta_i,phi_i,nv,theta_v,phi_v,&
           lad,xrs(1:thisnbands),xhc,xlai,rpl,xrl(1:thisnbands),xtl(1:thisnbands),allbrf(:,1:thisnbands))
! now unpack for each waveband from allbrf(1,1:nw)
      somebrf(1:nw) = allbrf(1,1:nw)
      do j=1,nbands
        thisbandpass(:) = bandpass(j,:)
        field = 0.
        result = unpack(somebrf(1:nw),wavesamples /= 0,field)
	! now multiply element by element by thisbandpass
        ! and sum up
        brf(1,j) = sum(thisbandpass * result)
      end do
end

subroutine msnadimbrfmd(ipt,nbands,nl,nv,lad,kab,kar,kw,km,n,s1,s2,s3,s4,theta_i,phi_i,theta_v,phi_v,xhc,xlai,rpl,brf,bandpass)
! multispectral version of nadimbrf
!      use dataspec_p5
      use  mo_nad
      implicit none
! information passed
      integer :: lad,nv,ipt,nbands,nl
      double precision :: brf(nv,nbands) !changing to variable dimension
      double precision :: theta_v(nv),phi_v(nv),theta_i,phi_i !variable dimension
      double precision :: xlai,xhc,rpl
      double precision :: kab,kar,kw,km,n,s1,s2,s3,s4
      double precision :: bandpass(nbands,nl)
! local variables
      double precision :: lrefl(nl)
      double precision :: ltrans(nl)
      double precision :: srefl(nl)
      integer :: j,thiswavelength,k,l0,l1
      double precision :: wavelength(nl)
      double precision xrl(nl),xtl(nl),xrs(nl)
      double precision :: thisbandpass(nl)
      integer :: thisnbands,wavesamples(nl)
      double precision :: allbrf(nv,nl),somebrf(nl),result(nl),field(nl)

! initialise local variables
      wavesamples = 0
      lrefl=0.
      ltrans=0.
      srefl=0.
      thiswavelength=0
      l0=0
      l1=0
      wavelength=0.
      xrl=0.
      xtl=0.
      xrs=0.
! start the code
! we have an int array of size [nl] called wavesamples
! this contains 1 if we need to simulate that wavelength and 0 otherwise
! The (normalised) bandpass functions are bandpass[nbands,nl]
! which are weightings

! so first, sort out wavesamples as a mask
      do j=1,nbands
        thisbandpass(:) = bandpass(j,:)
        where(thisbandpass .gt. 0) wavesamples = 1
      enddo
      thisnbands = sum(wavesamples)
! when we come to call the rt model, we only need to process
! those wavebands where wavesamples == 1
! start the code
! The output, brf is of dimension [nv,nbands]
! but we have to use thisnbands
! for calculating the brf of dimensions allbrf[nv,thisnbands]

      call leaf(wavelength,lrefl,ltrans, kab,kar,kw,km,n)
      call soil(wavelength,srefl,s1,s2,s3,s4)

! now pack up xrs, xrl,xtl
      xrs(1:thisnbands) = pack(srefl,wavesamples /= 0)
      xrl(1:thisnbands) = pack(lrefl,wavesamples /= 0)
      xtl(1:thisnbands) = pack(ltrans,wavesamples /= 0)
      print*,'xrs',xrs(1:thisnbands)
      print*,'xrl',xrl(1:thisnbands)
      print*,'xtl',xtl(1:thisnbands)
      print*,'control',thisnbands,ipt,theta_i,phi_i,nv,theta_v,phi_v
      print*,'values',lad,xhc,xlai,rpl
      
      nw = thisnbands
      call nadimbrfmd(ipt,theta_i,phi_i,nv,theta_v,phi_v,lad,&
                      xrs(1:thisnbands),xhc,xlai,rpl,xrl(1:thisnbands),xtl(1:thisnbands),allbrf(:,1:thisnbands))
! now unpack for each waveband from allbrf(1,1:nw)
      somebrf(1:nw) = allbrf(1,1:nw)
      print*,'somebrf',somebrf
      do j=1,nbands
        thisbandpass(:) = bandpass(j,:)
        field = 0.
        result = unpack(somebrf(1:nw),wavesamples /= 0,field)
        ! now multiply element by element by thisbandpass
        ! and sum up
        brf(1,j)  = sum(thisbandpass * result)
      end do
end

subroutine msnadimbrf_ad(ipt,nbands,nl,nv,lad,kab,kar,kw,km,n,s1,s2,s3,s4,&
theta_i,phi_i,theta_v,phi_v,xhc,xlai,rpl,thisbrf_ad&
,xlai_ad,xhc_ad,rpl_ad,xkab_ad,xkar_ad,xkw_ad,xkm_ad,xleafn_ad,xs1_ad,xs2_ad,xs3_ad,xs4_ad,bandpass)
      ! multispectral version of nadimbrf
!      use dataspec_p5
      use  mo_nad
      implicit none
! passed variables
      integer          :: lad,nv,ipt,nbands,nl
      double precision :: thisbrf_ad(nv,nbands)
      double precision :: theta_v(nv),phi_v(nv),theta_i,phi_i !variable dimension
      double precision :: xlai_ad,xhc_ad,rpl_ad
      double precision :: xkab_ad,xkar_ad,xkw_ad,xkm_ad,xleafn_ad
      double precision :: xs1_ad,xs2_ad,xs3_ad,xs4_ad
      double precision :: xlai,xhc,rpl
      double precision :: kab,kar,kw,km,n,s1,s2,s3,s4
      double precision :: bandpass(nbands,nl)

! local variables
      double precision :: lrl_ad(nl),ltl_ad(nl),lrs_ad(nl)
      double precision :: lhc_ad(nl),llai_ad(nl),lrpl_ad(nl)
      double precision :: xrl(nl),xtl(nl),xrs(nl),brfpack(nl)
      double precision :: lrefl(nl)
      double precision :: ltrans(nl)
      double precision :: srefl(nl)
      double precision :: drl_kab_ad(nbands),drl_kar_ad(nbands),drl_kw_ad(nbands),drl_km_ad(nbands),drl_n_ad(nbands)
      double precision :: dtl_kab_ad(nbands),dtl_kar_ad(nbands),dtl_kw_ad(nbands),dtl_km_ad(nbands),dtl_n_ad(nbands)
      double precision :: dsl_s1_ad(nbands),dsl_s2_ad(nbands),dsl_s3_ad(nbands),dsl_s4_ad(nbands)
      integer :: j,k,thiswavelength,l0,l1
      double precision :: wavelength(nl),nwave
      double precision :: drdkab(nl), drdkar(nl), drdkw(nl), drdm(nl)
      double precision :: dtdkab(nl), dtdkar(nl), dtdkw(nl), dtdm(nl),drdn(nl),dtdn(nl)
      double precision :: drefl_ds1(nl),drefl_ds2(nl),drefl_ds3(nl),drefl_ds4(nl)
      double precision :: brf_ad(nv,nbands) !changing to variable dimension
      double precision :: thisbandpass(nl),totalwavesamples(nbands)
      integer :: thisnbands,wavesamples(nl),sumwavesamples(nl)
      double precision :: this,allbrf_ad(nv,nl),somebrf(nl),result(nl),field(nl)
      double precision :: rsresult(nl),rlresult(nl),tlresult(nl),hcresult(nl),lairesult(nl),rplresult(nl)
      double precision :: s1result(nl),s2result(nl),s3result(nl),s4result(nl)
      double precision :: kabresult(nl),karresult(nl),kwresult(nl),kmresult(nl),knresult(nl)

! inittialise local variables
      allbrf_ad = 0
      wavesamples = 0
      sumwavesamples = 0
      drl_kab_ad=0.
      drl_kar_ad=0.
      drl_kw_ad=0.
      drl_km_ad=0.
      drl_n_ad=0.
      dtl_kab_ad=0.
      dtl_kar_ad=0.
      dtl_kw_ad=0.
      dtl_km_ad=0.
      dtl_n_ad=0.
      dsl_s1_ad=0.
      dsl_s2_ad=0.
      dsl_s3_ad=0.
      dsl_s4_ad=0.
      thiswavelength=0
      l0=0
      l1=0
!      wavelength=0
!      drdkab=0
!      drdkar=0
!      drdkw=0
!      drdm=0
!      dtdkab=0
!      dtdkar=0
!      dtdkw=0
!      dtdm=0
!      drdn=0
!      dtdn=0
!      drefl_ds1=0
!      drefl_ds2=0
!      drefl_ds3=0
!      drefl_ds4=0
      brf_ad=0
      lrl_ad=0.
      ltl_ad=0.
      lrs_ad=0.
      lhc_ad=0.
      llai_ad=0.
      lrpl_ad=0.
      xrl=0.
      xtl=0.
      xrs=0.
      lrefl=0.
      ltrans=0.
      srefl=0.
! start the code 
! we have an int array of size [nl] called wavesamples
! this contains 1 if we need to simulate that wavelength and 0 otherwise
! The (normalised) bandpass functions are bandpass[nbands,nl]
! which are weightings

! so first, sort out wavesamples as a mask
      do j=1,nbands
        thisbandpass(:) = bandpass(j,:)
        where(thisbandpass .gt. 0) 
		where (wavesamples .ne. 1) sumwavesamples = 1
		wavesamples = 1
	end where
! NB allbrf_ad is *input* below , so we have to load it up with a weighted version of thisbrf_ad(nv,nbands)
        allbrf_ad(1,:) = allbrf_ad(1,:) + thisbrf_ad(1,j)  *  thisbandpass
        totalwavesamples(j) = sum(sumwavesamples)
	sumwavesamples=0
      enddo
      thisnbands = sum(wavesamples)
! when we come to call the rt model, we only need to process
! those wavebands where wavesamples == 1
! We achieve this using wavesamples as a mask
! start the code
! The output, brf is of dimension [nv,nbands]
! but we have to use thisnbands
! for calculating the brf of dimensions allbrf_ad[nv,thisnbands]

      call leaf(wavelength,lrefl,ltrans, kab,kar,kw,km,n)
      call soil(wavelength,srefl,s1,s2,s3,s4)
      call leaf_ad(wavelength,drdkab,drdkar,drdkw,drdm,dtdkab,dtdkar,dtdkw,dtdm,drdn,dtdn,kab,kar,kw,km,n)
      call soil_ad(wavelength,drefl_ds1,drefl_ds2,drefl_ds3,drefl_ds4,s1,s2,s3,s4)
! now pack up xrs, xrl,xtl & allbrf_ad
      xrs(1:thisnbands) = pack(srefl,wavesamples /= 0)
      xrl(1:thisnbands) = pack(lrefl,wavesamples /= 0)
      xtl(1:thisnbands) = pack(ltrans,wavesamples /= 0)
      brfpack(1:thisnbands) = pack(allbrf_ad(1,:),wavesamples /= 0)
!	print*,brfpack(1:thisnbands)
!      allbrf_ad(1,1:thisnbands) = brfpack(1:thisnbands)
      nw = thisnbands
! NB allbrf_ad is *input* here, so we have to load it up with a weighted version of thisbrf_ad(nv,nbands)
      call nadimbrf_ad(ipt,theta_i,phi_i,nv,theta_v,phi_v,lad,xrs(1:nw),&
             lrs_ad(1:nw),xhc,lhc_ad(1:nw),&
             xlai,llai_ad(1:nw),rpl,lrpl_ad(1:nw),xrl(1:nw),lrl_ad(1:nw),xtl(1:nw),ltl_ad(1:nw),brfpack(1:nw))
!  now do the partials
! what we get out is e.g. lrs_ad : dr/drs
! what we want, is a bandpass-weighted 
! e.g. xs1_ad : dr/ds1 = sum ( bandpass(l) * dr(l)/drs * drs(l)/ds1 )
! where sum(bandpass(l)) = 1
      field=0
      rsresult = unpack(lrs_ad(1:nw),wavesamples /= 0,field)
      field = 0.
      rlresult = unpack(lrl_ad(1:nw),wavesamples /= 0,field)
      field = 0.
      tlresult = unpack(ltl_ad(1:nw),wavesamples /= 0,field)
      field = 0.
      hcresult = unpack(lhc_ad(1:nw),wavesamples /= 0,field)
      field = 0.
      lairesult = unpack(llai_ad(1:nw),wavesamples /= 0,field)
      field = 0.
      rplresult = unpack(lrpl_ad(1:nw),wavesamples /= 0,field)
      s1result = rsresult * drefl_ds1
      s2result = rsresult * drefl_ds2
      s3result = rsresult * drefl_ds3
      s4result = rsresult * drefl_ds4
      kabresult = drdkab*rlresult + dtdkab*tlresult
      karresult = drdkar*rlresult + dtdkar*tlresult
      kwresult  = drdkw*rlresult + dtdkw*tlresult
      kmresult  = drdm*rlresult + dtdm*tlresult
      knresult  = drdn*rlresult + drdn*tlresult
      do j=1,nbands
        thisbandpass(:) = bandpass(j,:)
	xs1_ad = xs1_ad + sum(s1result * thisbandpass)*totalwavesamples(j)
        xs2_ad = xs2_ad + sum(s2result * thisbandpass)*totalwavesamples(j)
        xs3_ad = xs3_ad + sum(s3result * thisbandpass)*totalwavesamples(j)
        xs4_ad = xs4_ad + sum(s4result * thisbandpass)*totalwavesamples(j)
	xkab_ad = xkab_ad + sum(kabresult*thisbandpass)*totalwavesamples(j)
        xkar_ad = xkar_ad + sum(karresult*thisbandpass)*totalwavesamples(j)
        xkw_ad = xkw_ad + sum(kwresult*thisbandpass)*totalwavesamples(j)
        xkm_ad = xkm_ad + sum(kmresult*thisbandpass)*totalwavesamples(j)
        xleafn_ad = xleafn_ad + sum(knresult*thisbandpass)*totalwavesamples(j)
	xhc_ad = xhc_ad + sum(thisbandpass * hcresult)*totalwavesamples(j)
        xlai_ad= xlai_ad + sum(thisbandpass * lairesult)*totalwavesamples(j)
        rpl_ad = rpl_ad + sum(thisbandpass * rplresult)*totalwavesamples(j)
      enddo
! zero the adjoint of brf
      thisbrf_ad(1,:) = 0.
end

