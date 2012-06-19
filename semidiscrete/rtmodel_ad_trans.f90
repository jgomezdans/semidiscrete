! -*- f90 -*-
      module mo_rtmodel
      implicit none
      integer nwmax
      double precision, parameter :: rtmodel_spectralInterval = 1.0
      double precision, parameter :: rtmodel_spectralMin = 400.0
      double precision, parameter :: rtmodel_spectralMax = 2500.0
      integer, parameter :: rtmodel_nparam = 13 
      integer, parameter :: rtmodel_version = 1
      integer, parameter :: rtmodel_subversion = 1
      integer :: lad
      double precision :: xlai
      double precision :: xhc
      double precision :: rpl
      double precision :: xs1,xs2,xs3,xs4
      double precision :: xkab,xkar,xkw,xkm,xleafn
      double precision :: xlai_ad
      double precision :: xhc_ad
      double precision :: rpl_ad
      double precision :: xs1_ad,xs2_ad,xs3_ad,xs4_ad
      double precision :: xkab_ad,xkar_ad,xkw_ad,xkm_ad,xleafn_ad
      double precision, parameter :: slai=-2.0
      double precision, parameter :: skab=-100.0
      double precision, parameter :: skar=-100.0
      double precision, parameter :: skw=(-1./50.)
      double precision, parameter :: skm=(-1./100.)
      character(2000),parameter :: rtmodel_param_names = "exp(-LAI/2.0):Canopy height:Leaf radius:&
           exp(-chlophyll/100.):leaf brown pigments:exp(-50. leaf water):exp(-100 leaf dry matter):&
           leaf layers (N):soil cmpt 1:&
           soil cmpt 2:soil cmpt 3:soil cmpt 4:lidf:"
      character(2000),parameter :: rtmodel_param_minmax = &
   "0.001;0.999:0.001;10:0.0001;1.0:0.001;0.999:0;1:0.001;0.999:0.001;0.999:0.5;2.5:0;1:0;1:0;1:0;1:0;5:"

! NB string list termined and separated by :
      contains 
      subroutine mo_rtzero_ad
      rpl_ad = 0.
      xhc_ad = 0.
      xlai_ad = 0.
      xs1_ad = 0.
      xs2_ad = 0.
      xs3_ad = 0.
      xs4_ad = 0.
      xkab_ad = 0.
      xkar_ad = 0.
      xkw_ad = 0.
      xkm_ad = 0.
      xleafn_ad = 0.
      end subroutine mo_rtzero_ad
      subroutine mo_rtsum_ad(x)
      implicit none
      double precision :: x(:)
      if(xlai_ad .ne. 0 .and. xlai .gt. 0 )x(1) = x(1) + xlai_ad * (slai/exp(xlai/slai))
      x(2) = x(2) + xhc_ad
      x(3) = x(3) + rpl_ad
      if(xkab_ad .ne. 0 .and. xkab .gt. 0) x(4) = x(4) + xkab_ad * (skab/exp(xkab/skab))
      if(xkar_ad .ne. 0 .and. xkar .gt. 0) x(5) = x(5) +  xkar_ad !* (skar/exp(xkar/skar))
      if(xkw_ad .ne. 0 .and. xkw .gt. 0 )x(6) = x(6) + xkw_ad * ((skw)/exp(xkw/skw))
      if(xkm_ad .ne. 0 .and. xkm .gt. 0)x(7) = x(7) + xkm_ad * ((skm)/exp(xkm/skm))
      x(8) = x(8) + xleafn_ad
      x(9) = x(9) + xs1_ad
      x(10) = x(10) + xs2_ad
      x(11) = x(11) + xs3_ad
      x(12) = x(12) + xs4_ad
! zero for lad
      x(13) = 0
      call mo_rtzero_ad()
      end subroutine mo_rtsum_ad

      subroutine mo_rtload(x)
      implicit none
      double precision :: x(:)
      xlai=0
      if(x(1) .gt. 0)xlai=slai*log(x(1))
      xhc=x(2)
      rpl=x(3)
      xkab=0
      if(x(4).gt.0)xkab=skab*log(x(4))
      xkar=0
      if(x(5) .gt.0)xkar=x(5) !skar*log(x(5))
      xkw=0.
      if(x(6) .gt. 0)xkw =skw*log(x(6))
      xkm =0
      if(x(7) .gt. 0 )xkm =skm*log(x(7))
      xleafn =x(8)
      xs1 =x(9)
      xs2 =x(10)
      xs3 =x(11)
      xs4 =x(12)
      lad = int(x(13)+0.5)
      end subroutine mo_rtload
! interface fn that doesnt require the module to be loaded
      end module mo_rtmodel

! this gets called before the rt model is ever called
! use it for any allocation / initialisation prep
subroutine rt_modelpre(nbands_to_use,bands_to_use)
	use mo_rtmodel
      use mo_nad
      implicit none
      integer nbands_to_use
      integer bands_to_use(nbands_to_use)
      nwmax = nbands_to_use
      nwmaxx = nbands_to_use
      call leafprep(nbands_to_use,bands_to_use)
      call nad_allocate(nbands_to_use)
end subroutine rt_modelpre

! use this eg for priniting anything interesting 
! or any additional deallocation etc
subroutine rt_modelpost()
      return
end subroutine rt_modelpost

subroutine rt_modeldpre(npt)
      implicit none
      integer npt,nw
! set the dimensions of multkl_multiple_dom
! and multl_multiple_dom
      call nad_mult_n(npt)
! make sure deallocated first
      call nad_mult_deallocate()
      call nad_mult_allocate()
end subroutine rt_modeldpre

subroutine rt_modeldpost()
	use dataspec_p5
      implicit none

      call deallocatespectra()
      call nad_mult_deallocate()
end subroutine rt_modeldpost

subroutine rt_model(nv,nbands,nl,ipt,x,theta_v,phi_v,theta_i,phi_i,brf,wb)
      use mo_rtmodel
!!!!!f2py intent(c) theta_v,phi_v,x,wb
      implicit none
      integer :: nv,nbands,nl
      integer :: ipt
      double precision, intent(out) :: brf(nv,nbands)
      double precision :: theta_v(nv),phi_v(nv),theta_i,phi_i
      double precision :: x(rtmodel_nparam)
      double precision :: wb(nbands,nl)

      call mo_rtload(x)
      brf = 0.
	!print*,size(wb(:,1)),nl,nv,nwmax
      call msnadimbrf(size(wb(:,1)),nl,nv,lad,xkab,xkar,xkw,xkm,xleafn,&
           xs1,xs2,xs3,xs4,theta_i,phi_i,theta_v,phi_v,xhc,xlai,rpl,brf,wb)
      return
end subroutine rt_model

subroutine rt_modelpred(nv,nbands,nl,ipt,x,theta_v,phi_v,theta_i,phi_i,brf,wb)
      use mo_rtmodel
!!!!!f2py intent(c) theta_v,phi_v,x,wb
      implicit none
      integer :: nv,nbands,nl
      integer :: ipt
      double precision, intent(out) :: brf(nv,nbands)
      double precision :: theta_v(nv),phi_v(nv),theta_i,phi_i
      double precision :: x(rtmodel_nparam)
      double precision :: wb(nbands,nl)

      call mo_rtload(x)
      brf = 0.
      call msnadimbrfmd(ipt,size(wb(:,1)),nl,nv,lad,xkab,xkar,xkw,xkm,xleafn,&
           xs1,xs2,xs3,xs4,theta_i,phi_i,theta_v,phi_v,xhc,xlai,rpl,brf,wb)
      return
end subroutine rt_modelpred

subroutine rt_modeld(nv,nbands,nl,ipt,x,theta_v,phi_v,theta_i,phi_i,brff_ad,wb,x_ad)
      use mo_rtmodel
      use mo_nad
!!!!!f2py intent(c) theta_v,phi_v,x,wb,x_ad
      implicit none
      integer :: nv,nbands,nl
      integer :: ipt
      double precision,intent(in) :: brff_ad(nbands)
      double precision,intent(out) :: x_ad(rtmodel_nparam)
      double precision :: x(rtmodel_nparam)
      double precision :: theta_v(nv),phi_v(nv),theta_i,phi_i
      double precision :: wb(nbands,nl)
      double precision :: brf_ad(nv,nbands)

      brf_ad(1,:) = brff_ad
      call mo_rtload(x)
      call mo_rtzero_ad()
      call nad_zero()
      call msnadimbrf_ad(ipt,size(wb(:,1)),nl,nv,lad,xkab,xkar,xkw,xkm,xleafn,xs1,xs2,xs3,xs4,theta_i,phi_i,theta_v,phi_v,&
           xhc,xlai,rpl,brf_ad,xlai_ad,xhc_ad,rpl_ad,xkab_ad,xkar_ad,xkw_ad,xkm_ad,xleafn_ad,&
           xs1_ad,xs2_ad,xs3_ad,xs4_ad,wb)
      call mo_rtsum_ad(x_ad)

end subroutine rt_modeld

subroutine rt_getnparams(n)
      use mo_rtmodel
      implicit none
      integer, intent(out) :: n
      n = rtmodel_nparam
end
subroutine rt_getspectralinterval(n)
      use mo_rtmodel
      implicit none
      double precision, intent(out) :: n
      n = rtmodel_spectralInterval
end
subroutine rt_getspectralmin(n)
      use mo_rtmodel
      implicit none
      double precision, intent(out) :: n
      n = rtmodel_spectralMin
end
subroutine rt_getspectralmax(n)
      use mo_rtmodel
      implicit none
      double precision, intent(out) :: n
      n = rtmodel_spectralMax
end

!subroutine rt_setnlmax(n)
!!      use mo_rtmodel
!      use mo_nad
!      implicit none
!      double precision, intent(in) :: n
!      nwmax = n 
!      nwmaxx = n
!end
subroutine rt_getnames(names)
      use mo_rtmodel
      use mo_nad
      implicit none
      character(2000),intent(out) :: names
      names = rtmodel_param_names
end

subroutine rt_getminmax(names)
      use mo_rtmodel
      use mo_nad
      implicit none
      character(2000),intent(out) :: names
      names = rtmodel_param_minmax
end



