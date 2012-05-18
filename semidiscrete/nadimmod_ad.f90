      module     multiple_dom_store_ad
      implicit none
      integer :: multkl_multiple_dom
      double precision, allocatable :: multkl_xi_1h(:,:,:)
      double precision, allocatable :: multkl_xi_2h(:,:,:)
      integer :: multl_multiple_dom
      double precision, allocatable :: multl_s_1h(:,:,:,:)
      double precision, allocatable :: multl_s_3h(:,:,:,:)
      double precision, allocatable :: multl_xi_2h(:,:,:,:)
      double precision, allocatable :: multl_xi_4h(:,:,:,:)

      end module     multiple_dom_store_ad

module mo_nad
      integer nw,nwmaxx
      integer, parameter :: nlm = 10   ! down from 50???          ! max number passes in multi-scattering loop
      double precision, parameter :: pi = 3.1415926535898

      double precision,allocatable :: lai_ad(:)
      double precision,allocatable :: rl_ad(:)
      double precision,allocatable  :: tl_ad(:)

      double precision,allocatable :: lai(:)
      double precision,allocatable :: rl(:)
      double precision,allocatable  :: tl(:)

      double precision :: points(32)
      double precision :: weights(32)

      integer :: ild
      integer :: number

      double precision,allocatable  :: xi1_ad(:)
      double precision,allocatable  :: ximt_ad(:)
      double precision,allocatable  :: xi1(:)
      double precision,allocatable  :: ximt(:)

      double precision :: phi_0
      double precision :: teta_0

      double precision,allocatable :: i0(:,:,:)
      double precision,allocatable :: xi1u(:,:,:)
      double precision,allocatable :: xif(:,:,:)

      double precision,allocatable :: i0_ad(:,:,:)
      double precision,allocatable :: xi1u_ad(:,:,:)
      double precision,allocatable :: xif_ad(:,:,:)

      double precision,allocatable :: rs(:)
      double precision,allocatable :: rs_ad(:)

      double precision,allocatable :: c1(:)
      double precision,allocatable  :: c1_ad(:)

      double precision,allocatable :: a_f(:)
      double precision,allocatable :: df(:)
      double precision,allocatable :: h_c(:)
      integer :: n_c
      double precision,allocatable :: r(:)
      double precision,allocatable :: x_ly(:)
      double precision,allocatable :: x_nf(:)

      double precision,allocatable  :: a_f_ad(:)
      double precision,allocatable :: df_ad(:)
      double precision,allocatable :: h_c_ad(:)
      double precision,allocatable :: r_ad(:)


      double precision :: ag
      double precision :: bg
      double precision :: cg
      double precision :: dg
end module mo_nad
      subroutine nad_allocate(nbands_to_use)
      use mo_nad
      use dataSpec_p5
      implicit none
	integer nbands_to_use
! all of these need to be allocated of size nwmax
! but when they are used, it will be just those wavebands we need
! so they all have to be qualified e.g. lai_ad(1:nw)
      call nad_deallocate()
	nwmaxx = nbands_to_use
!      print*,'ALLOC: nad_allocate',nbands_to_use
      allocate(lai_ad(nwmaxx))
      allocate(lai(nwmaxx))
      allocate(rl_ad(nwmaxx))
      allocate(tl_ad(nwmaxx))
      allocate(rs_ad(nwmaxx))
      allocate(rl(nwmaxx))
      allocate(tl(nwmaxx))
      allocate(rs(nwmaxx))
      allocate(xi1_ad(nwmaxx))
      allocate(ximt_ad(nwmaxx))
      allocate(xi1(nwmaxx))
      allocate(ximt(nwmaxx))
      allocate(i0(21,40,nwmaxx))
      allocate(xi1u(21,40,nwmaxx))
      allocate(xif(21,40,nwmaxx))
      allocate(i0_ad(21,40,nwmaxx))
      allocate(xi1u_ad(21,40,nwmaxx))
      allocate(xif_ad(21,40,nwmaxx))
      allocate(c1_ad(nwmaxx))
      allocate(c1(nwmaxx))
      allocate(a_f_ad(nwmaxx))
      allocate(a_f(nwmaxx))
      allocate(df_ad(nwmaxx))
      allocate(df(nwmaxx))
      allocate(h_c_ad(nwmaxx))
      allocate(h_c(nwmaxx))
      allocate(r_ad(nwmaxx))
      allocate(r(nwmaxx))
      allocate(x_nf(nwmaxx))
      allocate(x_ly(nwmaxx))

      end subroutine
      subroutine nad_deallocate()
        use mo_nad
        implicit none
!      print*,'DEALLOC: nad_deallocate'
      if(allocated(lai_ad))deallocate(lai_ad)
      if(allocated(lai))deallocate(lai)
      if(allocated(rl_ad))deallocate(rl_ad)
      if(allocated(tl_ad))deallocate(tl_ad)
      if(allocated(rs_ad))deallocate(rs_ad)
      if(allocated(rl))deallocate(rl)
      if(allocated(tl))deallocate(tl)
      if(allocated(rs))deallocate(rs)
      if(allocated(xi1_ad))deallocate(xi1_ad)
      if(allocated(ximt_ad))deallocate(ximt_ad)
      if(allocated(xi1))deallocate(xi1)
      if(allocated(ximt))deallocate(ximt)
      if(allocated(i0))deallocate(i0)
      if(allocated(xi1u))deallocate(xi1u)
      if(allocated(xif))deallocate(xif)
      if(allocated(i0_ad))deallocate(i0_ad)
      if(allocated(xi1u_ad))deallocate(xi1u_ad)
      if(allocated(xif_ad))deallocate(xif_ad)
      if(allocated(c1_ad))deallocate(c1_ad)
      if(allocated(c1))deallocate(c1)
      if(allocated(a_f_ad))deallocate(a_f_ad)
      if(allocated(a_f))deallocate(a_f)
      if(allocated(df_ad))deallocate(df_ad)
      if(allocated(df))deallocate(df)
      if(allocated(h_c_ad))deallocate(h_c_ad)
      if(allocated(h_c))deallocate(h_c)
      if(allocated(r_ad))deallocate(r_ad)
      if(allocated(r))deallocate(r)
      if(allocated(x_nf))deallocate(x_nf)
      if(allocated(x_ly))deallocate(x_ly)

      end subroutine

      subroutine nad_zero
        use mo_nad
        implicit none
        lai_ad=0
        lai=0
        rl_ad=0
        tl_ad=0
        rs_ad=0
        rl=0
        tl=0
        rs=0
        xi1_ad=0
        ximt_ad=0
        xi1=0
        ximt=0
        i0=0
        xi1u=0
        xif=0
        i0_ad=0
        xi1u_ad=0
        xif_ad=0
        c1_ad=0
        c1=0
        a_f_ad=0
        a_f=0
        df_ad=0
        df=0
        h_c_ad=0
        h_c=0
        r_ad=0
        r=0
        x_nf=0
        x_ly=0
      end subroutine nad_zero
      subroutine nad_mult_n(npt)
      use multiple_dom_store_ad
      use mo_nad
      implicit none
      integer npt
      multkl_multiple_dom = nlm*21*npt
      multl_multiple_dom = nlm*npt
      end subroutine nad_mult_n

      subroutine nad_mult_allocate()
      use multiple_dom_store_ad
      use dataSpec_p5
!	use mo_rtmodel
      use mo_nad
      implicit none
!      print*,'ALLOC: nad_mult_allocate',nw,multkl_multiple_dom,multl_multiple_dom
      allocate( multkl_xi_1h(40,multkl_multiple_dom,nwmaxx))
      allocate( multkl_xi_2h(40,multkl_multiple_dom,nwmaxx))
      allocate(  multl_s_1h(21,40,multl_multiple_dom,nwmaxx))
      allocate( multl_xi_2h(21,40,multl_multiple_dom,nwmaxx))
      allocate(  multl_s_3h(21,40,multl_multiple_dom,nwmaxx))
      allocate(multl_xi_4h(21,40,multl_multiple_dom,nwmaxx))
      end subroutine nad_mult_allocate

      subroutine nad_mult_deallocate()
      use multiple_dom_store_ad
      implicit none
!      print*,'DEALLOC: nad_mult_deallocate'
      if (allocated(multkl_xi_1h))deallocate( multkl_xi_1h )
      if (allocated(multkl_xi_2h))deallocate( multkl_xi_2h )
      if (allocated(multl_s_1h))deallocate( multl_s_1h )
      if (allocated(multl_xi_2h))deallocate( multl_xi_2h )
      if (allocated(multl_s_3h))deallocate( multl_s_3h )
      if (allocated(multl_xi_4h))deallocate( multl_xi_4h )
      end subroutine nad_mult_deallocate

