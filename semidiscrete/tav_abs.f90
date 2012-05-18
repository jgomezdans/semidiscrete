!     ******************************************************************
!     tav(teta0,ref0) computation of the transmittivity at the leaf 
!     surface for a given incidence solid angle. teta0 is the incidence
!     solid angle (in radian). the average angle that works in most 
!     cases is 59deg*pi/180. ref0 is the refaction index.
!
!     jacquemoud s., 1992
! 
!     ******************************************************************
!     stern f., 1964, transmission of isotropic radiation across an
!     interface between two dielectrics, appl.opt., vol.3, 1:111-113
!     allen w.a., 1973, transmission of isotropic light across a
!     dielectric surface in two and three dimensions, j.opt.soc.am.,
!     vol.63, 6:664-666
!	  feret et al. (2008), prospect-4 and 5: advances in the leaf optical
!	  properties model separating photosynthetic pigments, remote sensing of
!	  environment
!     ******************************************************************

subroutine tav_abs(op,teta0,ref0)
      use spectrum_width_p5
      double precision teta,ref,teta0,ref0(nl),op(nl)
      double precision ref2,tav_ab
      double precision a,b,b1,b2,k,pi
      double precision ts,tp,tp1,tp2,tp3,tp4,tp5
      integer count
pi	=	3.1415926

do count=1,nl
teta=   teta0

ref	=	ref0(count)

if (teta.eq.0.) then
	tav_ab=4.*ref/(ref+1.)**2
else

ref2=	ref**2
a	=	(ref+1.)**2/2.
k	=	-(ref2-1.)**2/4.
teta=	pi*teta/180.


if (teta.eq.pi/2.) then
	b1=0.
else
	b1=dsqrt((sin(teta)*sin(teta)-(ref2+1.)/2.)**2+k)
endif

b2	=	sin(teta)*sin(teta)-(ref2+1.)/2.
b	=	b1-b2
ts	=	(k**2/(6.*b**3)+k/b-b/2.)-(k**2/(6.*a**3)+k/a-a/2.)
tp1	=	-2.*ref2*(b-a)/(ref2+1.)**2
tp2	=	-2.*ref2*(ref2+1.)*dlog(b/a)/(ref2-1.)**2
tp3	=	ref2*(1./b-1./a)/2.
tp4	=	16.*ref2**2*(ref2**2+1.)*dlog((2.*(ref2+1.)*b-(ref2-1.)**2)/ &
(2.*(ref2+1.)*a-(ref2-1.)**2))/((ref2+1.)**3*(ref2-1.)**2)
tp5	=	16.*ref2**3*(1./(2.*(ref2+1.)*b-(ref2-1.)**2)-1./(2.*(ref2&
+1.)*a-(ref2-1.)**2))/(ref2+1.)**3
tp	=	tp1+tp2+tp3+tp4+tp5
tav_ab	=	(ts+tp)/(2.*sin(teta)*sin(teta))
endif
op(count) = tav_ab
enddo
return
end
