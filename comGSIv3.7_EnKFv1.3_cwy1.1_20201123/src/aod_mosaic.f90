module aod_mosaic

use gridmod            , only: istart, jstart, nlon, nlat, lon1, nsig
use constants          , only: zero, one, r1000, grav
use radiance_mod       , only: n_actual_aerosols, aerosol_names, n_aerosols_jac, aerosol_names_jac
use aeroinfo           , only: nsigaerojac, aerojacindxs
use kinds              , only: r_kind,r_single,i_kind
use obsmod             , only: rmin_value, sf_jac, nsize_mosaic_wrf, n_cvmas, n_cvnum, n_cvaer
use jfunc              , only: jiter, iter

implicit none

public  call_aod_mosaic

logical                                       ::  ini_fit=.true.
!!! cwy mosaic ------------------------------------------------------------
integer(i_kind),parameter                     ::  nswbands=5
!!! cwy mosaic ------------------------------------------------------------
integer(i_kind), parameter                    ::  prefr=7, prefi=7, nsiz=200, ncoef=50
real(r_kind), parameter                       ::  pie=4.*atan(1.)
integer                                       ::  nrefr, nrefi
real, dimension(ncoef,prefr,prefi,nswbands)   ::  extpsw, scapsw, abspsw, asmpsw
real, dimension(prefr,nswbands)               ::  refrtabsw
real, dimension(prefi,nswbands)               ::  refitabsw

contains
!=========================================================================
subroutine call_aod_mosaic( obstime, data_s, nchanl, nreal, layer_od, zhgt, jacobian_aero, typename )

use guess_grids       , only: ges_prsl, ges_tsen, hrdifsig, nfldsig, geop_hgti
use gsi_bundlemod     , only: gsi_bundlegetpointer
use gsi_metguess_mod  , only: gsi_metguess_bundle
use gsi_chemguess_mod , only: gsi_chemguess_bundle
use control_vectors   , only: cvars => nrf_var      ! all var names in control_vectors:: in anavinfo
use berror            , only: dssv
use mpimod            , only: mype
use mpeu_util         , only: getindex

implicit none

real(r_kind),intent(in)                                ::  obstime
real(r_kind),dimension(nchanl+nreal),intent(in)        ::  data_s
integer(i_kind),intent(in)                             ::  nchanl, nreal
real(r_kind),dimension(nsig,nchanl),intent(out)        ::  layer_od
real(r_kind),dimension(nsig),intent(out)               ::  zhgt
real(r_kind),dimension(nsigaerojac,nchanl),intent(out) ::  jacobian_aero
character( len=* ), intent(in)                         ::  typename

real(r_kind),dimension(nsig)                           ::  ugkg_gcm3, nmkg_ncm3 !!!, zhgt

integer(i_kind)                                        ::  indx, nn, knc
integer(i_kind)                                        ::  ix1, iy1, ix, iy, ixp, iyp, m1, ier, istatus
integer(i_kind)                                        ::  j0, kk, itsig, itsigp, ksw, ksz, kaer, kch, kctr
real(r_kind)                                           ::  w00, w01, w10, w11, ww, dsw1, dsw0
real(r_kind)                                           ::  dx, dy, delx, dely, delx1, dely1
real(r_kind)                                           ::  dtsig, dtsigp, totaod_channel, binaod_channel
real(r_kind),dimension(nsig)                           ::  prsl, temp, shum, den_air
real(r_kind),dimension(nsig+1)                         ::  hegt

real(r_kind),pointer,dimension(:,:,:)                  ::  qges_itsig=>NULL()
real(r_kind),pointer,dimension(:,:,:)                  ::  qges_itsigp=>NULL()
real(r_kind),pointer,dimension(:,:,:)                  ::  aeroges_itsig=>NULL()

real(r_kind)                        ::  dlo_sect, dhi_sect, num_aer_lo, num_aer_hi
real(r_kind)                        ::  dp_wet, vol_dry, vol_wet
complex                             ::  ri_ave
real(r_kind), dimension(nswbands)   ::  wavmidsw

real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  num_aer, rad_wet
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  vol_so4
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  vol_no3
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  vol_nh4
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  vol_ocx
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  vol_bcx
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  vol_clx
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  vol_nax
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  vol_oin
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  vol_h2o

real(r_kind)                        ::  den_so4
real(r_kind)                        ::  den_no3
real(r_kind)                        ::  den_nh4
real(r_kind)                        ::  den_ocx
real(r_kind)                        ::  den_bcx
real(r_kind)                        ::  den_clx
real(r_kind)                        ::  den_nax
real(r_kind)                        ::  den_oin
real(r_kind)                        ::  den_h2o

real(r_kind), dimension(nswbands)   ::  refrsw_so4, refisw_so4
real(r_kind), dimension(nswbands)   ::  refrsw_no3, refisw_no3
real(r_kind), dimension(nswbands)   ::  refrsw_nh4, refisw_nh4
real(r_kind), dimension(nswbands)   ::  refrsw_ocx, refisw_ocx
real(r_kind), dimension(nswbands)   ::  refrsw_bcx, refisw_bcx
real(r_kind), dimension(nswbands)   ::  refrsw_clx, refisw_clx
real(r_kind), dimension(nswbands)   ::  refrsw_nax, refisw_nax
real(r_kind), dimension(nswbands)   ::  refrsw_oin, refisw_oin
real(r_kind), dimension(nswbands)   ::  refrsw_h2o, refisw_h2o

complex, dimension(nswbands)        ::  swref_index_so4
complex, dimension(nswbands)        ::  swref_index_no3
complex, dimension(nswbands)        ::  swref_index_nh4
complex, dimension(nswbands)        ::  swref_index_ocx
complex, dimension(nswbands)        ::  swref_index_bcx
complex, dimension(nswbands)        ::  swref_index_clx
complex, dimension(nswbands)        ::  swref_index_nax
complex, dimension(nswbands)        ::  swref_index_oin
complex, dimension(nswbands)        ::  swref_index_h2o

complex, dimension(nsig,nsize_mosaic_wrf,nswbands) ::  swrefindx_wet, swrefindx_dry
real(r_kind),dimension(nsig,nsize_mosaic_wrf,nswbands)      ::  weighte, weights

integer(i_kind)                     ::  knr, kni, itab, jtab
real(r_kind)                        ::  rmin, rmax, refrmin, refrmax, refimin, refimax
real(r_kind)                        ::  drefr, drefi, bma, bpa, xr
real(r_kind)                        ::  xrmin, xrmax, thesum, xlog, xrad
real(r_kind)                        ::  pext, psca, pasm
complex*16                          ::  crefd, refc
real*8                              ::  thesize
real*8, dimension(nsiz)             ::  qext, qsca, gqsc
real, dimension(nsiz)               ::  qext4, qsca4, qabs4, asymm
real, dimension(nsiz)               ::  rs
real, dimension(ncoef)              ::  cext, csca, casm, ch
real, dimension(nsig,nswbands)      ::  swtauaer, swextaer, swssaaer, swgymaer

real(r_kind),dimension(nsig,n_actual_aerosols)  ::  aero

! miev0
real                                ::  refr, refi, ttab, utab
logical                             ::  prefct=.false., anyang=.false. 
real*8, dimension(0:7,1)            ::  pmom
complex*16                          ::  sforw, sback, tforw(2), tback(2)
complex*16, dimension(1)            ::  s1, s2
real*8                              ::  mimcut=0.
real*8, dimension(1)                ::  xmu=(/1./)
integer                             ::  numang=0, nmom=7, ipolzen=0, momdim=7
logical, dimension(2)               ::  prnt=(/.false.,.false./)

!!! total_aod ----------------------------------------
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  rad_dry
real(r_kind),dimension(nsig,nsize_mosaic_wrf)      ::  vol_lay_dryaer, vol_lay_wetaer
real, dimension(nsig,nsize_mosaic_wrf,nswbands)    ::  betaR, betaI
real, dimension(nsig,nsize_mosaic_wrf,nswbands)    ::  eext_lay

real(r_kind)                        ::  pabs, dp_dry
real                                ::  Iup_Ilow, Rup_Rlow, uu, vv
real                                ::  a00, a01, a10, a11
real                                ::  rden
real                                ::  beta2, diff_aero, alpha, aa, bb

real, dimension(nsig,nsize_mosaic_wrf)         ::  rden_mix
!!! total_aod ----------------------------------------
real                                ::  armin, armax, axx

!----------------------------
integer                             ::  ilat, ilon
! different from the set in chemmod.f90
ilat=4
ilon=3

! COMMENTS IN WRF-CHEM chem/module_optical_averaging.F
!
! Define refractive indicies
! * assume na and cl are the same as nacl
! * assume so4, no3, and nh4 are the same as nh4no3
! * assume ca and co3 are the same as caco3
! * assume msa is just msa
! Further work:
! * to be more precise, need to compute electrolytes to apportion
!   so4, no3, nh4, na, cl, msa, ca, co3 among various componds
!   as was done previously in module_mosaic_therm.F

!----------------------------
! 300, 400, 600, 999 nm in WRF-Chem
!!! cwy mosaic ------------------------------------------------------------
!              440nm   , 550nm   , 675nm   , 870nm   , 1020nm
wavmidsw   =(/ 0.44e-4 , 0.55e-4 , 0.675e-4, 0.87e-4 , 1.02e-4 /) ! cm

refrsw_so4 =(/ 1.540   , 1.530   , 1.525   , 1.520   , 1.520   /)
refisw_so4 =(/ 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  /)

refrsw_no3 =(/ 1.540   , 1.530   , 1.525   , 1.520   , 1.520   /)
refisw_no3 =(/ 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  /)

refrsw_nh4 =(/ 1.540   , 1.530   , 1.525   , 1.520   , 1.520   /)
refisw_nh4 =(/ 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  /)

refrsw_ocx =(/ 1.550   , 1.550   , 1.550   , 1.550   , 1.550   /)
refisw_ocx =(/ 0.001   , 0.001   , 0.001   , 0.001   , 0.001   /)

refrsw_bcx =(/ 1.950   , 1.950   , 1.950   , 1.950   , 1.950   /)
refisw_bcx =(/ 0.790   , 0.790   , 0.790   , 0.790   , 0.790   /)

refrsw_clx =(/ 1.560   , 1.550   , 1.540   , 1.530   , 1.530   /)
refisw_clx =(/ 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  /)

refrsw_nax =(/ 1.560   , 1.550   , 1.540   , 1.530   , 1.530   /)
refisw_nax =(/ 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  /)

refrsw_oin =(/ 1.530   , 1.530   , 1.530   , 1.530   , 1.530   /) ! dust
refisw_oin =(/ 0.003   , 0.002   , 0.0015  , 0.0010  , 0.0010  /) ! dust

refrsw_h2o =(/ 1.337   , 1.335   , 1.332   , 1.330   , 1.328   /)
refisw_h2o =(/ 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  , 1.0e-7  /)

! g/cm3
den_so4=1.8
den_no3=1.8
den_nh4=1.8
den_ocx=1.4
den_bcx=1.8
den_clx=2.2
den_nax=2.2
den_oin=2.6
den_h2o=1.0
!!! cwy mosaic ------------------------------------------------------------

!----------------------------
do ksw=1, nswbands
  swref_index_so4(ksw)=cmplx( refrsw_so4(ksw), refisw_so4(ksw) )
  swref_index_no3(ksw)=cmplx( refrsw_no3(ksw), refisw_no3(ksw) )
  swref_index_nh4(ksw)=cmplx( refrsw_nh4(ksw), refisw_nh4(ksw) )
  swref_index_ocx(ksw)=cmplx( refrsw_ocx(ksw), refisw_ocx(ksw) )
  swref_index_bcx(ksw)=cmplx( refrsw_bcx(ksw), refisw_bcx(ksw) )
  swref_index_clx(ksw)=cmplx( refrsw_clx(ksw), refisw_clx(ksw) )
  swref_index_nax(ksw)=cmplx( refrsw_nax(ksw), refisw_nax(ksw) )
  swref_index_oin(ksw)=cmplx( refrsw_oin(ksw), refisw_oin(ksw) )
  swref_index_h2o(ksw)=cmplx( refrsw_h2o(ksw), refisw_h2o(ksw) )
end do

!----------------------------
rmin=0.005e-4
rmax=50.e-4
!=========================================== ini_fit start
if( ini_fit ) then
  ini_fit=.false.

  refrmin= 1.50; refrmax=refrmin
  refimin=-1.50; refimax=refimin
  do ksw=1, nswbands
    do kaer=1, 9
      if( kaer==1 ) then; refr=refrsw_so4(ksw); refi=-1*refisw_so4(ksw); endif
      if( kaer==2 ) then; refr=refrsw_no3(ksw); refi=-1*refisw_no3(ksw); endif
      if( kaer==3 ) then; refr=refrsw_nh4(ksw); refi=-1*refisw_nh4(ksw); endif
      if( kaer==4 ) then; refr=refrsw_ocx(ksw); refi=-1*refisw_ocx(ksw); endif
      if( kaer==5 ) then; refr=refrsw_bcx(ksw); refi=-1*refisw_bcx(ksw); endif
      if( kaer==6 ) then; refr=refrsw_clx(ksw); refi=-1*refisw_clx(ksw); endif
      if( kaer==7 ) then; refr=refrsw_nax(ksw); refi=-1*refisw_nax(ksw); endif
      if( kaer==8 ) then; refr=refrsw_oin(ksw); refi=-1*refisw_oin(ksw); endif
      if( kaer==9 ) then; refr=refrsw_h2o(ksw); refi=-1*refisw_h2o(ksw); endif
      refrmin=min(refrmin,refr); refrmax=max(refrmax,refr)
      refimin=min(refimin,refi); refimax=max(refimax,refi)
    end do

    drefr=refrmax-refrmin
    if( drefr>1.e-4 ) then
      nrefr=prefr
      drefr=drefr/(nrefr-1)
    else
      nrefr=1
    end if

    drefi=refimax-refimin
    if( drefi>1.e-4 ) then
      nrefi=prefi
      drefi=drefi/(nrefi-1)
    else
      nrefi=1
    end if

    bma=0.5*log(rmax/rmin)
    bpa=0.5*log(rmax*rmin)

    do knr=1, nrefr
    do kni=1, nrefi
      refrtabsw(knr,ksw)=refrmin+(knr-1)*drefr
      refitabsw(kni,ksw)=refimin/0.2*(0.2**real(kni))
      if( kni==nrefi ) refitabsw(kni,ksw)=-1.0e-20
      crefd=cmplx( refrtabsw(knr,ksw),refitabsw(kni,ksw) )

      do ksz=1, nsiz
        xr=cos(pie*(float(ksz)-0.5)/float(nsiz))
        rs(ksz)=exp(xr*bma+bpa)
        thesize=2.*pie*rs(ksz)/wavmidsw(ksw)      ! wavmidsw
        thesize=min(thesize,10000.d0)

        call miev0( thesize,crefd,prefct,mimcut,anyang,numang,xmu,nmom,ipolzen,momdim,prnt,&
                    qext(ksz),qsca(ksz),gqsc(ksz),pmom,sforw,sback,s1,s2,tforw,tback ) 
        
        ! qext: mass extinction coef, L^2 M^-1
        qext4(ksz)=qext(ksz)
        qsca4(ksz)=min(qsca(ksz),qext(ksz))
        qabs4(ksz)=qext4(ksz)-qsca4(ksz)
        asymm(ksz)=gqsc(ksz)/qsca4(ksz)

        qext4(ksz)=max(qext4(ksz),rmin_value)
        qsca4(ksz)=max(qsca4(ksz),rmin_value)
        qabs4(ksz)=max(qabs4(ksz),rmin_value)
        asymm(ksz)=max(asymm(ksz),rmin_value)
      end do

      call fitcurv( rs,qext4,extpsw(1,knr,kni,ksw),ncoef,nsiz )
      call fitcurv( rs,qsca4,scapsw(1,knr,kni,ksw),ncoef,nsiz )
      call fitcurv( rs,qabs4,abspsw(1,knr,kni,ksw),ncoef,nsiz )
      call fitcurv( rs,asymm,asmpsw(1,knr,kni,ksw),ncoef,nsiz )
    end do
    end do
  end do
end if
!=========================================== ini_fit end

!----------------------------
m1=mype+1

dx=data_s(ilat)
dy=data_s(ilon)

if( dx<0 .or. dy<0 ) then
  layer_od=zero
  jacobian_aero=zero
  return
end if

ix1=dx
ix1=max(1,min(ix1,nlat))
delx=dx-ix1
delx=max(zero,min(delx,one))
ix=ix1-istart(m1)+2
ixp=ix+1
if(ix1==nlat) then
  ixp=ix
end if
delx1=one-delx

iy1=dy
dely=dy-iy1
iy=iy1-jstart(m1)+2
if(iy<1) then
  iy1=iy1+nlon
  iy=iy1-jstart(m1)+2
end if
if(iy>lon1+1) then
  iy1=iy1-nlon
  iy=iy1-jstart(m1)+2
end if
iyp=iy+1
dely1=one-dely

w00=delx1*dely1; w10=delx*dely1; w01=delx1*dely; w11=delx*dely

if( ix<0 .or. iy<0 ) then
print *, 'cwy warning: mype, ns-ix, we-iy, dx, dy, istart, jstart=', mype, ix, iy, dx, dy, istart(m1), jstart(m1)
stop
end if

if( obstime>hrdifsig(1) .and. obstime<hrdifsig(nfldsig) ) then
  do j0=1, nfldsig-1
    if( obstime>hrdifsig(j0) .and. obstime<=hrdifsig(j0+1) ) then
      itsig=j0
      itsigp=j0+1
      dtsig=((hrdifsig(j0+1)-obstime)/(hrdifsig(j0+1)-hrdifsig(j0)))
    end if
  end do
else if( obstime<=hrdifsig(1) ) then
  itsig=1
  itsigp=1
  dtsig=one
else
  itsig=nfldsig
  itsigp=nfldsig
  dtsig=one
end if
dtsigp=one-dtsig

!----------------------------
call gsi_bundlegetpointer(gsi_metguess_bundle(itsig ),'q',qges_itsig ,istatus)
call gsi_bundlegetpointer(gsi_metguess_bundle(itsigp),'q',qges_itsigp,istatus)

do kk=1, nsig
  prsl(kk)=( ges_prsl(ix ,iy ,kk,itsig )*w00+ &         ! 0.1*hPa
             ges_prsl(ixp,iy ,kk,itsig )*w10+ &
             ges_prsl(ix ,iyp,kk,itsig )*w01+ &
             ges_prsl(ixp,iyp,kk,itsig )*w11 )*dtsig+ &
           ( ges_prsl(ix ,iy ,kk,itsigp)*w00+ &
             ges_prsl(ixp,iy ,kk,itsigp)*w10+ &
             ges_prsl(ix ,iyp,kk,itsigp)*w01+ &
             ges_prsl(ixp,iyp,kk,itsigp)*w11 )*dtsigp
  prsl(kk)=prsl(kk)*1000.                               ! Pa

  temp(kk)=( ges_tsen(ix ,iy ,kk,itsig )*w00+ &         ! K
             ges_tsen(ixp,iy ,kk,itsig )*w10+ &
             ges_tsen(ix ,iyp,kk,itsig )*w01+ &
             ges_tsen(ixp,iyp,kk,itsig )*w11 )*dtsig+ &
           ( ges_tsen(ix ,iy ,kk,itsigp)*w00+ &
             ges_tsen(ixp,iy ,kk,itsigp)*w10+ &
             ges_tsen(ix ,iyp,kk,itsigp)*w01+ &
             ges_tsen(ixp,iyp,kk,itsigp)*w11 )*dtsigp

  shum(kk)=( qges_itsig (ix ,iy ,kk)*w00+ &             ! kg/kg
             qges_itsig (ixp,iy ,kk)*w10+ &
             qges_itsig (ix ,iyp,kk)*w01+ &
             qges_itsig (ixp,iyp,kk)*w11 )*dtsig+ &
           ( qges_itsigp(ix ,iy ,kk)*w00+ &
             qges_itsigp(ixp,iy ,kk)*w10+ &
             qges_itsigp(ix ,iyp,kk)*w01+ &
             qges_itsigp(ixp,iyp,kk)*w11 )*dtsigp
end do
do kk=1, nsig+1
  hegt(kk)=( geop_hgti(ix ,iy ,kk,itsig )*w00+ &         ! m
             geop_hgti(ixp,iy ,kk,itsig )*w10+ &
             geop_hgti(ix ,iyp,kk,itsig )*w01+ &
             geop_hgti(ixp,iyp,kk,itsig )*w11 )*dtsig+ &
           ( geop_hgti(ix ,iy ,kk,itsigp)*w00+ &
             geop_hgti(ixp,iy ,kk,itsigp)*w10+ &
             geop_hgti(ix ,iyp,kk,itsigp)*w01+ &
             geop_hgti(ixp,iyp,kk,itsigp)*w11 )*dtsigp
end do
do kk=1, nsig
  zhgt(kk)=hegt(kk+1)-hegt(kk)
end do

do kk=1, nsig
  den_air(kk)=prsl(kk)/((1+0.608*shum(kk))*temp(kk)*287.058)      ! kg/m^3
  ugkg_gcm3(kk)=den_air(kk)*1.e-12                                ! ug/kg -> g/cm^3
  nmkg_ncm3(kk)=den_air(kk)*1.e-6                                 !  #/kg -> #/cm^3
end do

!----------------------------
do kaer=1, n_actual_aerosols
  call gsi_bundlegetpointer(gsi_chemguess_bundle(1),aerosol_names(kaer),aeroges_itsig,ier)

  do kk=1, nsig
    if( index(aerosol_names(kaer),'num_a')>0 ) then
      aero(kk,kaer)=( aeroges_itsig(ix ,iy ,kk)*w00+ &         ! num  in  #/kg -> #/cm^3
                      aeroges_itsig(ixp,iy ,kk)*w10+ &
                      aeroges_itsig(ix ,iyp,kk)*w01+ &
                      aeroges_itsig(ixp,iyp,kk)*w11 )*nmkg_ncm3(kk)
    else if( index(aerosol_names(kaer),'mas_a')==0 ) then
      aero(kk,kaer)=( aeroges_itsig(ix ,iy ,kk)*w00+ &         ! conc in ug/kg -> g/cm^3
                      aeroges_itsig(ixp,iy ,kk)*w10+ &
                      aeroges_itsig(ix ,iyp,kk)*w01+ &
                      aeroges_itsig(ixp,iyp,kk)*w11 )*ugkg_gcm3(kk)
    else
      cycle
    end if
    if( aero(kk,kaer)<0 ) then
      print *, 'cwy warning aero<0 --- aero=', aero(kk,kaer), trim(aerosol_names(kaer))
      stop
    end if
  end do

  select case( trim(aerosol_names(kaer)) )
  case( 'so4_a01' ); vol_so4(:,1)=aero(:,kaer)/den_so4 ! g/cm3 / g/cm3
  case( 'so4_a02' ); vol_so4(:,2)=aero(:,kaer)/den_so4
  case( 'so4_a03' ); vol_so4(:,3)=aero(:,kaer)/den_so4
  case( 'so4_a04' ); vol_so4(:,4)=aero(:,kaer)/den_so4
  case( 'so4_a05' ); vol_so4(:,5)=aero(:,kaer)/den_so4
  case( 'so4_a06' ); vol_so4(:,6)=aero(:,kaer)/den_so4
  case( 'so4_a07' ); vol_so4(:,7)=aero(:,kaer)/den_so4
  case( 'so4_a08' ); vol_so4(:,8)=aero(:,kaer)/den_so4

  case( 'no3_a01' ); vol_no3(:,1)=aero(:,kaer)/den_no3
  case( 'no3_a02' ); vol_no3(:,2)=aero(:,kaer)/den_no3
  case( 'no3_a03' ); vol_no3(:,3)=aero(:,kaer)/den_no3
  case( 'no3_a04' ); vol_no3(:,4)=aero(:,kaer)/den_no3
  case( 'no3_a05' ); vol_no3(:,5)=aero(:,kaer)/den_no3
  case( 'no3_a06' ); vol_no3(:,6)=aero(:,kaer)/den_no3
  case( 'no3_a07' ); vol_no3(:,7)=aero(:,kaer)/den_no3
  case( 'no3_a08' ); vol_no3(:,8)=aero(:,kaer)/den_no3

  case( 'nh4_a01' ); vol_nh4(:,1)=aero(:,kaer)/den_nh4
  case( 'nh4_a02' ); vol_nh4(:,2)=aero(:,kaer)/den_nh4
  case( 'nh4_a03' ); vol_nh4(:,3)=aero(:,kaer)/den_nh4
  case( 'nh4_a04' ); vol_nh4(:,4)=aero(:,kaer)/den_nh4
  case( 'nh4_a05' ); vol_nh4(:,5)=aero(:,kaer)/den_nh4
  case( 'nh4_a06' ); vol_nh4(:,6)=aero(:,kaer)/den_nh4
  case( 'nh4_a07' ); vol_nh4(:,7)=aero(:,kaer)/den_nh4
  case( 'nh4_a08' ); vol_nh4(:,8)=aero(:,kaer)/den_nh4

  case( 'oc_a01'  ); vol_ocx(:,1)=aero(:,kaer)/den_ocx
  case( 'oc_a02'  ); vol_ocx(:,2)=aero(:,kaer)/den_ocx
  case( 'oc_a03'  ); vol_ocx(:,3)=aero(:,kaer)/den_ocx
  case( 'oc_a04'  ); vol_ocx(:,4)=aero(:,kaer)/den_ocx
  case( 'oc_a05'  ); vol_ocx(:,5)=aero(:,kaer)/den_ocx
  case( 'oc_a06'  ); vol_ocx(:,6)=aero(:,kaer)/den_ocx
  case( 'oc_a07'  ); vol_ocx(:,7)=aero(:,kaer)/den_ocx
  case( 'oc_a08'  ); vol_ocx(:,8)=aero(:,kaer)/den_ocx

  case( 'bc_a01'  ); vol_bcx(:,1)=aero(:,kaer)/den_bcx
  case( 'bc_a02'  ); vol_bcx(:,2)=aero(:,kaer)/den_bcx
  case( 'bc_a03'  ); vol_bcx(:,3)=aero(:,kaer)/den_bcx
  case( 'bc_a04'  ); vol_bcx(:,4)=aero(:,kaer)/den_bcx
  case( 'bc_a05'  ); vol_bcx(:,5)=aero(:,kaer)/den_bcx
  case( 'bc_a06'  ); vol_bcx(:,6)=aero(:,kaer)/den_bcx
  case( 'bc_a07'  ); vol_bcx(:,7)=aero(:,kaer)/den_bcx
  case( 'bc_a08'  ); vol_bcx(:,8)=aero(:,kaer)/den_bcx

  case( 'cl_a01'  ); vol_clx(:,1)=aero(:,kaer)/den_clx
  case( 'cl_a02'  ); vol_clx(:,2)=aero(:,kaer)/den_clx
  case( 'cl_a03'  ); vol_clx(:,3)=aero(:,kaer)/den_clx
  case( 'cl_a04'  ); vol_clx(:,4)=aero(:,kaer)/den_clx
  case( 'cl_a05'  ); vol_clx(:,5)=aero(:,kaer)/den_clx
  case( 'cl_a06'  ); vol_clx(:,6)=aero(:,kaer)/den_clx
  case( 'cl_a07'  ); vol_clx(:,7)=aero(:,kaer)/den_clx
  case( 'cl_a08'  ); vol_clx(:,8)=aero(:,kaer)/den_clx

  case( 'na_a01'  ); vol_nax(:,1)=aero(:,kaer)/den_nax
  case( 'na_a02'  ); vol_nax(:,2)=aero(:,kaer)/den_nax
  case( 'na_a03'  ); vol_nax(:,3)=aero(:,kaer)/den_nax
  case( 'na_a04'  ); vol_nax(:,4)=aero(:,kaer)/den_nax
  case( 'na_a05'  ); vol_nax(:,5)=aero(:,kaer)/den_nax
  case( 'na_a06'  ); vol_nax(:,6)=aero(:,kaer)/den_nax
  case( 'na_a07'  ); vol_nax(:,7)=aero(:,kaer)/den_nax
  case( 'na_a08'  ); vol_nax(:,8)=aero(:,kaer)/den_nax

  case( 'oin_a01' ); vol_oin(:,1)=aero(:,kaer)/den_oin
  case( 'oin_a02' ); vol_oin(:,2)=aero(:,kaer)/den_oin
  case( 'oin_a03' ); vol_oin(:,3)=aero(:,kaer)/den_oin
  case( 'oin_a04' ); vol_oin(:,4)=aero(:,kaer)/den_oin
  case( 'oin_a05' ); vol_oin(:,5)=aero(:,kaer)/den_oin
  case( 'oin_a06' ); vol_oin(:,6)=aero(:,kaer)/den_oin
  case( 'oin_a07' ); vol_oin(:,7)=aero(:,kaer)/den_oin
  case( 'oin_a08' ); vol_oin(:,8)=aero(:,kaer)/den_oin

  case( 'num_a01' ); num_aer(:,1)=aero(:,kaer)            ! #/cm^3
  case( 'num_a02' ); num_aer(:,2)=aero(:,kaer)
  case( 'num_a03' ); num_aer(:,3)=aero(:,kaer)
  case( 'num_a04' ); num_aer(:,4)=aero(:,kaer)
  case( 'num_a05' ); num_aer(:,5)=aero(:,kaer)
  case( 'num_a06' ); num_aer(:,6)=aero(:,kaer)
  case( 'num_a07' ); num_aer(:,7)=aero(:,kaer)
  case( 'num_a08' ); num_aer(:,8)=aero(:,kaer)

  case( 'water_a01' ); vol_h2o(:,1)=aero(:,kaer)/den_h2o
  case( 'water_a02' ); vol_h2o(:,2)=aero(:,kaer)/den_h2o
  case( 'water_a03' ); vol_h2o(:,3)=aero(:,kaer)/den_h2o
  case( 'water_a04' ); vol_h2o(:,4)=aero(:,kaer)/den_h2o
  case( 'water_a05' ); vol_h2o(:,5)=aero(:,kaer)/den_h2o
  case( 'water_a06' ); vol_h2o(:,6)=aero(:,kaer)/den_h2o
  case( 'water_a07' ); vol_h2o(:,7)=aero(:,kaer)/den_h2o
  case( 'water_a08' ); vol_h2o(:,8)=aero(:,kaer)/den_h2o

  case default; cycle
  end select
end do

dlo_sect=3.90625e-6
dhi_sect=10.0e-4

do kk=1, nsig
  do ksz=1, nsize_mosaic_wrf
    vol_dry=vol_so4(kk,ksz) + vol_no3(kk,ksz) + &
            vol_nh4(kk,ksz) + vol_ocx(kk,ksz) + &
            vol_bcx(kk,ksz) + vol_clx(kk,ksz) + &
            vol_nax(kk,ksz) + vol_oin(kk,ksz)
    vol_wet=vol_h2o(kk,ksz) + vol_dry

    !!! total_aod ----------------------------------------
    vol_lay_dryaer(kk,ksz)=vol_dry
    vol_lay_wetaer(kk,ksz)=vol_wet
    !!! total_aod ----------------------------------------

    num_aer_lo=1.90985*vol_dry/(dlo_sect**3)
    num_aer_hi=1.90985*vol_dry/(dhi_sect**3)

    if( num_aer(kk,ksz)>num_aer_lo ) then
      num_aer(kk,ksz)=num_aer_lo
    else if( num_aer(kk,ksz)<num_aer_hi ) then
      num_aer(kk,ksz)=num_aer_hi
    end if

    !!! total_aod ----------------------------------------
    if( vol_dry<=1.e-15 .and. num_aer(kk,ksz)<=1.e-10 ) then
      dp_wet=dhi_sect
      rad_wet(kk,ksz)=dp_wet/2.     ! cm
    else
      dp_wet=(1.90985*vol_wet/num_aer(kk,ksz))**0.3333333
      rad_wet(kk,ksz)=dp_wet/2.     ! cm
      dp_dry=(1.90985*vol_dry/num_aer(kk,ksz))**0.3333333
      rad_dry(kk,ksz)=dp_dry/2.     ! cm
    end if
    !!! total_aod ----------------------------------------

    do ksw=1, nswbands
      ri_ave=swref_index_so4(ksw)*vol_so4(kk,ksz) + swref_index_no3(ksw)*vol_no3(kk,ksz) + &
             swref_index_nh4(ksw)*vol_nh4(kk,ksz) + swref_index_ocx(ksw)*vol_ocx(kk,ksz) + &
             swref_index_bcx(ksw)*vol_bcx(kk,ksz) + swref_index_clx(ksw)*vol_clx(kk,ksz) + &
             swref_index_nax(ksw)*vol_nax(kk,ksz) + swref_index_oin(ksw)*vol_oin(kk,ksz) + &
             swref_index_h2o(ksw)*vol_h2o(kk,ksz)
      ri_ave=ri_ave/vol_wet
      swrefindx_wet(kk,ksz,ksw)=ri_ave

      ri_ave=swref_index_so4(ksw)*vol_so4(kk,ksz) + swref_index_no3(ksw)*vol_no3(kk,ksz) + &
             swref_index_nh4(ksw)*vol_nh4(kk,ksz) + swref_index_ocx(ksw)*vol_ocx(kk,ksz) + &
             swref_index_bcx(ksw)*vol_bcx(kk,ksz) + swref_index_clx(ksw)*vol_clx(kk,ksz) + &
             swref_index_nax(ksw)*vol_nax(kk,ksz) + swref_index_oin(ksw)*vol_oin(kk,ksz)
      ri_ave=ri_ave/vol_dry
      swrefindx_dry(kk,ksz,ksw)=ri_ave
    end do

    rden_mix(kk,ksz)= den_so4*vol_so4(kk,ksz) + den_no3*vol_no3(kk,ksz) + &
                      den_nh4*vol_nh4(kk,ksz) + den_ocx*vol_ocx(kk,ksz) + &
                      den_bcx*vol_bcx(kk,ksz) + den_clx*vol_clx(kk,ksz) + &
                      den_nax*vol_nax(kk,ksz) + den_oin*vol_oin(kk,ksz)
    rden_mix(kk,ksz)= rden_mix(kk,ksz)/vol_dry

  end do
end do

!----------------------------
swtauaer=0.
swssaaer=0.
swgymaer=0.

xrmin=log(rmin)
xrmax=log(rmax)
dokk: do kk=1, nsig
  thesum=0.
  do ksz=1, nsize_mosaic_wrf
    thesum=thesum+num_aer(kk,ksz)
  end do
  if( thesum<1.e-21 ) cycle
  doksw: do ksw=1, nswbands
    doksz: do ksz=1, nsize_mosaic_wrf

      refc=swrefindx_wet(kk,ksz,ksw)
      refr=real(refc)
      refi=-imag(refc)
      if( rad_wet(kk,ksz)<rmin ) rad_wet(kk,ksz)=rmin
      if( rad_wet(kk,ksz)>rmax ) rad_wet(kk,ksz)=rmax
      !!! total_aod ----------------------------------------
      if( rad_dry(kk,ksz)<rmin ) rad_dry(kk,ksz)=rmin
      if( rad_dry(kk,ksz)>rmax ) rad_dry(kk,ksz)=rmax
      !!! total_aod ----------------------------------------
      xlog=log( rad_wet(kk,ksz) )
      xrad=(2*xlog-xrmax-xrmin)/(xrmax-xrmin)

      itab=0
      !            table(km,im,jm)   km    im    jm    x    y    xtab(im)         ytab(jm)         ix   jy   t    tu   out(km)
      call binterp(extpsw(1,1,1,ksw),ncoef,nrefr,nrefi,refr,refi,refrtabsw(1,ksw),refitabsw(1,ksw),itab,jtab,ttab,utab,cext)
      call binterp(scapsw(1,1,1,ksw),ncoef,nrefr,nrefi,refr,refi,refrtabsw(1,ksw),refitabsw(1,ksw),itab,jtab,ttab,utab,csca)
      call binterp(asmpsw(1,1,1,ksw),ncoef,nrefr,nrefi,refr,refi,refrtabsw(1,ksw),refitabsw(1,ksw),itab,jtab,ttab,utab,casm)

      ! chebyshev polynomials
      ch(1)=1.
      ch(2)=xrad
      do knc=3, ncoef
        ch(knc)=2.*xrad*ch(knc-1)-ch(knc-2)
      end do

      pext=0.5*cext(1)
      do knc=2, ncoef
        pext=pext+ch(knc)*cext(knc)
      end do
      pext=exp(pext)

      psca=0.5*csca(1)
      do knc=2, ncoef
        psca=psca+ch(knc)*csca(knc)
      end do
      psca=exp(psca)

      psca=min(psca,pext)
      pabs=pext-psca

      !!! cwy mosaic ------------------------------------------------------------
      Rup_Rlow=refrtabsw(itab+1,ksw)-refrtabsw(itab,ksw)
      Iup_Ilow=refitabsw(jtab+1,ksw)-refitabsw(jtab,ksw)

      ! for assimilating absorbing aod
      a00=0.; a01=0.; a10=0.; a11=0.
      uu=(refr-refrtabsw(itab,ksw))/Rup_Rlow
      if( uu>1 ) print *, 'cwy warning uu>1, refr=', refr, refrtabsw(itab,ksw), Rup_Rlow
      do kch=1, ncoef
        a00=a00+ch(kch)*(extpsw(kch,itab  ,jtab  ,ksw)-scapsw(kch,itab  ,jtab  ,ksw))
        a01=a01+ch(kch)*(extpsw(kch,itab  ,jtab+1,ksw)-scapsw(kch,itab  ,jtab+1,ksw))
        a10=a10+ch(kch)*(extpsw(kch,itab+1,jtab  ,ksw)-scapsw(kch,itab+1,jtab  ,ksw))
        a11=a11+ch(kch)*(extpsw(kch,itab+1,jtab+1,ksw)-scapsw(kch,itab+1,jtab+1,ksw))
      end do
      a00=a00*pabs; a01=a01*pabs
      a10=a10*pabs; a11=a11*pabs
      betaI(kk,ksz,ksw)=((uu-1)*a00-uu*a01+(1-uu)*a10+uu*a11)/Iup_Ilow

      ! for assimilating scattering aod
      a00=0.; a01=0.; a10=0.; a11=0.
      vv=(refi-refitabsw(jtab,ksw))/Iup_Ilow
      if( vv>1 ) print *, 'cwy warning vv>1, refi=', refi, refitabsw(jtab,ksw), Iup_Ilow
      do kch=1, ncoef
        a00=a00+ch(kch)*scapsw(kch,itab  ,jtab  ,ksw)
        a01=a01+ch(kch)*scapsw(kch,itab  ,jtab+1,ksw)
        a10=a10+ch(kch)*scapsw(kch,itab+1,jtab  ,ksw)
        a11=a11+ch(kch)*scapsw(kch,itab+1,jtab+1,ksw)
      end do
      a00=a00*psca; a01=a01*psca
      a10=a10*psca; a11=a11*psca
      betaR(kk,ksz,ksw)=((vv-1)*a00+(1-vv)*a01-vv*a10+vv*a11)/Rup_Rlow

      eext_lay(kk,ksz,ksw)=(psca+pabs)*pie*exp(xlog)**2 ! cm^2
      !!! cwy mosaic ------------------------------------------------------------

      pasm=0.5*casm(1)
      do knc=2, ncoef
        pasm=pasm+ch(knc)*casm(knc)
      end do
      pasm=exp(pasm)

      psca=min(psca,pext)
      weighte(kk,ksz,ksw)=pext*pie*exp(xlog)**2       ! cm^2, pext: extinction coef. weighte: extinction cross-section, L^2 #^-1
      weights(kk,ksz,ksw)=psca*pie*exp(xlog)**2       ! weighte: scattering cross section

      ! swtauaer in L^-1, volume extinction coefficient; num_aer in # L^-3
      swtauaer(kk,ksw)=swtauaer(kk,ksw)+weighte(kk,ksz,ksw)*num_aer(kk,ksz) ! cm^2 * #/cm^3 ->  cm^-1
      swssaaer(kk,ksw)=swssaaer(kk,ksw)+weights(kk,ksz,ksw)*num_aer(kk,ksz)
      swgymaer(kk,ksw)=swgymaer(kk,ksw)+weights(kk,ksz,ksw)*num_aer(kk,ksz)*pasm
    end do doksz

    swgymaer(kk,ksw)=swgymaer(kk,ksw)/swssaaer(kk,ksw)
    swextaer(kk,ksw)=swtauaer(kk,ksw)*1.e+5           ! km^-1
    swssaaer(kk,ksw)=swssaaer(kk,ksw)/swtauaer(kk,ksw)
  end do doksw
end do dokk

do ksw=1, nswbands
do kk=1, nsig
  swtauaer(kk,ksw)=swtauaer(kk,ksw)*zhgt(kk)*100.     ! zhgt in meter, swtauaer in dimensionless
end do
end do

!----------------------------
!real(r_kind),dimension(nsig,nchanl),intent(out)        ::  layer_od
!real(r_kind),dimension(nsigaerojac,nchanl),intent(out) ::  jacobian_aero
layer_od=zero
jacobian_aero=zero

!!! cwy mosaic ------------------------------------------------------------
dokch: do kch=1, nchanl
  if( data_s(nreal+kch)>zero ) then
    if( typename=='ce318_aod' ) then
      if( kch==1 ) then; layer_od(:,kch)=swtauaer(:,1); ksw=1; endif
      if( kch==2 ) then; layer_od(:,kch)=swtauaer(:,3); ksw=3; endif
      if( kch==3 ) then; layer_od(:,kch)=swtauaer(:,4); ksw=4; endif
      if( kch==4 ) then; layer_od(:,kch)=swtauaer(:,5); ksw=5; endif
    end if
    if( typename=='modis_aod' ) then
      if( kch==4 ) then; layer_od(:,kch)=swtauaer(:,2); ksw=2; endif
    end if

    totaod_channel=sum(layer_od(:,kch))

    dokctr: do kctr=1, n_aerosols_jac     ! aerosoljacnames(1:njac)==aerosol_names_jac(1:njac)
      indx=aerojacindxs(kctr)             ! jacobian aerosols should be in the order shown in aerojacindxs
      ksz=-1
      if( index(aerosol_names_jac(kctr),'a01')>0 ) ksz=1
      if( index(aerosol_names_jac(kctr),'a02')>0 ) ksz=2
      if( index(aerosol_names_jac(kctr),'a03')>0 ) ksz=3
      if( index(aerosol_names_jac(kctr),'a04')>0 ) ksz=4
      if( index(aerosol_names_jac(kctr),'a05')>0 ) ksz=5
      if( index(aerosol_names_jac(kctr),'a06')>0 ) ksz=6
      if( index(aerosol_names_jac(kctr),'a07')>0 ) ksz=7
      if( index(aerosol_names_jac(kctr),'a08')>0 ) ksz=8

      if( n_cvaer/=0 ) then
      if( index(aerosol_names_jac(kctr),'so4')>0 ) then; refr=refrsw_so4(ksw); refi=refisw_so4(ksw)*-1; rden=den_so4; endif
      if( index(aerosol_names_jac(kctr),'no3')>0 ) then; refr=refrsw_no3(ksw); refi=refisw_no3(ksw)*-1; rden=den_no3; endif
      if( index(aerosol_names_jac(kctr),'nh4')>0 ) then; refr=refrsw_nh4(ksw); refi=refisw_nh4(ksw)*-1; rden=den_nh4; endif
      if( index(aerosol_names_jac(kctr),'oc' )>0 ) then; refr=refrsw_ocx(ksw); refi=refisw_ocx(ksw)*-1; rden=den_ocx; endif
      if( index(aerosol_names_jac(kctr),'bc' )>0 ) then; refr=refrsw_bcx(ksw); refi=refisw_bcx(ksw)*-1; rden=den_bcx; endif
      if( index(aerosol_names_jac(kctr),'oin')>0 ) then; refr=refrsw_oin(ksw); refi=refisw_oin(ksw)*-1; rden=den_oin; endif
      if( index(aerosol_names_jac(kctr),'cl' )>0 ) then; refr=refrsw_clx(ksw); refi=refisw_clx(ksw)*-1; rden=den_clx; endif
      if( index(aerosol_names_jac(kctr),'na' )>0 ) then; refr=refrsw_nax(ksw); refi=refisw_nax(ksw)*-1; rden=den_nax; endif
      end if

      if( ksz>0 ) then
        dolev: do kk=1, nsig
          !ugkg_gcm3(kk)=den_air(kk)*1.e-12                       ! ug/kg -> g/cm^3
          !nmkg_ncm3(kk)=den_air(kk)*1.e-6                        !  #/kg -> #/cm^3
          !
          ! jac for num_a0x
          if( n_cvnum/=0 .and. index(aerosol_names_jac(kctr),'num_a')>0 ) then
            if( totaod_channel>rmin_value .and. vol_lay_wetaer(kk,ksz)>rmin_value ) then
              jacobian_aero(indx+kk,kch)=eext_lay(kk,ksz,ksw)
              jacobian_aero(indx+kk,kch)=jacobian_aero(indx+kk,kch)*nmkg_ncm3(kk)*zhgt(kk)/sum(zhgt)*1.e+8
            end if
          end if
          ! jac for mas_a0x
          if( n_cvmas/=0 .and. index(aerosol_names_jac(kctr),'mas_a')>0 ) then
            if( totaod_channel>rmin_value .and. vol_lay_wetaer(kk,ksz)>rmin_value ) then
              refr= real(swrefindx_wet(kk,ksz,ksw))
              refi=-imag(swrefindx_wet(kk,ksz,ksw))
              rden=rden_mix(kk,ksz)

              jacobian_aero(indx+kk,kch)=3*eext_lay(kk,ksz,ksw)/(4*pie*rad_dry(kk,ksz)**3*rden)
              jacobian_aero(indx+kk,kch)=jacobian_aero(indx+kk,kch)+ &
              & num_aer(kk,ksz)*pie*rad_wet(kk,ksz)**2*refr/(rden*vol_lay_wetaer(kk,ksz))*betaR(kk,ksz,ksw)
              jacobian_aero(indx+kk,kch)=jacobian_aero(indx+kk,kch)+ &
              & num_aer(kk,ksz)*pie*rad_wet(kk,ksz)**2*refi/(rden*vol_lay_wetaer(kk,ksz))*betaI(kk,ksz,ksw)
              jacobian_aero(indx+kk,kch)=jacobian_aero(indx+kk,kch)*ugkg_gcm3(kk)*zhgt(kk)/sum(zhgt)*1.e+8
            end if
          end if
          ! jac for aer_a0x
          if( n_cvaer/=0 .and. index(aerosol_names_jac(kctr),'num_a')==0 ) then
            if( totaod_channel>rmin_value .and. vol_lay_wetaer(kk,ksz)>rmin_value ) then
              jacobian_aero(indx+kk,kch)=3*eext_lay(kk,ksz,ksw)/(4*pie*rad_dry(kk,ksz)**3*rden)
              jacobian_aero(indx+kk,kch)=jacobian_aero(indx+kk,kch)+&
              & num_aer(kk,ksz)*pie*rad_wet(kk,ksz)**2*refr/(rden*vol_lay_wetaer(kk,ksz))*betaR(kk,ksz,ksw)
              jacobian_aero(indx+kk,kch)=jacobian_aero(indx+kk,kch)+&
              & num_aer(kk,ksz)*pie*rad_wet(kk,ksz)**2*refi/(rden*vol_lay_wetaer(kk,ksz))*betaI(kk,ksz,ksw)
              jacobian_aero(indx+kk,kch)=jacobian_aero(indx+kk,kch)*ugkg_gcm3(kk)*zhgt(kk)/sum(zhgt)*1.e+8
            end if
          end if
        end do dolev
      else
        write(*,*) 'cwy: ksz<0; stop; ksz, aerosol_names_jac(kctr)=',ksz,aerosol_names_jac(kctr)
        stop
      end if
    end do dokctr
  end if
end do dokch
!!! cwy mosaic ------------------------------------------------------------
end subroutine call_aod_mosaic

!=========================================================================
include './sub_mieaer.f90'

end module aod_mosaic
