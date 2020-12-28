module set_crtm_aerosolmod
!$$$ module documentation block
!           .      .    .                                       .
! module:   set_crtm_aerosolmod
!  prgmmr: todling          org: gmao                date: 2011-06-01
!
! abstract: module providing interface to set-crtm-aerosol procedures
!
! program history log:
!   2011-06-01  todling
!   2011-09-20  hclin   - separate na and na_crtm for p25 handling
!
! subroutines included:
!   sub Set_CRTM_Aerosol_
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

private

public Set_CRTM_Aerosol

contains

!!!interface Set_CRTM_Aerosol
!!!  subroutine Set_CRTM_Aerosol_ ( km, na, na_crtm, aero_name, aero_conc, rh, aerosol)
!!! add ---------------------------
  subroutine Set_CRTM_Aerosol ( km, na, na_crtm, aero_name, aero_conc, rh, tk, aerosol)
  use kinds, only: i_kind,r_kind
  use constants, only: tiny_r_kind
  use mpimod, only: mype
  use CRTM_Aerosol_Define, only: CRTM_Aerosol_type
  use mpeu_util, only: getindex
  use crtm_module, only: SULFATE_AEROSOL,BLACK_CARBON_AEROSOL,ORGANIC_CARBON_AEROSOL,&
      DUST_AEROSOL,SEASALT_SSAM_AEROSOL,SEASALT_SSCM1_AEROSOL,SEASALT_SSCM2_AEROSOL,SEASALT_SSCM3_AEROSOL
  implicit none
  integer(i_kind) , intent(in)    :: km                ! number of levels
  integer(i_kind) , intent(in)    :: na                ! number of aerosols
  integer(i_kind) , intent(in)    :: na_crtm           ! number of aerosols seen by CRTM
  character(len=*), intent(in)    :: aero_name(na)     ! [na]    GOCART aerosol names: du0001, etc.
  real(r_kind),     intent(inout) :: aero_conc(km,na)  ! [km,na] aerosol concentration (Kg/m2)
  type(CRTM_Aerosol_type), intent(inout) :: aerosol(na_crtm)! [na]   CRTM Aerosol object

  real(r_kind),     intent(in)    :: rh(km)            ! [km]    relative humdity [0,1]
  real(r_kind),     intent(in)    :: tk(km)            ! [km]    temperature [K]
  
  integer       :: kna, hygro_B, kk
  real(r_kind)  :: rhoh2o, zB, zmw, zR, zdsig, zda, zdb, zdc, zdd, zde, zq4, zq3, zq2, zq1, zq0
  real(r_kind)  :: zdrd, zdrsig, zdrw, zreff, rhval

  do kna=1, na_crtm
    aerosol(kna)%is_allocated=.true.
    aerosol(kna)%max_layers=km
    aerosol(kna)%n_layers=km

    hygro_B=0
    if( trim(aero_name(kna))=='sulf' ) then; aerosol(kna)%type=SULFATE_AEROSOL;        hygro_B=1; zdrd=0.0695; zdrsig=2.03; endif
    if( trim(aero_name(kna))=='bc1'  ) then; aerosol(kna)%type=BLACK_CARBON_AEROSOL;   hygro_B=2; zdrd=0.0118; zdrsig=2.00; endif
    if( trim(aero_name(kna))=='bc2'  ) then; aerosol(kna)%type=BLACK_CARBON_AEROSOL;   hygro_B=2; zdrd=0.0118; zdrsig=2.00; endif
    if( trim(aero_name(kna))=='oc1'  ) then; aerosol(kna)%type=ORGANIC_CARBON_AEROSOL; hygro_B=2; zdrd=0.0212; zdrsig=2.20; endif
    if( trim(aero_name(kna))=='oc2'  ) then; aerosol(kna)%type=ORGANIC_CARBON_AEROSOL; hygro_B=2; zdrd=0.0212; zdrsig=2.20; endif

    if( trim(aero_name(kna))=='dust1') aerosol(kna)%type=DUST_AEROSOL
    if( trim(aero_name(kna))=='dust2') aerosol(kna)%type=DUST_AEROSOL
    !------
    if( trim(aero_name(kna))=='dust3') aerosol(kna)%type=DUST_AEROSOL
    if( trim(aero_name(kna))=='dust4') aerosol(kna)%type=DUST_AEROSOL
    if( trim(aero_name(kna))=='dust5') aerosol(kna)%type=DUST_AEROSOL
    !------

    if( trim(aero_name(kna))=='seas1') aerosol(kna)%type=SEASALT_SSAM_AEROSOL
    if( trim(aero_name(kna))=='seas2') aerosol(kna)%type=SEASALT_SSAM_AEROSOL
    !------
    if( trim(aero_name(kna))=='seas3') aerosol(kna)%type=SEASALT_SSAM_AEROSOL
    if( trim(aero_name(kna))=='seas4') aerosol(kna)%type=SEASALT_SSAM_AEROSOL
    !------

    aerosol(kna)%concentration(:)=aero_conc(1:km,kna)

    !------ hygroscoicity
    rhoh2o=1.                                             !!! water density in g/cm3
    if( hygro_B/=0 ) then
      do kk=1, km
        rhval=rh(kk)
        if( rhval>0.99 ) rhval=0.99
        if( rhval<0.20 ) rhval=0.20

        select case( hygro_B )
        case( 1 )
          zB=2.42848-3.85261*rhval+1.88159*rhval*rhval !!! rh: 0 to 1
        case( 2 )
          zB=5.e-7
        end select

        zmw=18.*1.e-3                                     !!! kg/mol
        zR=8.3144621                                      !!! J/mol/K
        zdsig=0.0761-1.55*1.e-4*(tk(kk)-273)              !!! J/m2, water surface tension
        zda=log(rhval)                                    !!! rh: 0 to 1
        zdb=2*zmw*zdsig*1.e-6                             !!! *1.e-6 for drzw in m
        zdc=zR*tk(kk)*rhoh2o*1.e+3                        !!! rhoh2o*1.e+3 for water density in kg/m3
        zdd=zB*zdrd**3                                    !!! zdrd in um
        zde=zdrd**3

        zq4=zda*zdc
        zq3=-1*zdb
        zq2=0.
        zq1=zdc*zdd-zda*zdc*zde
        zq0=zdb*zde
       
        call quartic( zq4, zq3, zq2, zq1, zq0, zdrd, zdrw )
        zreff=zdrw*exp(2.5*(log(zdrsig))**2)               !!! um
        aerosol(kna)%effective_radius(kk)=zreff*1.e-3      !!! um -> mm
      end do
    else
      aerosol(kna)%effective_radius(1:km)=2.5*1.e-3        !!! um -> mm
    end if
  end do
!!! add ---------------------------

  end subroutine Set_CRTM_Aerosol
!!!  end subroutine Set_CRTM_Aerosol_
!!!end interface

subroutine quartic( q4, q3, q2, q1, q0, x0, xx )
implicit none
double precision, intent( in )          ::      q4, q3, q2, q1, q0, x0
double precision, intent( out )         ::      xx  
double precision                        ::      func, fdev, xn, xn1, eps 
integer                                 ::      knum, xsta
eps=1.e-4
xn=x0
knum=0.
xsta=0.
do
  fdev=4*q4*xn**3+3*q3*xn**2+2*q2*xn+q1
  if( fdev>0 ) then
    xsta=xsta+1
    xn=xsta
    cycle
  end if
  func=q4*xn**4+q3*xn**3+q2*xn**2+q1*xn+q0
  xn1=xn-func/fdev
  if( abs(xn1-xn)<=eps .and. xn1>0 .and. xn1<100 ) exit
  xn=(xn1+xn)/2.
  knum=knum+1
  if( knum>500 ) then
    write(*,*) 'cwy warning: too many iterations. knum>xxx. RH range may not be valid.'
    write(*,*) 'cwy warning: STOP: q4|q3|q2|q1|q0=', q4, q3, q2, q1, q0
    write(*,*) 'cwy warning: xn, xn1=', xn, xn1; stop 99
  end if
end do
xx=xn
end subroutine quartic

end module set_crtm_aerosolmod
