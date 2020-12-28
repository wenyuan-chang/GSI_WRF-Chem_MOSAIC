module aerosol_water_awa
use kinds, only : r_kind, i_kind

contains
!!! cwy ----------------------- aerosol water content -------------------
real function aerosol_water( rso4, rno3, rclx, rocx, rbcx, roin, rhum )
real(r_kind), intent( in )          ::  rso4, rno3, rclx, rocx, rbcx, roin, rhum
integer(i_kind)                     ::  kaer
real(r_kind)                        ::  aw, aw_min, xm, bin_molality, b_zsr
real(r_kind)                        ::  dum, tmpa, aermol, aermas, aerden, kappa
real(r_kind), dimension(6)          ::  a_zsr

! a_zsr, b_zsr, aw_min data from WRF_Chem/chem/module_mosaic_therm.F
dum=0.; tmpa=0.
do kaer=1, 6
  if( kaer==1 ) then ! so4 follows (nh4)2so4, nmol/m3
    aw_min=0.1
    a_zsr(1)=  1.30894
    a_zsr(2)= -7.09922
    a_zsr(3)= 20.62831
    a_zsr(4)=-32.19965
    a_zsr(5)= 25.17026
    a_zsr(6)= -7.81632
    b_zsr   = 28.0811
    aermol  =rso4
  end if
  if( kaer==2 ) then ! no3 follows nh4no3, nmol/m3
    aw_min=0.1
    a_zsr(1)=  0.43507
    a_zsr(2)=  6.38220
    a_zsr(3)=-30.19797
    a_zsr(4)= 53.36470
    a_zsr(5)=-43.44203
    a_zsr(6)= 13.46158
    b_zsr   = 33.4049
    aermol  =rno3
  end if
  if( kaer==3 ) then ! clx follows nh4cl, nmol/m3
    aw_min=0.1
    a_zsr(1)=  0.45309
    a_zsr(2)=  2.65606
    a_zsr(3)=-14.7730
    a_zsr(4)= 26.2936
    a_zsr(5)=-20.5735
    a_zsr(6)=  5.94255
    b_zsr   = 30.8888
    aermol  =rclx
  end if
  if( kaer==4 ) then ! ocx, aermas: ng/m3
    aw_min=0.1
    aermas=rocx
    aerden=1.0       ! g/cm3
    kappa =1.e-4
  end if
  if( kaer==5 ) then ! bcx, aermas: ng/m3
    aw_min=0.1
    aermas=rbcx
    aerden=1.7       ! g/cm3
    kappa =1.e-6
  end if
  if( kaer==6 ) then ! oin, aermas: ng/m3
    aw_min=0.1
    aermas=roin
    aerden=2.6       ! g/cm3
    kappa =0.068
  end if

  aw=max(rhum,aw_min)
  aw=min(aw,0.999999)

  if( kaer<=3 ) then
    if( aw<0.97 ) then
      xm=a_zsr(1)
      xm=xm+a_zsr(2)*aw
      xm=xm+a_zsr(3)*aw**2
      xm=xm+a_zsr(4)*aw**3
      xm=xm+a_zsr(5)*aw**4
      xm=xm+a_zsr(6)*aw**5
      bin_molality=55.509*xm/(1.-xm)     ! bin_molality: mol/kg
    else
      bin_molality=-b_zsr*log(aw)        ! bin_molality: mol/kg
    end if
    dum=dum+aermol/bin_molality        ! ug/m3
  else
    tmpa=tmpa+aermas/aerden*kappa        ! n(cm3)/m3=ncc/m3
  end if
end do
dum=dum+1.e-3*tmpa*aw/(1.-aw)            ! ug/m3, tmpa*water_density(=1g/cm3)*1.e-3
aerosol_water=dum*1.e-9                  ! kg(water)/m^3(air)

end function aerosol_water
!!! cwy ----------------------- aerosol water content -------------------

end module aerosol_water_awa

