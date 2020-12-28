module stppm2_5mod
  
!$$$ module documentation block
!           .      .    .                                       .
! module:   stppm2_5mod    module for stppm2_5 and its tangent linear stppm2_5_tl
!  pgrmmr:
!
! abstract: module for stppm2_5 and its tangent linear stppm2_5_tl
!
! program history log:
!   2010-10-15  pagowski - initial code for pm2_5
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub stppm2_5
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
  use m_obsNode  , only: obsNode
  use m_pm2_5Node, only: pm2_5Node
  use m_pm2_5Node, only: pm2_5Node_typecast
  use m_pm2_5Node, only: pm2_5Node_nextcast
  implicit none

  private
  public stppm2_5
  
contains
  
  subroutine stppm2_5(pm2_5head,rval,sval,out,sges,nstep)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    stppm2_5        calcuate penalty and stepsize from q
!                            with addition of nonlinear qc.
!   prgmmr: derber           org: np23                date: 1991-02-26
!
! abstract: calculate penalty and contribution to stepsize from pm2_5
!           using nonlinear qc.
!
! program history log:
!   2010-10-15  pagowski - initial code for pm2_5
!
!   input argument list:
!     pm2_5head
!     rpm2_5       - search direction for pm2_5
!     spm2_5       - analysis increment for pm2_5
!     sges     - stepsize estimates (nstep)
!     nstep    - number of stepsize estimates (== 0 means use outer iteration values)
!
!   output argument list:
!     out(1:nstep)   - contribution of penalty from pm2_5 sges(1:nstep)
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
    use kinds, only: r_kind,i_kind,r_quad
    use qcmod, only: nlnqc_iter,varqc_iter
    use constants, only: half,one,two,tiny_r_kind,cg_term,zero_quad,max_varname_length,zero
    use gsi_bundlemod, only: gsi_bundle
    use gsi_bundlemod, only: gsi_bundlegetpointer
    use gridmod, only: cmaq_regional,wrf_mass_regional
    USE chemmod, ONLY: d_2_5,s_2_5,nh4_mfac,oc_mfac,laeroana_gocart ,&
!!! cwy mosaic4 ---------------------------------------------------- chem
                       laeroana_mosaic4
!!! cwy mosaic4 ---------------------------------------------------- chem
    implicit none
    
! declare passed variables
    class(obsNode), pointer             ,intent(in   ) :: pm2_5head
    integer(i_kind)                     ,intent(in   ) :: nstep
    real(r_quad),dimension(max(1,nstep)),intent(inout) :: out
    type(gsi_bundle)                    ,intent(in   ) :: rval,sval
    real(r_kind),dimension(max(1,nstep)),intent(in   ) :: sges
    
! declare local variables
    integer(i_kind) ier,istatus
    integer(i_kind) j1,j2,j3,j4,j5,j6,j7,j8,kk
    real(r_kind) cg_pm2_5,val,val2,wgross,wnotgross,pm2_5_pg
    real(r_kind),dimension(max(1,nstep))::pen
    real(r_kind) w1,w2,w3,w4,w5,w6,w7,w8,qq
    real(r_kind),pointer,dimension(:):: rpm2_5,spm2_5
    type(pm2_5Node), pointer :: pm2_5ptr

    character(len=max_varname_length) :: aeroname    

    out=zero_quad

!   If no pm2_5 data return
    if(.not. associated(pm2_5head))return

!!! cwy mosaic4 ---------------------------------------------------- chem
    if (cmaq_regional .or. (wrf_mass_regional .and. .not.laeroana_gocart .and. .not.laeroana_mosaic4)) then
!!! cwy mosaic4 ---------------------------------------------------- chem
    
       ier=0
       call gsi_bundlegetpointer(sval,'pm2_5',spm2_5,istatus);ier=ier+istatus
       call gsi_bundlegetpointer(rval,'pm2_5',rpm2_5,istatus);ier=ier+istatus

       if(ier/=0) return

       !pm2_5ptr => pm2_5head
       pm2_5ptr => pm2_5Node_typecast(pm2_5head)
       do while (associated(pm2_5ptr))
          if(pm2_5ptr%luse)then
             if(nstep > 0)then
                j1=pm2_5ptr%ij(1)
                j2=pm2_5ptr%ij(2)
                j3=pm2_5ptr%ij(3)
                j4=pm2_5ptr%ij(4)
                j5=pm2_5ptr%ij(5)
                j6=pm2_5ptr%ij(6)
                j7=pm2_5ptr%ij(7)
                j8=pm2_5ptr%ij(8)
                w1=pm2_5ptr%wij(1)
                w2=pm2_5ptr%wij(2)
                w3=pm2_5ptr%wij(3)
                w4=pm2_5ptr%wij(4)
                w5=pm2_5ptr%wij(5)
                w6=pm2_5ptr%wij(6)
                w7=pm2_5ptr%wij(7)
                w8=pm2_5ptr%wij(8)


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)
                val2=w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)-pm2_5ptr%res


                do kk=1,nstep
                   qq=val2+sges(kk)*val
                   pen(kk)=qq*qq*pm2_5ptr%err2
                end do
             else
                pen(1)=pm2_5ptr%res*pm2_5ptr%res*pm2_5ptr%err2
             end if

!  modify penalty term if nonlinear qc

             if (nlnqc_iter .and. pm2_5ptr%pg > tiny_r_kind .and. &
                  pm2_5ptr%b  > tiny_r_kind) then
                pm2_5_pg=pm2_5ptr%pg*varqc_iter
                cg_pm2_5=cg_term/pm2_5ptr%b
                wnotgross= one-pm2_5_pg
                wgross = pm2_5_pg*cg_pm2_5/wnotgross
                do kk=1,max(1,nstep)
                   pen(kk)= -two*log((exp(-half*pen(kk))+wgross)/(one+wgross))
                end do
             endif

             out(1) = out(1)+pen(1)*pm2_5ptr%raterr2
             do kk=2,nstep
                out(kk) = out(kk)+(pen(kk)-pen(1))*pm2_5ptr%raterr2
             end do
          end if

          !pm2_5ptr => pm2_5ptr%llpoint
          pm2_5ptr => pm2_5Node_nextcast(pm2_5ptr)

       end do

    endif

    if (wrf_mass_regional .and. laeroana_gocart) then

       !pm2_5ptr => pm2_5head
       pm2_5ptr => pm2_5Node_typecast(pm2_5head)

       do while (associated(pm2_5ptr))
          if(pm2_5ptr%luse)then
             if(nstep > 0)then
                j1=pm2_5ptr%ij(1)
                j2=pm2_5ptr%ij(2)
                j3=pm2_5ptr%ij(3)
                j4=pm2_5ptr%ij(4)
                j5=pm2_5ptr%ij(5)
                j6=pm2_5ptr%ij(6)
                j7=pm2_5ptr%ij(7)
                j8=pm2_5ptr%ij(8)
                w1=pm2_5ptr%wij(1)
                w2=pm2_5ptr%wij(2)
                w3=pm2_5ptr%wij(3)
                w4=pm2_5ptr%wij(4)
                w5=pm2_5ptr%wij(5)
                w6=pm2_5ptr%wij(6)
                w7=pm2_5ptr%wij(7)
                w8=pm2_5ptr%wij(8)

                istatus=0

                val2=-pm2_5ptr%res
                val=zero

                aeroname='bc1'

                call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(456)
                endif

                call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)

                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(457)
                endif


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)+val
                val2=w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)+val2


                nullify(spm2_5,rpm2_5)


                aeroname='bc2'

                call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(456)
                endif

                call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)

                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(457)
                endif


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)+val
                val2=w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)+val2


                nullify(spm2_5,rpm2_5)


                aeroname='sulf'

                call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(456)
                endif

                call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)

                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(457)
                endif


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)+val
                val2=(w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)) * nh4_mfac +val2


                nullify(spm2_5,rpm2_5)


                aeroname='p25'

                call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(456)
                endif

                call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)

                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(457)
                endif


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)+val
                val2=w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)+val2


                nullify(spm2_5,rpm2_5)


                aeroname='oc1'

                call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(456)
                endif

                call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)

                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(457)
                endif


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)+val
                val2=(w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)) * oc_mfac +val2


                nullify(spm2_5,rpm2_5)

                aeroname='oc2'

                call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(456)
                endif

                call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)

                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(457)
                endif


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)+val
                val2=(w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)) * oc_mfac +val2


                nullify(spm2_5,rpm2_5)


                aeroname='seas1'

                call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(456)
                endif

                call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)

                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(457)
                endif


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)+val
                val2=w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)+val2


                nullify(spm2_5,rpm2_5)


                aeroname='seas2'

                call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(456)
                endif

                call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)

                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(457)
                endif


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)+val
                val2=(w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)) * s_2_5 +val2


                nullify(spm2_5,rpm2_5)


                aeroname='dust1'

                call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(456)
                endif

                call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)

                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(457)
                endif


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)+val
                val2=w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)+val2


                nullify(spm2_5,rpm2_5)


                aeroname='dust2'

                call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(456)
                endif

                call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)

                if(istatus /= 0) then
                   write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for ',&
                        aeroname
                   call stop2(457)
                endif


                val= w1* rpm2_5(j1)+w2* rpm2_5(j2)+w3* rpm2_5(j3)+w4* rpm2_5(j4)+ &
                     w5* rpm2_5(j5)+w6* rpm2_5(j6)+w7* rpm2_5(j7)+w8* rpm2_5(j8)+val
                val2=(w1* spm2_5(j1)+w2* spm2_5(j2)+w3* spm2_5(j3)+w4* spm2_5(j4)+ &
                     w5* spm2_5(j5)+w6* spm2_5(j6)+w7* spm2_5(j7)+w8* spm2_5(j8)) * d_2_5 +val2


                nullify(spm2_5,rpm2_5)

                do kk=1,nstep
                   qq=val2+sges(kk)*val
                   pen(kk)=qq*qq*pm2_5ptr%err2
                end do
             else
                pen(1)=pm2_5ptr%res*pm2_5ptr%res*pm2_5ptr%err2
             end if

!  modify penalty term if nonlinear qc

             if (nlnqc_iter .and. pm2_5ptr%pg > tiny_r_kind .and. &
                  pm2_5ptr%b  > tiny_r_kind) then
                pm2_5_pg=pm2_5ptr%pg*varqc_iter
                cg_pm2_5=cg_term/pm2_5ptr%b
                wnotgross= one-pm2_5_pg
                wgross = pm2_5_pg*cg_pm2_5/wnotgross
                do kk=1,max(1,nstep)
                   pen(kk)= -two*log((exp(-half*pen(kk))+wgross)/(one+wgross))
                end do
             endif
             
             out(1) = out(1)+pen(1)*pm2_5ptr%raterr2
             do kk=2,nstep
                out(kk) = out(kk)+(pen(kk)-pen(1))*pm2_5ptr%raterr2
             end do
          end if
          
          !pm2_5ptr => pm2_5ptr%llpoint
          pm2_5ptr => pm2_5Node_nextcast(pm2_5ptr)

       end do

    endif

!!! cwy mosaic4 -------------------------------------------------------------- chem
    if( wrf_mass_regional .and. laeroana_mosaic4 ) then
      pm2_5ptr=>pm2_5Node_typecast(pm2_5head)
      
      do while( associated(pm2_5ptr) )
        if( pm2_5ptr%luse ) then
          if( nstep>0 ) then
            j1=pm2_5ptr%ij(1)
            j2=pm2_5ptr%ij(2)
            j3=pm2_5ptr%ij(3)
            j4=pm2_5ptr%ij(4)
            j5=pm2_5ptr%ij(5)
            j6=pm2_5ptr%ij(6)
            j7=pm2_5ptr%ij(7)
            j8=pm2_5ptr%ij(8)
            w1=pm2_5ptr%wij(1)
            w2=pm2_5ptr%wij(2)
            w3=pm2_5ptr%wij(3)
            w4=pm2_5ptr%wij(4)
            w5=pm2_5ptr%wij(5)
            w6=pm2_5ptr%wij(6)
            w7=pm2_5ptr%wij(7)
            w8=pm2_5ptr%wij(8)

            istatus=0
            val=zero
            val2=-pm2_5ptr%res    !!! log(ug/kg)

            if( .true. ) then

            aeroname='num_a01'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val =(w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8))*&
                  pm2_5ptr%dpmx_dnum(1)+val
            val2=(w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8))*&
                  pm2_5ptr%dpmx_dnum(1)+val2
            nullify(spm2_5,rpm2_5)

            aeroname='num_a02'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val =(w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8))*&
                  pm2_5ptr%dpmx_dnum(2)+val
            val2=(w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8))*&
                  pm2_5ptr%dpmx_dnum(2)+val2
            nullify(spm2_5,rpm2_5)

            aeroname='num_a03'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val =(w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8))*&
                  pm2_5ptr%dpmx_dnum(3)+val
            val2=(w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8))*&
                  pm2_5ptr%dpmx_dnum(3)+val2
            nullify(spm2_5,rpm2_5)

            else

            aeroname='so4_a01'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='so4_a02'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='so4_a03'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='no3_a01'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='no3_a02'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='no3_a03'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='nh4_a01'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='nh4_a02'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='nh4_a03'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='oc_a01'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='oc_a02'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='oc_a03'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='bc_a01'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='bc_a02'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='bc_a03'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='cl_a01'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='cl_a02'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='cl_a03'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='na_a01'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='na_a02'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='na_a03'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='oin_a01'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='oin_a02'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            aeroname='oin_a03'
            call gsi_bundlegetpointer(sval,trim(aeroname),spm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(456)
            end if
            call gsi_bundlegetpointer(rval,trim(aeroname),rpm2_5,istatus)
            if( istatus/=0 ) then
              write(6,*) 'error gsi_bundlegetpointer in stppm2_5 for '//aeroname
              call stop2(457)
            end if
            val = w1*rpm2_5(j1)+w2*rpm2_5(j2)+w3*rpm2_5(j3)+w4*rpm2_5(j4)+ &
                  w5*rpm2_5(j5)+w6*rpm2_5(j6)+w7*rpm2_5(j7)+w8*rpm2_5(j8)+ val
            val2= w1*spm2_5(j1)+w2*spm2_5(j2)+w3*spm2_5(j3)+w4*spm2_5(j4)+ &
                  w5*spm2_5(j5)+w6*spm2_5(j6)+w7*spm2_5(j7)+w8*spm2_5(j8)+ val2
            nullify(spm2_5,rpm2_5)

            end if

            do kk=1, nstep
              qq=val2+sges(kk)*val
              pen(kk)=qq*qq*pm2_5ptr%err2
            end do
          else
            pen(1)=pm2_5ptr%res*pm2_5ptr%res*pm2_5ptr%err2
          end if

          out(1)=out(1)+pen(1)*pm2_5ptr%raterr2
          do kk=2, nstep
            out(kk)=out(kk)+(pen(kk)-pen(1))*pm2_5ptr%raterr2
          end do
        end if

        pm2_5ptr=>pm2_5Node_nextcast(pm2_5ptr)
      end do
    end if
!!! cwy mosaic4 -------------------------------------------------------------- chem
       
    return
  end subroutine stppm2_5
  
end module stppm2_5mod
