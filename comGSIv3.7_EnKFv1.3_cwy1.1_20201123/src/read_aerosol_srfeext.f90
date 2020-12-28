subroutine read_aerosol_srfeext(nread,ndata,nodata,infile,obstype,gstime,lunout, &
           twind,sis,nobs)
!
!   output argument list:
!     nread    - number of modis aerosol observations read
!     ndata    - number of modis aerosol profiles retained for further processing
!     nodata   - number of modis aerosol observations retained for further processing
!     nobs     - array of observations on each subdomain for each processor
!
  use kinds,     only: r_kind, r_double, i_kind
  use gridmod,   only: nlat, nlon, regional, tll2xy, rlats, rlons
  use chemmod,   only: aod_qa_limit, luse_deepblue
  use constants, only: deg2rad, zero, two, three, four, r60inv
  use obsmod,    only: rmiss_single
  use gsi_4dvar, only: l4dvar,l4densvar,iwinbgn,winlen,thin4d
  use satthin,   only: itxmax,makegrids,destroygrids,checkob, &
      finalcheck,map2tgrid,score_crit
  use mpimod, only: npe
  use guess_grids, only: hrdifsig
  !!! cwy mosaic -----------------------------------------------
  use guess_grids, only: ges_prsi
  use obsmod, only: limit_srfeext
  use gridmod, only: nsig
  use convthin, only: make3grids
  !!! cwy mosaic -----------------------------------------------
  implicit none
!
! Declare local parameters
  real(r_kind), parameter :: r360 = 360.0_r_kind
!
! Declare passed variables
!
  character(len=*),intent(in)    :: obstype, infile
  character(len=20),intent(in)   :: sis
  integer(i_kind), intent(in)    :: lunout
  integer(i_kind), intent(inout) :: nread
  integer(i_kind),dimension(npe), intent(inout) :: nobs
  integer(i_kind), intent(inout) :: ndata, nodata
  real(r_kind),    intent(in)    :: gstime, twind
!
! Declare local variables
!
  logical :: outside, iuse
  
  character (len= 8) :: subset
  character (len=10) :: date

  integer(i_kind) :: naerodat
  integer(i_kind) :: idate, jdate, ksatid, iy, iret, im, ihh, idd
  integer(i_kind) :: lunin = 10
  integer(i_kind) :: nmind, i, n
  integer(i_kind) :: k, ilat, ilon, nreal, nchanl
  integer(i_kind) :: kidsat
  integer(i_kind), dimension(5) :: idate5
!
!!! cwy ----------------------------------------------- wavnum
  real(r_kind) :: sprs
  character(len=23) :: aerostr  = 'OPTD1 OPTD2 OPTD3'
!!! cwy ----------------------------------------------- wavnum
  character(len=46) :: aerogstr = 'SAID CLATH CLONH YEAR MNTH DAYS HOUR MINU PRSS'

  integer(i_kind) :: itx, itt, ksw

  real(r_kind) :: tdiff, sstime, dlon, dlat, t4dv, timedif, crit1, dist1
  real(r_kind) :: slons0, slats0, rsat, solzen, azimuth, dlat_earth, dlon_earth
  real(r_kind) :: dlat_earth_deg, dlon_earth_deg
  real(r_kind) :: styp, dbcf, qaod

  real(r_kind), allocatable, dimension(:,:) :: aeroout
  real(r_kind), allocatable, dimension(:)   :: dataaod
  real(r_double), dimension( 9) :: hdraerog
!!! cwy ----------------------------------------------- wavnum
  real(r_double), dimension( 3) :: srf_eext
!!! cwy ----------------------------------------------- wavnum

  integer, parameter                        :: ce318max=10000  !!! maximum ce318 obs

!**************************************************************************
! Set constants.  Initialize variables
  rsat=999._r_kind
  ! output position of LON and LAT
  ilon=3
  ilat=4
  nread = 0
  ndata = 0
  nodata = 0

  if ( obstype == 'srf_eext' ) then
     open(lunin,file=trim(infile),form='unformatted')
!!! cwy ------------------------------------------------- specified table
!!!     open(99,file=trim(infile)//'.table')
!!!     print *, 'cwy srf.eext.bufr.table=', trim(infile)//'.table'
!!!     call openbf(lunin,'IN',99)
!!! cwy ------------------------------------------------- specified table
     call openbf(lunin,'IN',lunin)
     call datelen(10)
     call readmg(lunin,subset,idate,iret)

     if ( iret == 0 ) then
        if (subset == 'NC008041') then
           write(6,*)'READ_AEROSOL_SRF_EEXT: SRF_EEXT data type, subset = ',subset
           nreal=6
           nchanl=20
           naerodat=nreal+nchanl
           allocate (aeroout(naerodat,ce318max))    !!! maximum ce318 obs
           allocate (dataaod(nchanl))
           aeroout=0.

           iy = 0
           im = 0
           idd= 0
           ihh= 0
           write(date,'( i10)') idate
           read (date,'(i4,3i2)') iy,im,idd,ihh
           write(6,'(''READ_AEROSOL_SRFEEXT: aerosol bufr file '',a,''  date is '',i4,4i2.2,a)')trim(infile),iy,im,idd,ihh

           read_loop: do
              call readsb(lunin,iret)
              if (iret/=0) then
                 call readmg(lunin,subset,jdate,iret)
                 if (iret/=0) exit read_loop
                 cycle read_loop
              endif
     
              call ufbint(lunin,hdraerog,9,1,iret,aerogstr)
              rsat = hdraerog(1); ksatid=rsat

              slats0= hdraerog(2)
              slons0= hdraerog(3)
              if(slons0< zero) slons0=slons0+r360
              if(slons0>=r360) slons0=slons0-r360
              dlat_earth_deg = slats0
              dlon_earth_deg = slons0
              dlat_earth = slats0 * deg2rad
              dlon_earth = slons0 * deg2rad

              call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
              if(outside) cycle read_loop

              !    Convert observation time to relative time
              idate5(1) = hdraerog(4)  !year
              idate5(2) = hdraerog(5)  !month
              idate5(3) = hdraerog(6)  !day
              idate5(4) = hdraerog(7)  !hour
              idate5(5) = hdraerog(8)  !minute

              call ufbint(lunin,srf_eext,3,1,iret,aerostr)

              call w3fs21(idate5,nmind)
              t4dv=real((nmind-iwinbgn),r_kind)*r60inv
              sstime=real(nmind,r_kind)
              tdiff=(sstime-gstime)*r60inv
              if (l4dvar.or.l4densvar) then
                 if(t4dv<zero .OR. t4dv>winlen) cycle read_loop
              else
                 !!! cwy mosaic --------------------- window for eext, default=0.5
                 if ( abs(tdiff)>twind .or. abs(tdiff)>0.5 ) then
                    cycle read_loop
                 end if
                 !!! cwy mosaic --------------------- window for eext, default=0.5
              end if

              if( limit_srfeext==1 ) then
              sprs=hdraerog(9)/100.  ! sprs, hPa
              if( sprs<100000. .and. sprs>ges_prsi(dlat,dlon,2,1)*10. ) then
              else
                cycle read_loop
              end if
              end if

              ! if any wavelength aod is missing, cycle
!!! cwy ----------------------------------------------- wavnum
              do ksw=2, 2
!!! cwy ----------------------------------------------- wavnum
              if ( srf_eext(ksw)>=10000_r_double .or. srf_eext(ksw)<=0. ) cycle read_loop
              end do

              nread = nread + 1   !nread = nread + nchanl

              dataaod = rmiss_single
!!! cwy ----------------------------------------------- wavnum
!!! NOTE: 1-3 for absorbing aod, 4-7 for scattering aod
              dataaod(1) = srf_eext(1)
              dataaod(2) = srf_eext(2)
              dataaod(3) = srf_eext(3)
              print '(a17,3f10.3)', 'cwy obs srf_eext=', dataaod(1), dataaod(2), dataaod(3)
!!! cwy ----------------------------------------------- wavnum

              itx=nread
              if( itx>ce318max ) then
                write(*,*) 'cwy stop error: ce318obs>ce318max in read_aerosol_srfeext.f90', itx
                stop
              end if

              aeroout( 1,itx) = rsat
              aeroout( 2,itx) = tdiff
              aeroout( 3,itx) = dlon               ! grid relative longitude
              aeroout( 4,itx) = dlat               ! grid relative latitude
              aeroout( 5,itx) = dlon_earth_deg     ! earth relative longitude (degrees)
              aeroout( 6,itx) = dlat_earth_deg     ! earth relative latitude (degrees)
              do i = 1, nchanl
                 aeroout(i+nreal,itx) = dataaod(i)
              end do

              ndata=ndata+1
       
           end do read_loop

           do n = 1, ndata
              do i = 1, nchanl
                 if ( aeroout(i+nreal,n) > rmiss_single ) nodata = nodata + 1
              end do
           end do
           call count_obs(ndata,naerodat,ilat,ilon,aeroout,nobs)
           write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
           write(lunout) ((aeroout(k,n),k=1,naerodat),n=1,ndata)

           ! Deallocate local arrays
           deallocate(aeroout)
           deallocate(dataaod)

        else       ! subset /= NC008041
           write(6,*)'READ_AEROSOL_SRFEEXT:  *** WARNING: unknown aerosol data type, subset=',subset
           write(6,*)' infile=',infile, ', lunin=',lunin, ', obstype=',obstype
           write(6,*)' SKIP PROCESSING OF THIS MODIS FILE'
        endif

     else          ! read subset iret /= 0
        write(6,*)'READ_AEROSOL_SRFEEXT:  *** WARNING: read subset error, obstype=',obstype,', iret=',iret
     end if
     call closbf(lunin)
     close(lunin)
  else             ! obstype /= 'modis'
     write(6,*)'READ_AEROSOL_SRFEEXT:  *** WARNING: unknown aerosol input type, obstype=',obstype
  endif

  ! Deallocate satthin arrays
  call destroygrids

end subroutine read_aerosol_srfeext
