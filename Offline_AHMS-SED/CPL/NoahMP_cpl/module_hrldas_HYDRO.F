module module_HRLDAS_HYDRO

! NDHMS  module   
#ifdef MPP_LAND
    use module_mpp_land, only: global_nx, global_ny, decompose_data_real, &
                               write_io_real, my_id, mpp_land_bcast_real1, IO_id, &
                               mpp_land_bcast_int1
#endif
    use module_HYDRO_drv, only: HYDRO_ini, HYDRO_exe, HYDRO_rst_out
    use module_rt_data,   only: rt_domain
    use module_namelist,  only: nlst_rt
    use module_gw_gw2d_data, only: gw2d
    use module_rt_rt2d_data, only: rt2d  ! Qian

    implicit none
    integer :: begg, endg
    integer :: numg, numl, numc, nump
    INTEGER, PARAMETER :: double=8
    real(kind=double), pointer :: r2p(:,:) , r1p(:)
    integer ::  begl, endl, begc, endc, begp, endp
    integer :: nn
    integer :: open_unit_status
    real, allocatable, dimension(:,:) :: vg_test

CONTAINS

    subroutine hrldas_cpl_HYDRO(STC,     SMC,       SH2OX,    &
                                infxsrt, sfcheadrt, soldrain, &
                                solinfl, hgwrt,               &
                                ht,      prcp,      evap,     &
                                ii,      jj,        kk,       &
                                xland,   IRRIG,     K1_FACTOR,& 
                                slope,   GINDEX,    LAI,      &
                                IRFRACT )
       implicit none

       integer,intent(in) :: ii     ! grid number in west-west direction
       integer,intent(in) :: jj     ! grid number in south-north direction
       integer,intent(in) :: kk     ! number of soil layers

       real,   intent(inout) :: STC(ii,jj,kk)
       real,   intent(inout) :: SMC(ii,jj,kk)
       real,   intent(inout) :: SH2OX(ii,jj,kk)
       real,   intent(inout) :: infxsrt(ii,jj)
       real,   intent(inout) :: sfcheadrt(ii,jj)
       real,   intent(inout) :: soldrain(ii,jj)
       real,   intent(inout) :: solinfl(ii,jj)
       real,   intent(inout) :: hgwrt(ii,jj)! groundwater head, m
       real,   intent(in)    :: ht(ii,jj)   ! terrain height, m
       real,   intent(in)    :: prcp(ii,jj) ! precipitation, mm
       real,   intent(in)    :: evap(ii,jj) ! evapotranspiration, mm
       real,   intent(inout) :: IRRIG(ii,jj) ! irrigation, mm
       real,   intent(in)    :: K1_FACTOR(ii,jj) !
       real,   intent(in)    :: slope(ii,jj) !
       real,   intent(in)    :: GINDEX(ii,jj) !
       real,   intent(in)    :: LAI(ii,jj) !
       real,   intent(in)    :: IRFRACT(ii,jj) !
       real,   intent(inout) :: xland(ii,jj)! index of land (=1) and water (=2)

!lpr add 2014-06-24
#ifdef WRF_HYDRO_RAPID
#ifdef MPP_LAND
        real, dimension(ii,jj) :: g_runoff
#endif
        real, dimension(ii,jj) :: runoff
#endif
!lpr add end
        integer :: did
        integer :: ntime
        integer :: i, j, k
        real*8  :: t1, t2, dact
        save dact

        ! output flux and state variable
        did = 1
        ntime = 1
!--------LPR add 2014-06-24---------------------------
!    it is rapid model
#ifdef WRF_HYDRO_RAPID
        if(nlst_rt(did)%channel_option .eq. 4) then
                runoff  = infxsrt + soldrain
                call hrldas_RAPID_ini(ntime)    
#ifdef MPP_LAND
                call write_io_real(runoff,g_runoff)
                if(my_id.eq.io_id) &
                call hrldas_RAPID_exe(runoff,g_global_nx,global_ny)     
#else
                call hrldas_RAPID_exe(runoff,ii,jj)     
#endif
                sfcheadrt = 0.0
                return
        endif
#endif
!--------LPR add end----------------------------------


        if(  nlst_rt(did)%SUBRTSWCRT  .eq. 0   .and. &
             nlst_rt(did)%OVRTSWCRT   .eq. 0   .and. &
             nlst_rt(did)%GWBASESWCRT .eq. 0   .and. &
             nlst_rt(did)%OVCHRTSWCRT .eq. 0   .and. &   ! Qian, 2016.08
             nlst_rt(did)%INTERGWSFW  .eq. 0         &
          )  return

        ! decompose the hrldas 1-d data into routing domain
        do k = 1, kk 
           RT_DOMAIN(did)%STC(:,:,k)   = STC(1:ii,1:jj,k)
           RT_DOMAIN(did)%SMC(:,:,k)   = SMC(1:ii,1:jj,k)
           RT_DOMAIN(did)%SH2OX(:,:,k) = SH2OX(1:ii,1:jj,k)
        end do
        RT_DOMAIN(did)%infxsrt  = infxsrt(1:ii,1:jj)  ! surface runoff, mm
        RT_DOMAIN(did)%soldrain = soldrain(1:ii,1:jj) ! unsaturated soil drainage, mm
        rt_domain(did)%SOLINFL  = solinfl(1:ii,1:jj)  ! net flux crossing ground surface, mm
        rt_domain(did)%PRECIPRT = prcp(1:ii,1:jj)     ! timestep precipitation, mm
        rt_domain(did)%EVAPRT   = evap(1:ii,1:jj)     ! timetep evapotranspiration, mm
        rt_domain(did)%IRRIG    = IRRIG(1:ii,1:jj)    ! timestep irrigation, mm cjiang,water use, use factorfrom papera jia, 2006
        rt_domain(did)%K1_FACTOR= K1_FACTOR(1:ii,1:jj)!
        rt_domain(did)%slope    = slope(1:ii,1:jj)    !
        rt_domain(did)%GINDEX   = GINDEX(1:ii,1:jj)   !
        rt_domain(did)%LAI      = LAI(1:ii,1:jj)      !
        rt_domain(did)%IRFRACT  = IRFRACT(1:ii,1:jj)  !
#ifdef MPP_LAND
        if(my_id .eq. IO_id) then
           call  time_seconds(t1)
        endif
#endif

        ! step hydrological routing
        call HYDRO_exe(did)

#ifdef MPP_LAND
        if(my_id .eq. IO_id) then
           call  time_seconds(t2)
           dact = dact + t2 - t1
           write(6,*) "accumulated time (s): ",dact
        endif
#endif

           ! Qian, 2018.08
           ! feedback land-lake mask for LSM, by XLAND
!           IF (nlst_rt(did)%INTERGWSFW .eq. 1) THEN
!              do i = 1, ii
!                 do j = 1, jj
!                    if (rt2d(did)%hend(i,j) - rt2d(did)%elev(i,j) .ge. -1.E-3) then
!                       xland(i,j) = 2.0
!                    endif
!                 end do
!              end do
!           ENDIF

        ! add for update the HRLDAS state variable.
        do k = 1, kk
           STC(1:ii,1:jj,k)   = rt_domain(did)%STC(:,:,k)
           SMC(1:ii,1:jj,k)   = rt_domain(did)%SMC(:,:,k)
           SH2OX(1:ii,1:jj,k) = rt_domain(did)%SH2OX(:,:,k)
        end do
        sfcheadrt(1:ii,1:jj)  = rt_domain(did)%sfcheadrt

       ! feedback groundwater head for land surface model
       IF (nlst_rt(did)%INTERGWSFW .eq. 1) THEN
          hgwrt(1:ii,1:jj) = ht(1:ii,1:jj) &         ! transfer depth of groundwater head
                           - (rt2d(did)%elev - rt2d(did)%hend)
       ELSEIF (nlst_rt(did)%GWBASESWCRT .eq. 3 ) THEN
          hgwrt(1:ii,1:jj) = ht(1:ii,1:jj) &         ! transfer depth of groundwater head
                           - (gw2d(did)%elev - gw2d(did)%h)
       ENDIF

      ! constrained by available river water

       IRRIG(1:ii,1:jj) = rt_domain(did)%IRRIG(:,:)

       print*, rt_domain(did)%IRRIG(85,4), "rt_domain(did)%IRRIG(85,4)"

!? not sure for the following
!           grid%xice(its:ite,jts:jte) = rt_domain(did)%sice

     end subroutine hrldas_cpl_HYDRO

    subroutine hrldas_cpl_HYDRO_ini(STC,SMC,SH2OX,infxsrt,sfcheadrt,soldrain,hgwrt,&
                                    ht,ii,jj,kk,dt, olddate,zsoil,isltyp,ivgtyp,xland)
       implicit none

       integer,intent(in) :: ii     ! grid number in west-west direction
       integer,intent(in) :: jj     ! grid number in south-north direction
       integer,intent(in) :: kk     ! number of soil layers

       real,   intent(in) :: dt     ! NOAH_TIMESTEP, s
       real,   intent(inout) :: STC(ii,jj,kk)
       real,   intent(inout) :: SMC(ii,jj,kk)
       real,   intent(inout) :: SH2OX(ii,jj,kk)
       real,   intent(inout) :: infxsrt(ii,jj)
       real,   intent(inout) :: sfcheadrt(ii,jj)
       real,   intent(inout) :: soldrain(ii,jj)
       real,   intent(inout) :: hgwrt(ii,jj)
       real,   intent(in)    :: ht(ii,jj)
       real,   intent(in)    :: zsoil(kk)
       integer,intent(inout) :: isltyp(ii,jj)
       integer,intent(inout) :: ivgtyp(ii,jj)
       real,   intent(inout) :: xland(ii,jj)  ! Qian added
       character(len = *)    :: olddate
     
       ! LOCAL VARIABLES
       integer :: did
       integer :: ntime
       integer :: i,j,k

       ! output flux and state variable
        did = 1

        if(.not. RT_DOMAIN(did)%initialized) then  
           nlst_rt(did)%dt = dt
           nlst_rt(did)%olddate(1:19)   = olddate(1:19)
           nlst_rt(did)%startdate(1:19) = olddate(1:19)

           nlst_rt(did)%nsoil = kk
#ifdef MPP_LAND
           call mpp_land_bcast_int1(nlst_rt(did)%nsoil)
#endif
           allocate(nlst_rt(did)%zsoil8(nlst_rt(did)%nsoil))
           allocate(nlst_rt(did)%dzsoil(nlst_rt(did)%nsoil))

           ! Qian
           if(zsoil(1) <  0.) then
              nlst_rt(did)%zsoil8(1:nlst_rt(did)%nsoil) = zsoil(1:nlst_rt(did)%nsoil)
           else
              nlst_rt(did)%zsoil8(1:nlst_rt(did)%nsoil) = -1.*zsoil(1:nlst_rt(did)%nsoil)
           endif

           nlst_rt(did)%dzsoil(1) = -1.*nlst_rt(did)%zsoil8(1)
           do k = 2, kk
              nlst_rt(did)%dzsoil(k) = -1.*(nlst_rt(did)%zsoil8(k) - nlst_rt(did)%zsoil8(k-1)) 
           end do           

!           call HYDRO_ini(ntime,did,ix0=1,jx0=1)
           call HYDRO_ini(ntime,did,ix0=ii,jx0=jj,vegtyp=ivgtyp, soltyp=isltyp)  ! Qian

           if (nlst_rt(did)%sys_cpl .ne. 1) then
              write(6,*) "Error: sys_cpl should be 1."
              call hydro_stop()
           endif

           RT_DOMAIN(did)%SLDPTH(1:nlst_rt(did)%nsoil) = nlst_rt(did)%dzsoil(1:nlst_rt(did)%nsoil)

           RT_DOMAIN(did)%initialized = .true.

#ifdef WRF_HYDRO_RAPID
!--------LPR add 2014-06-24---------------------------
         if(nlst_rt(did)%channel_option .eq. 4) then
!  it  is rapid model.
             return
         endif
!--------LPR add 2014-06-24---------------------------
#endif

        write(6,*) "Time Interval Determination in hrldas_cpl_HYDRO"
        write(6,*) "Noah  interval:   nlst_rt(did)%dt = ", nlst_rt(did)%dt
        write(6,*) "HYDRO interval: nlst_rt(did)%dtrt = ", nlst_rt(did)%dtrt
        write(6,*) "2D Terrain Routing interval: nlst_rt(did)%DT_2DTHRT = ", nlst_rt(did)%DT_2DTHRT

        write(6,*) "Soil Layers Determination in hrldas_cpl_HYDRO"
        write(6,*) "nlst_rt(did)%nsoil  = ", nlst_rt(did)%nsoil
        write(6,*) "nlst_rt(did)%dzsoil = ", nlst_rt(did)%dzsoil
        write(6,*) "nlst_rt(did)%zsoil8 = ", nlst_rt(did)%zsoil8

        if(  nlst_rt(did)%SUBRTSWCRT  .eq. 0   .and. &
             nlst_rt(did)%OVRTSWCRT   .eq. 0   .and. &
             nlst_rt(did)%GWBASESWCRT .eq. 0   .and. &
             nlst_rt(did)%OVCHRTSWCRT .eq. 0   .and. &   ! Qian, 2016.08
             nlst_rt(did)%INTERGWSFW  .eq. 0         &
          )  return

#ifdef MPP_LAND
           call mpp_land_bcast_real1(nlst_rt(did)%dt)
#endif

           sfcheadrt(1:ii,1:jj) = &
           rt_domain(did)%sfcheadrt(1:rt_domain(did)%ix,1:rt_domain(did)%jx)

           ! feedback groundwater head for land surface model
           IF (nlst_rt(did)%INTERGWSFW .eq. 1) THEN
              hgwrt(1:ii,1:jj) = ht(1:ii,1:jj) &         
                                 - (rt2d(did)%elev - rt2d(did)%hend)
           ELSEIF (nlst_rt(did)%GWBASESWCRT .eq. 3 ) THEN
              hgwrt(1:ii,1:jj) = ht(1:ii,1:jj) &         
                                 - (gw2d(did)%elev - gw2d(did)%h)
           ENDIF

           ! Qian, 2018.08
           ! feedback lake mask for LSM, by XLAND
!           IF (nlst_rt(did)%INTERGWSFW .eq. 1) THEN
!              do i = 1, ii 
!                 do j = 1, jj
!                    if (rt2d(did)%hend(i,j) - rt2d(did)%elev(i,j) .ge. -1.E-3) then
!                       xland(i,j) = 2.0
!                    endif
!                 end do
!              end do
!           ENDIF

           ! The soil temperature is not updated in hydrological model now 
!           STC(1:ii,1:jj,1:kk)  = rt_domain(did)%STC(1:rt_domain(did)%ix,1:rt_domain(did)%jx,1:nlst_rt(did)%nsoil)
           SMC(1:ii,1:jj,1:kk)  = rt_domain(did)%SMC(1:rt_domain(did)%ix,1:rt_domain(did)%jx,1:nlst_rt(did)%nsoil)
           SH2OX(1:ii,1:jj,1:kk)= rt_domain(did)%SH2OX(1:rt_domain(did)%ix,1:rt_domain(did)%jx,1:nlst_rt(did)%nsoil)

!           if(nlst_rt(did)%rst_typ .eq. 1) then  ! restart
!              STC   = rt_domain(did)%STC
!              SMC   = rt_domain(did)%SMC
!              SH2OX = rt_domain(did)%SH2OX
!           else
!              if(nlst_rt(did)%sys_cpl .eq. 1) then
!                 where( abs(STC)   .gt. 500) stc = 282
!                 where( abs(SMC)   .gt. 500) SMC = 0.25
!                 where( abs(SH2OX) .gt. 500) SH2OX = 0.25
!              endif
!           endif

        endif
     end subroutine hrldas_cpl_HYDRO_ini

     subroutine open_print_mpp(iunit)
       implicit none
       integer iunit
       character(len=16) fileout

       if(open_unit_status == 999) return
         open_unit_status = 999 

#ifdef MPP_LAND  
       write(fileout,'(a11,i0.5)') "diag_hydro.",my_id    
#else
       write(fileout,'(a11,i0.5)') "diag_hydro.",0    
#endif
       open(iunit,file=fileout,form="formatted")

     end subroutine open_print_mpp

end module module_HRLDAS_HYDRO
