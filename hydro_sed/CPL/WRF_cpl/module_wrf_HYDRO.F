module module_WRF_HYDRO

#ifdef MPP_LAND
    use module_mpp_land, only: global_nx, global_ny, decompose_data_real, &
                               write_io_real, my_id, mpp_land_bcast_real1, &
                               IO_id, mpp_land_bcast_real, mpp_land_bcast_int1
    use module_CPL_LAND, only: CPL_LAND_INIT, cpl_outdate
#endif
    use module_HYDRO_drv,    only: HYDRO_ini, HYDRO_exe
    use module_rt_data,      only: rt_domain
    use module_gw_gw2d_data, only: gw2d
    use module_rt_rt2d_data, only: rt2d   ! Qian, 2016.09
!   use module_CPL_LAND,     only: CPL_LAND_INIT, cpl_outdate !Qian, 2016.05.04
    use module_namelist,     only: nlst_rt
    USE module_domain,       only: domain, domain_clock_get
    USE module_configure,    only: grid_config_rec_type
!    USE module_configure,    only: config_flags
    USE module_configure,    only: model_config_rec

    implicit none
     
    !yw  added for check soil moisture and soiltype
    integer ::  checkSOIL_flag

#ifndef MPP_LAND
    character(len=19) :: cpl_outdate
#endif

    ! added to consider the adaptive time step from WRF model. 
    real    :: dtrt0  
    integer :: mm0
    real, save :: time_sim   ! simulated time from beginning time

CONTAINS

! -----------------------------------------------
! wrf_cpl_HYDRO is used to couple WRF with HYDRO
! wrf_cpl_HYDRO will not call the off-line lsm 
! -----------------------------------------------
 
!   subroutine wrf_cpl_HYDRO(HYDRO_dt,grid, config_flags, its,ite,jts,jte) ! yw
    subroutine wrf_cpl_HYDRO(HYDRO_dt, grid, its, ite, jts, jte)

       implicit none
       TYPE ( domain ), INTENT(INOUT) :: grid
!       TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags ! ywGW
       real,    intent(inout) :: HYDRO_dt
       integer, intent(in)    :: its
       integer, intent(in)    :: ite
       integer, intent(in)    :: jts
       integer, intent(in)    :: jte

       integer :: k, ix, jx, mm, nn
       integer :: did
       integer :: ntime
       integer :: i, j

       ! output flux and state variable
        did   = 1
        ix    = ite - its + 1
        jx    = jte - jts + 1
        ntime = 1

        if (HYDRO_dt .le. 0) then
           write(6,*) "Warning: HYDRO_dt <= 0 from land input. set it to be 1 seconds."
           HYDRO_dt = 1.
        endif

        nlst_rt(did)%dt = HYDRO_dt 
  
ifini:  if (.not. RT_DOMAIN(did)%initialized) then

           !nlst_rt(did)%nsoil = config_flags%num_soil_layers     ! yw
           !nlst_rt(did)%nsoil = model_config_rec%num_metgrid_soil_levels
           nlst_rt(did)%nsoil  = grid%num_soil_layers
         
#ifdef MPP_LAND
           call mpp_land_bcast_int1(nlst_rt(did)%nsoil)
#endif
           allocate(nlst_rt(did)%zsoil8(nlst_rt(did)%nsoil))
           allocate(nlst_rt(did)%dzsoil(nlst_rt(did)%nsoil))  ! Qian, 2016.12
           nlst_rt(did)%dzsoil(1:nlst_rt(did)%nsoil) = grid%DZS(1:nlst_rt(did)%nsoil)
           if(grid%zs(1) <  0.) then
              nlst_rt(did)%zsoil8(1:nlst_rt(did)%nsoil) = grid%zs(1:nlst_rt(did)%nsoil)
           else
              nlst_rt(did)%zsoil8(1:nlst_rt(did)%nsoil) = -1.*grid%zs(1:nlst_rt(did)%nsoil)
           endif

           CALL domain_clock_get(grid, current_timestr=cpl_outdate)
           nlst_rt(did)%startdate(1:19) = cpl_outdate(1:19)
           nlst_rt(did)%olddate(1:19)   = cpl_outdate(1:19)

#ifdef MPP_LAND
           call CPL_LAND_INIT(its,ite,jts,jte)
#endif

#ifdef HYDRO_D
           write(6,*) "---------------------------------------"
           write(6,*) "Start HYDRO_ini ......" 
           write(6,*) "sf_surface_physics is ", grid%sf_surface_physics
#endif

           if(grid%sf_surface_physics .eq. 5) then    
               ! coupled with CLM4 
               call HYDRO_ini(ntime, did=did, ix0=1,  jx0=1)
           else
               call HYDRO_ini(ntime, did,     ix0=ix, jx0=jx,      &
                              vegtyp=grid%IVGTYP(its:ite,jts:jte), & 
                              soltyp=grid%isltyp(its:ite,jts:jte)  )
           endif

           RT_DOMAIN(did)%SLDPTH(1:nlst_rt(did)%nsoil) = nlst_rt(did)%dzsoil(1:nlst_rt(did)%nsoil)
#ifdef HYDRO_D
           write(6,*) "Finish HYDRO initilization"
           write(6,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
#endif

           if(nlst_rt(did)%sys_cpl .ne. 2) then
               write(6,*) "Error: coupled with WRF "
               write(6,*) "Error: sys_cpl should be 2 !"
               call hydro_stop()
           endif

           nlst_rt(did)%startdate(1:19) = cpl_outdate(1:19)
           nlst_rt(did)%olddate(1:19)   = cpl_outdate(1:19)

           ! Hydro timestep setting
           nlst_rt(did)%dt = HYDRO_dt
           dtrt0 = nlst_rt(did)%dtrt 
           if(nlst_rt(did)%dtrt .ge. HYDRO_dt) then
              !qxia nlst_rt(did)%dtrt = HYDRO_dt
              !qxia mm0 = 1

              mm = dtrt0/HYDRO_dt
              if (mm*HYDRO_dt .lt. dtrt0) nlst_rt(did)%dtrt = mm*HYDRO_dt
              mm0 = mm
           else
              mm = HYDRO_dt/nlst_rt(did)%dtrt
              if(mm*nlst_rt(did)%dtrt .lt. HYDRO_dt) nlst_rt(did)%dtrt = HYDRO_dt/mm
              mm0 = mm
           endif
           dtrt0 = nlst_rt(did)%dtrt 

           write(6,*) "Timestep Interval Determination in wrf_cpl_HYDRO"
           write(6,*) "WRF interval:   nlst_rt(did)%dt   = ", nlst_rt(did)%dt
           write(6,*) "HYDRO interval: nlst_rt(did)%dtrt = ", nlst_rt(did)%dtrt
           write(6,*) "2D Terrain Routing interval: nlst_rt(did)%DT_2DTHRT = ", nlst_rt(did)%DT_2DTHRT
           write(6,*) "Ratio of HYDRO interval to WRF interval = ", mm0

           ! initialize, qxia 
           time_sim = 0.
           rt_domain(did)%soldrain = 0.
           rt_domain(did)%infxsrt  = 0.
           rt_domain(did)%SOLINFL  = 0.
           rt_domain(did)%PRECIPRT = 0.
           rt_domain(did)%EVAPRT   = 0.  
 
        endif ifini

        ! if WRF model time step changed
        if((mm0*nlst_rt(did)%dtrt) .ne. HYDRO_dt) then 
           if (dtrt0 .ge. HYDRO_dt) then
              !qxia nlst_rt(did)%dtrt = HYDRO_dt
              !qxia mm0 = 1
            
              mm = dtrt0/HYDRO_dt
              if (mm*HYDRO_dt .lt. dtrt0) nlst_rt(did)%dtrt = mm*HYDRO_dt
              mm0 = mm  
           else
              mm = HYDRO_dt/dtrt0
              if(mm*dtrt0 .lt. HYDRO_dt) nlst_rt(did)%dtrt = HYDRO_dt/mm
              mm0 = mm
           endif
        endif

        if(  nlst_rt(did)%SUBRTSWCRT  .eq. 0   .and. &
             nlst_rt(did)%OVRTSWCRT   .eq. 0   .and. &
             nlst_rt(did)%GWBASESWCRT .eq. 0   .and. &
             nlst_rt(did)%OVCHRTSWCRT .eq. 0   .and. &   ! Qian, 2016.08
             nlst_rt(did)%INTERGWSFW  .eq. 0         &
          )  return

        ! simulated time in WRF
        time_sim = time_sim + nlst_rt(did)%dt
        ! number of soil layers  
        nn = nlst_rt(did)%nsoil

        ! get the data from WRF 
       if( (.not. RT_DOMAIN(did)%initialized) .and. (nlst_rt(did)%rst_typ .eq. 1) ) then
#ifdef HYDRO_D
          write(6,*) "restart initial data from offline file"
#endif
       else
          do k = 1, nlst_rt(did)%nsoil
             RT_DOMAIN(did)%STC(:,:,k)   = grid%TSLB(its:ite,k,jts:jte) 
             RT_DOMAIN(did)%smc(:,:,k)   = grid%smois(its:ite,k,jts:jte) 
             RT_DOMAIN(did)%sh2ox(:,:,k) = grid%sh2o(its:ite,k,jts:jte) 
          end do
          rt_domain(did)%sfcheadrt = grid%sfcheadrt(its:ite,jts:jte)                             ! surface runoff head at current time step, mm
          ! accumulated value at every hydrological timestep
          rt_domain(did)%soldrain  = rt_domain(did)%soldrain  + grid%soldrain(its:ite,jts:jte)   ! time-step soil drainage to groundwater, mm
          rt_domain(did)%infxsrt   = rt_domain(did)%infxsrt   + grid%infxsrt(its:ite,jts:jte)    ! time-step surface runoff, mm 
          rt_domain(did)%SOLINFL   = rt_domain(did)%SOLINFL   + grid%solinfl(its:ite,jts:jte)    ! time-step net infiltration, mm
          rt_domain(did)%PRECIPRT  = rt_domain(did)%PRECIPRT  + &                                ! time-step precipitation, mm
                                     grid%rainshv(its:ite,jts:jte) + & 
                                     grid%rainncv(its:ite,jts:jte) + &
                                     grid%raincv(its:ite,jts:jte) 
          rt_domain(did)%EVAPRT    = rt_domain(did)%EVAPRT + grid%qfx(its:ite,jts:jte)*HYDRO_dt  ! time-step evapotranspiration, mm
       endif 

      ! running HYDRO at hydrological timestep
      if (mod(time_sim,nlst_rt(did)%dtrt) .eq. 0.) then
#ifdef HYDRO_D
         write(6,*) "---------------------------------------"
         write(6,*) "Start HYDRO_exe ......"
#endif
         ! step hydrological routing
         call HYDRO_exe(did)
#ifdef HYDRO_D
         write(6,*) "Finish HYDRO_exe"
         write(6,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
#endif

         ! updated variables in HYDRO to WRF
         ! soil layers 
         ! soil moisture is updated in HYDRO if nlst_rt(did)%INTERGWSFW == 1
         IF (nlst_rt(did)%INTERGWSFW .eq. 1) THEN  
            do k = 1, nlst_rt(did)%nsoil
               ! grid%TSLB(its:ite,k,jts:jte) = RT_DOMAIN(did)%STC(:,:,k)
               grid%smois(its:ite,k,jts:jte)  = RT_DOMAIN(did)%smc(:,:,k) 
               grid%sh2o(its:ite,k,jts:jte)   = RT_DOMAIN(did)%sh2ox(:,:,k)
            end do
         ENDIF

         ! terrain water condition 
         ! surface water head is updated in HYDRO if nlst_rt(did)%OVCHRTSWCRT .eq. 1
         IF (nlst_rt(did)%OVCHRTSWCRT .eq. 1) THEN
            grid%sfcheadrt(its:ite,jts:jte) = rt_domain(did)%sfcheadrt
         ENDIF

         ! groundwater head
         ! currently only CLMS & NoahMP 
         IF ( (grid%sf_surface_physics .eq. 5) .or. & ! CLM4   land surface model 
            (grid%sf_surface_physics .eq. 4)      & ! NoahMP land surface model
            ) THEN

            IF (nlst_rt(did)%INTERGWSFW .eq. 1) THEN
               grid%hgwrt(its:ite,jts:jte) = grid%ht(its:ite,jts:jte) &         ! transfer depth of groundwater head
                                           - (rt2d(did)%elev - rt2d(did)%hend)
            ELSEIF (nlst_rt(did)%GWBASESWCRT .eq. 3 ) THEN
               grid%hgwrt(its:ite,jts:jte) = grid%ht(its:ite,jts:jte) &         ! transfer depth of groundwater head
                                           - (gw2d(did)%elev - gw2d(did)%h)
            ENDIF

         ENDIF

         ! zero accumulated variable 
         rt_domain(did)%soldrain = 0.  ! accumulated soil drainage (mm) during hdrological timestep, DTRT
         rt_domain(did)%infxsrt  = 0.
         rt_domain(did)%SOLINFL  = 0.
         rt_domain(did)%PRECIPRT = 0.
         rt_domain(did)%EVAPRT   = 0.

      endif 
          !print*,'QIAN, HYDRO ---> WRF: '
          !print*, 'GRID%SMC(25,25)  = ', grid%smois(25,:,25)
          !print*, 'GRID%SH2O(25,25) = ', grid%sh2o(25,:,25)
          !print*, 'GRID%sfcheadrt   = ', grid%sfcheadrt(25,25)

       ! provide groundwater soil flux to WRF for fully coupled simulations (FERSCH 09/2014)            
!       if(nlst_rt(did)%GWBASESWCRT .eq. 3 ) then
          ! yw, comment the following two lines. Not ready for WRF3.7 release
!          grid%qsgw(its:ite,jts:jte) = gw2d(did)%qsgw           ! yw
!          config_flags%gwsoilcpl     = nlst_rt(did)%gwsoilcpl   ! yw
!       end if

       ! yw, not sure for the following
       ! grid%xice(its:ite,jts:jte) = rt_domain(did)%sice

       RT_DOMAIN(did)%initialized = .true.

#ifdef HYDRO_D
           write(6,*) "###########################################################"
#endif

     end subroutine wrf_cpl_HYDRO

!program drive rtland
! This subroutine will be used if the 4-layer Noah lsm is not used.
      subroutine wrf2lsm (z1,v1,kk1,z,vout,ix,jx,kk,vegtyp)
!  input: z1,v1,kk1,z,ix,jx,kk
!  output: vout
!  interpolate based on soil layer: z1 and z 
!  z :  soil layer of output variable.
!  z1: array of soil layers of input variable.
         implicit none
         integer:: i,j,k
         integer:: kk1, ix,jx,kk, vegtyp(ix,jx)
         real :: z1(kk1), z(kk), v1(ix,kk1,jx),vout(ix,jx,kk)

       
         do j = 1, jx
            do i = 1, ix
                do k = 1, kk
                  call interpLayer(Z1,v1(i,1:kk1,j),kk1,Z(k),vout(i,j,k)) 
                end do
            end do
         end do
      end subroutine wrf2lsm

! This subroutine will be used if the 4-layer Noah lsm is not used.
      subroutine lsm2wrf (z1,v1,kk1,z,vout,ix,jx,kk,vegtyp)
!  input: z1,v1,kk1,z,ix,jx,kk
!  output: vout
!  interpolate based on soil layer: z1 and z 
!  z :  soil layer of output variable.
!  z1: array of soil layers of input variable.
         implicit none
         integer:: i,j,k
         integer:: kk1, ix,jx,kk, vegtyp(ix,jx)
         real :: z1(kk1), z(kk), v1(ix,jx,kk1),vout(ix,kk,jx)

       
         do j = 1, jx
            do i = 1, ix
                 do k = 1, kk
                    call interpLayer(Z1,v1(i,j,1:kk1),kk1,Z(k),vout(i,k,j)) 
                 end do
            end do
         end do
      end subroutine lsm2wrf

      subroutine interpLayer(inZ,inV,inK,outZ,outV)
         implicit none
         integer:: k, k1, k2
         integer :: inK
         real:: inV(inK),inZ(inK)
         real:: outV, outZ, w1, w2

         if(outZ .le. inZ(1)) then
             w1 = (inZ(2)-outZ)/(inZ(2)-inZ(1))
             w2 = (inZ(1)-outZ)/(inZ(2)-inZ(1))
             outV = inV(1)*w1-inV(2)*w2
             return
         elseif(outZ .ge. inZ(inK)) then
             w1 = (outZ-inZ(inK-1))/(inZ(inK)-inZ(inK-1)) 
             w2 = (outZ-inZ(inK))  /(inZ(inK)-inZ(inK-1))
             outV = inV(inK)*w1 -inV(inK-1)* w2
             return
         else  
            do k = 2, inK
             if((inZ(k) .ge. outZ).and.(inZ(k-1) .le. outZ) ) then
                k1  = k-1
                k2 = k
                w1 = (outZ-inZ(k1))/(inZ(k2)-inZ(k1))
                w2 = (inZ(k2)-outZ)/(inZ(k2)-inZ(k1))
                outV = inV(k2)*w1 + inV(k1)*w2
                return 
             end if 
            end do
         endif
      end subroutine interpLayer

      subroutine lsm_wrf_input(did,vegtyp,soltyp,ix,jx)
         implicit none
         integer did, leng
         parameter(leng=100)
         integer :: i,j, nn, ix,jx
         integer, dimension(ix,jx) :: soltyp, vegtyp
         real, dimension(leng) :: xdum1, MAXSMC,refsmc,wltsmc


         where(soltyp == 14) VEGTYP = 16
         where(VEGTYP == 16 ) soltyp = 14

         RT_DOMAIN(did)%VEGTYP = vegtyp

!      input OV_ROUGH from OVROUGH.TBL
#ifdef MPP_LAND
       if(my_id .eq. IO_id) then
#endif

       open(71,file="HYDRO.TBL", form="formatted")
!read OV_ROUGH first
          read(71,*) nn
          read(71,*)
          do i = 1, nn
             read(71,*) RT_DOMAIN(did)%OV_ROUGH(i)
          end do
!read parameter for LKSAT
          read(71,*) nn
          read(71,*)
          do i = 1, nn
             read(71,*) xdum1(i), MAXSMC(i),refsmc(i),wltsmc(i)
          end do
       close(71)

#ifdef MPP_LAND
       endif
       call mpp_land_bcast_real(leng,RT_DOMAIN(did)%OV_ROUGH)
       call mpp_land_bcast_real(leng,xdum1)
       call mpp_land_bcast_real(leng,MAXSMC)
       call mpp_land_bcast_real(leng,refsmc)
       call mpp_land_bcast_real(leng,wltsmc)
#endif

       rt_domain(did)%lksat = 0.0
       do j = 1, RT_DOMAIN(did)%jx
             do i = 1, RT_DOMAIN(did)%ix
                rt_domain(did)%lksat(i,j) = xdum1(soltyp(i,j) ) * 1000.0
                IF(rt_domain(did)%VEGTYP(i,j) == 1 ) THEN   ! urban
                    rt_domain(did)%SMCMAX1(i,j) = 0.45
                    rt_domain(did)%SMCREF1(i,j) = 0.42
                    rt_domain(did)%SMCWLT1(i,j) = 0.40
                else
                    rt_domain(did)%SMCMAX1(i,j) = MAXSMC(soltyp(I,J))
                    rt_domain(did)%SMCREF1(i,j) = refsmc(soltyp(I,J))
                    rt_domain(did)%SMCWLT1(i,j) = wltsmc(soltyp(I,J))
                ENDIF
             end do
       end do


      end subroutine lsm_wrf_input

      subroutine  checkSoil(did) 
          implicit none
          integer :: did
          where(rt_domain(did)%smc(:,:,1) <=0) RT_DOMAIN(did)%VEGTYP = 16
          where(rt_domain(did)%sh2ox(:,:,1) <=0) RT_DOMAIN(did)%VEGTYP = 16
          where(rt_domain(did)%smc(:,:,1) >=100) RT_DOMAIN(did)%VEGTYP = 16
          where(rt_domain(did)%sh2ox(:,:,1) >=100) RT_DOMAIN(did)%VEGTYP = 16
      end subroutine checkSoil

end module module_wrf_HYDRO
