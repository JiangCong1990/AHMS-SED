module module_HYDRO_drv
#ifdef MPP_LAND 
   use module_HYDRO_io, only:  output_rt, mpp_output_chrt, mpp_output_lakes, mpp_output_chrtgrd, &
                               restart_out_bi, restart_in_bi
   USE module_mpp_land
#else
   use module_HYDRO_io, only:  output_rt, output_chrt, output_lakes
#endif
   use module_HYDRO_io, only: sub_output_gw, restart_out_nc, restart_in_nc,  &
        get_file_dimension ,get2d_lsm_real, get2d_lsm_vegtyp, get2d_lsm_soltyp, &
        output_lsm, &
        sub_output_rt2d                  ! Qian, 2016.09
   use module_rt_data, only: rt_domain
   use module_GW_baseflow
   use module_gw_gw2d
   use module_gw_gw2d_data, only: gw2d
   use module_channel_routing, only: drive_channel
   use module_namelist, only: nlst_rt, read_rt_nlst
   use module_routing, only: getChanDim, landrt_ini
   use module_HYDRO_utils
!   use module_namelist
   use module_lsm_forcing, only: geth_newdate
   use module_rt_rt2d                    ! Qian, 2016.06
   use module_rt_rt2d_data, ONLY: rt2d   ! Qian, 2016.06 
   
   implicit none

   contains

   subroutine HYDRO_rst_out(did)
      implicit none

      integer :: did
      integer :: rst_out  
      integer :: outflag
      character(len=19) :: out_date
#ifdef MPP_LAND
      character(len=19) :: str_tmp
#endif
      rst_out = -99

#ifdef MPP_LAND
   if(IO_id .eq. my_id) then
#endif
     if(nlst_rt(did)%dt .gt. nlst_rt(did)%rst_dt*60) then
        call geth_newdate(out_date, nlst_rt(did)%startdate, nint(nlst_rt(did)%dt*rt_domain(did)%rst_counts))
     else
        call geth_newdate(out_date, nlst_rt(did)%startdate, nint(nlst_rt(did)%rst_dt*60*rt_domain(did)%rst_counts))
     endif

     if ( (nlst_rt(did)%rst_dt .gt. 0) .and. (out_date(1:19) == nlst_rt(did)%olddate(1:19)) ) then
           rst_out = 99
           rt_domain(did)%rst_counts = rt_domain(did)%rst_counts + 1
     endif

     ! restart every month automatically.
     if ( (nlst_rt(did)%olddate(9:10) == "01") .and. (nlst_rt(did)%olddate(12:13) == "00") .and. &
          (nlst_rt(did)%olddate(15:16) == "00").and. (nlst_rt(did)%olddate(18:19) == "00") .and. &
          (nlst_rt(did)%rst_dt .le. 0)  ) rst_out = 99

#ifdef MPP_LAND
   endif
     call mpp_land_bcast_int1(rst_out)
#endif


    if(rst_out .gt. 0) then

 
#ifdef MPP_LAND
      if(nlst_rt(did)%rst_bi_out .eq. 1) then
             if(my_id .lt. 10) then
                  write(str_tmp,'(I1)') my_id 
             else if(my_id .lt. 100) then 
                  write(str_tmp,'(I2)') my_id 
             else if(my_id .lt. 1000) then
                  write(str_tmp,'(I3)') my_id 
             else if(my_id .lt. 10000) then
                  write(str_tmp,'(I4)') my_id 
             else if(my_id .lt. 100000) then 
                  write(str_tmp,'(I5)') my_id 
             else
                continue
             endif
             call mpp_land_bcast_char(16,nlst_rt(did)%olddate(1:16))
             ! option of nlst_rt(did)%RESTART_OUTNAME is added by Qian, 2016.11
             if (nlst_rt(did)%RESTART_OUTNAME .eq. "") then                           
                call RESTART_OUT_bi(trim("HYDRO_RST."//nlst_rt(did)%olddate(1:16)   &
                     //"_DOMAIN"//trim(nlst_rt(did)%hgrid)//"."//trim(str_tmp)),  did)
             else
                 call RESTART_OUT_bi(trim(nlst_rt(did)%RESTART_OUTNAME)//&
                                     trim(nlst_rt(did)%hgrid)//"_"//trim(str_tmp)//&
                                     "_"//nlst_rt(did)%olddate(1:16), did)
             endif
      else
#endif
             if (nlst_rt(did)%RESTART_OUTNAME .eq. "") then
                call RESTART_OUT_nc(trim("HYDRO_RST."//nlst_rt(did)%olddate(1:16)   &
                     //"_DOMAIN"//trim(nlst_rt(did)%hgrid)),  did)
             else
                call RESTART_OUT_nc(trim(nlst_rt(did)%RESTART_OUTNAME)//&
                                    trim(nlst_rt(did)%hgrid)//&
                                    "_"//nlst_rt(did)%olddate(1:16), did)
             endif
#ifdef MPP_LAND
      endif
#endif

#ifdef MPP_LAND
      if(IO_id .eq. my_id) then
#endif
#ifdef HYDRO_D
        if (nlst_rt(did)%RESTART_OUTNAME .eq. "") then
           write(6,*) "restartFile = ",  "RESTART."//nlst_rt(did)%olddate(1:16)// &
                      "_DOMAIN"//trim(nlst_rt(did)%hgrid)
        else
           write(6,*) "restartFile = ", trim(nlst_rt(did)%RESTART_OUTNAME)//&
                                        trim(nlst_rt(did)%hgrid)//&
                                        "_"//nlst_rt(did)%olddate(1:16)
        endif
#endif
#ifdef MPP_LAND
      endif
#endif

    endif

   end subroutine HYDRO_rst_out

! ---------------------------------------------------------
   subroutine HYDRO_out(did)
      implicit none

      integer :: did
      integer :: outflag, rtflag
      character(len=19) out_date
      integer :: Kt, ounit

!    real, dimension(RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx):: soilmx_tmp, &
!           runoff1x_tmp, runoff2x_tmp, runoff3x_tmp,etax_tmp, &
!           EDIRX_tmp,ECX_tmp,ETTX_tmp,RCX_tmp,HX_tmp,acrain_tmp, &
!           ACSNOM_tmp, esnow2d_tmp, drip2d_tmp,dewfall_tmp, fpar_tmp, &
!           qfx_tmp, prcp_out_tmp, etpndx_tmp

      outflag = -99

#ifdef MPP_LAND
      if(IO_id .eq. my_id) then
#endif
      
      if(nlst_rt(did)%olddate(1:19) .eq. nlst_rt(did)%startdate(1:19) .and. rt_domain(did)%his_out_counts .eq. 0) then
#ifdef HYDRO_D
         write(6,*) "output hydrology at time : ",nlst_rt(did)%olddate(1:19), rt_domain(did)%his_out_counts
#endif
        outflag = 99
      else
         if (nlst_rt(did)%dt .gt. nlst_rt(did)%out_dt*60) then
            call geth_newdate(out_date, nlst_rt(did)%startdate, nint(nlst_rt(did)%dt*rt_domain(did)%out_counts))
         else
            call geth_newdate(out_date, nlst_rt(did)%startdate, nint(nlst_rt(did)%out_dt*60*rt_domain(did)%out_counts))
         endif

         if (out_date(1:19) == nlst_rt(did)%olddate(1:19) ) then
#ifdef HYDRO_D
            write(6,*) "output hydrology at time : ",nlst_rt(did)%olddate(1:19), rt_domain(did)%out_counts
#endif
            outflag = 99
         endif
      endif

#ifdef MPP_LAND
      endif
      call mpp_land_bcast_int1(outflag)
#endif
      call HYDRO_rst_out(did) 

      if (outflag .lt. 0) return

      rt_domain(did)%out_counts     = rt_domain(did)%out_counts + 1
      rt_domain(did)%his_out_counts = rt_domain(did)%his_out_counts + 1

      if(nlst_rt(did)%out_dt*60 .gt. nlst_rt(did)%DT) then
        kt = rt_domain(did)%his_out_counts*nlst_rt(did)%out_dt*60/nlst_rt(did)%DT
      else
        kt = rt_domain(did)%his_out_counts
      endif

      ! jump the ouput for the initial time when it has restart file from routing.
      rtflag = -99
#ifdef MPP_LAND
      if(IO_id .eq. my_id) then
#endif
      if ((trim(nlst_rt(did)%restart_file) /= "") .and. &
          (nlst_rt(did)%startdate(1:19) == nlst_rt(did)%olddate(1:19)) .and. &
          (nlst_rt(did)%rst_bi_in == 1) &  ! added by Qian, 2016.11
         ) then
         print*, "yyyywww restart_file = ", trim(nlst_rt(did)%restart_file) 
!Qian         rtflag = 1 
      endif
#ifdef MPP_LAND
      endif  
      call mpp_land_bcast_int1(rtflag)
#endif

      !yw keep the initial time otuput for debug
      if(rtflag == 1) return ! jump the initial time output for routing restart

      if(nlst_rt(did)%LSMOUT_DOMAN .eq. 1)  &
         call output_lsm(trim(nlst_rt(did)%olddate(1:4)//nlst_rt(did)%olddate(6:7)//nlst_rt(did)%olddate(9:10)  &
                         //nlst_rt(did)%olddate(12:13)//nlst_rt(did)%olddate(15:16)//  &
                         ".LSMOUT_DOMAIN"//trim(nlst_rt(did)%hgrid)), did)

      ! routing switch
      if(     nlst_rt(did)%SUBRTSWCRT  .gt. 0 &
         .or. nlst_rt(did)%OVRTSWCRT   .gt. 0 &
         .or. nlst_rt(did)%GWBASESWCRT .gt. 0 &
         .or. nlst_rt(did)%OVCHRTSWCRT .eq. 1 &
         .or. nlst_rt(did)%INTERGWSFW  .eq. 1 & ! Qian, 2016.08
        ) then

          if(nlst_rt(did)%RTOUT_DOMAIN .eq. 1)  &
            call output_rt(    &
                 nlst_rt(did)%igrid,        nlst_rt(did)%split_output_count, &
                 RT_DOMAIN(did)%ixrt,       RT_DOMAIN(did)%jxrt, &
                 nlst_rt(did)%nsoil,  &
!                nlst_rt(did)%startdate, nlst_rt(did)%olddate, RT_DOMAIN(did)%QSUBRT,&
                 nlst_rt(did)%sincedate,    nlst_rt(did)%olddate, RT_DOMAIN(did)%QSUBRT,&
                 RT_DOMAIN(did)%ZWATTABLRT, RT_DOMAIN(did)%SMCRT,&
                 RT_DOMAIN(did)%SUB_RESID,       &
                 RT_DOMAIN(did)%q_sfcflx_x, RT_DOMAIN(did)%q_sfcflx_y,&
                 RT_DOMAIN(did)%soxrt,      RT_DOMAIN(did)%soyrt,&
                 RT_DOMAIN(did)%QSTRMVOLRT, RT_DOMAIN(did)%SFCHEADSUBRT, &
                 nlst_rt(did)%geo_finegrid_flnm, nlst_rt(did)%DT,&
                 RT_DOMAIN(did)%SLDPTH,     RT_DOMAIN(did)%LATVAL,&
                 RT_DOMAIN(did)%LONVAL,     RT_DOMAIN(did)%dist, &
                 nlst_rt(did)%RTOUT_DOMAIN, RT_DOMAIN(did)%QBDRYRT )

          if(nlst_rt(did)%GWBASESWCRT .eq. 3) then
            if(nlst_rt(did)%output_gw  .eq. 1)  &
              call sub_output_gw(    &
                   nlst_rt(did)%igrid,     nlst_rt(did)%split_output_count, &
                   RT_DOMAIN(did)%ixrt,    RT_DOMAIN(did)%jxrt,    &
                   nlst_rt(did)%nsoil,                             &
!                  nlst_rt(did)%startdate, nlst_rt(did)%olddate,   &
                   nlst_rt(did)%sincedate, nlst_rt(did)%olddate,   &
                   gw2d(did)%h,            RT_DOMAIN(did)%SMCRT,   &
                   gw2d(did)%convgw,       gw2d(did)%excess,       &
                   gw2d(did)%qsgwrt,       gw2d(did)%qgw_chanrt,   &
                   nlst_rt(did)%geo_finegrid_flnm,nlst_rt(did)%DT, &
                   RT_DOMAIN(did)%SLDPTH,  RT_DOMAIN(did)%LATVAL,  &
                   RT_DOMAIN(did)%LONVAL,  rt_domain(did)%dist,    &
                   nlst_rt(did)%output_gw, nlst_rt(did)%GW2D_OUTNAME)
          endif

          ! output rt2d results, Qian, 2016.09
          if ((nlst_rt(did)%OVCHRTSWCRT .eq. 1) .or. &
              (nlst_rt(did)%INTERGWSFW  .eq. 1)) then   
            if(nlst_rt(did)%output_rt2d .eq. 1) then
              call sub_output_rt2d( &
                   nlst_rt(did)%igrid,       nlst_rt(did)%split_output_count, &
                   RT_DOMAIN(did)%ixrt,      RT_DOMAIN(did)%jxrt,    &
                   nlst_rt(did)%nsoil,       nlst_rt(did)%sincedate, &
                   nlst_rt(did)%olddate,     rt2d(did)%dt,           &
                   nlst_rt(did)%output_rt2d, nlst_rt(did)%RT2D_OUTNAME, &
                   rt2d(did)%hend,                                   &
                   rt2d(did)%runoff,         rt2d(did)%rsat,         &
                   rt2d(did)%sinfil,         rt2d(did)%tinfil,       &
                   rt2d(did)%fwater,         rt2d(did)%fsat,         &
                   rt2d(did)%cuiflux,        rt2d(did)%cgiflux,      &
                   rt2d(did)%cglflux,        rt2d(did)%hlake,        &
                   rt2d(did)%streamflux,     rt2d(did)%istreamflux,  &
                   rt2d(did)%w_us,           rt2d(did)%w_gw,         &           
                   rt2d(did)%w_sf,           rt2d(did)%ac_pmewsf,    &
                   rt2d(did)%ac_infil_grnd,  rt2d(did)%ac_drain_grnd,&
                   rt2d(did)%ac_roff_grnd,   rt2d(did)%ac_rsat_grnd, &
                   rt2d(did)%ac_precip,      rt2d(did)%ac_evap,      &
                   rt2d(did)%ac_cgi,         rt2d(did)%ac_cgl,       &
                   rt2d(did)%ac_cui,         rt2d(did)%ac_ocngw,     &
                   rt2d(did)%ac_ocnsf,       rt2d(did)%ac_fixgw,     &
                   rt2d(did)%ac_fixsf,       rt2d(did)%ltype,        &
                   rt2d(did)%elev,           rt2d(did)%bank,         &
                   rt2d(did)%width,                                  &
                   rt2d(did)%bot,            rt2d(did)%fbed,         &
                   rt2d(did)%topo,           rt2d(did)%irg,          &
                   rt2d(did)%sabs,           rt2d(did)%ac_precevap,  &
                   rt2d(did)%ac_irgrw,       rt2d(did)%ac_irggw,     &
                   rt2d(did)%veloc,                                  &
                   rt2d(did)%qsKR,           nlst_rt(did)%NumFracts,    &
                   rt2d(did)%SusSed,         rt2d(did)%DepSed,       &
                   rt2d(did)%ErodSed,        rt2d(did)%SedFlux,      &
                   rt2d(did)%t_SedFlux,      rt2d(did)%qsEH,         &
                   rt2d(did)%SedConc,        rt2d(did)%T_SedConc)
            endif

          endif
     
          if(nlst_rt(did)%CHANRTSWCRT.eq.1 .or. nlst_rt(did)%CHANRTSWCRT.eq.2) then 
            if(nlst_rt(did)%CHRTOUT_DOMAIN .eq. 1)  &
#ifdef MPP_LAND
               call mpp_output_chrt(rt_domain(did)%gnlinks,rt_domain(did)%map_l2g, &
#else
               call output_chrt(  &
#endif
               nlst_rt(did)%igrid, nlst_rt(did)%split_output_count, &
               RT_DOMAIN(did)%NLINKS,RT_DOMAIN(did)%ORDER, &
!              nlst_rt(did)%startdate,nlst_rt(did)%olddate,RT_DOMAIN(did)%CHLON,&
               nlst_rt(did)%sincedate,nlst_rt(did)%olddate,RT_DOMAIN(did)%CHLON,&
               RT_DOMAIN(did)%CHLAT, &
               RT_DOMAIN(did)%HLINK, RT_DOMAIN(did)%ZELEV, &
               RT_DOMAIN(did)%QLINK,nlst_rt(did)%DT,Kt, &
               RT_DOMAIN(did)%STRMFRXSTPTS,nlst_rt(did)%order_to_write)

      
#ifdef MPP_LAND
            if(nlst_rt(did)%CHRTOUT_GRID  .eq. 1)  &
              call mpp_output_chrtgrd(nlst_rt(did)%igrid, nlst_rt(did)%split_output_count, &
                   RT_DOMAIN(did)%ixrt,RT_DOMAIN(did)%jxrt, RT_DOMAIN(did)%NLINKS,   &
                   RT_DOMAIN(did)%GCH_NETLNK, &
                   nlst_rt(did)%startdate, nlst_rt(did)%olddate, &
                   RT_DOMAIN(did)%qlink, nlst_rt(did)%dt, nlst_rt(did)%geo_finegrid_flnm,   &
                   RT_DOMAIN(did)%gnlinks,RT_DOMAIN(did)%map_l2g,                   &
                   RT_DOMAIN(did)%g_ixrt,RT_DOMAIN(did)%g_jxrt )
#endif

            if(RT_DOMAIN(did)%NLAKES .gt. 0)  then
              if(nlst_rt(did)%outlake  .eq. 1)  &
#ifdef MPP_LAND
              call mpp_output_lakes( RT_DOMAIN(did)%lake_index, &
#else
              call output_lakes(  &
#endif
                   nlst_rt(did)%igrid, nlst_rt(did)%split_output_count, &
                   RT_DOMAIN(did)%NLAKES, &
                   trim(nlst_rt(did)%sincedate), trim(nlst_rt(did)%olddate), &
                   RT_DOMAIN(did)%LATLAKE,RT_DOMAIN(did)%LONLAKE, &
                   RT_DOMAIN(did)%ELEVLAKE,RT_DOMAIN(did)%QLAKEI, &
                   RT_DOMAIN(did)%QLAKEO, &
                   RT_DOMAIN(did)%RESHT,nlst_rt(did)%DT,Kt)
            endif
          endif

#ifdef HYDRO_D
            write(6,*) "end calling output functions"
#endif
          endif  ! end of routing switch

#ifdef HYDRO_D
          write(6,*) "finish HYDRO_out"
#endif

      end subroutine HYDRO_out

! ---------------------------------------------------------
      subroutine HYDRO_rst_in(did)
        integer :: did
        integer:: flag 

   flag = -1
#ifdef MPP_LAND
   if(my_id.eq.IO_id) then
#endif
     if (trim(nlst_rt(did)%restart_file) /= "") then
        flag = 99
        rt_domain(did)%timestep_flag = 99   ! continue run
     endif 
#ifdef MPP_LAND
   endif 
   call mpp_land_bcast_int1(flag)
#endif

   nlst_rt(did)%sincedate = nlst_rt(did)%startdate
   
   if (flag.eq.99) then
#ifdef MPP_LAND
     if(my_id.eq.IO_id) then
#endif
#ifdef HYDRO_D
       write(6,*) "*** read restart data: ",trim(nlst_rt(did)%restart_file)
#endif
#ifdef MPP_LAND
     endif 
#endif

#ifdef MPP_LAND
     if(nlst_rt(did)%rst_bi_in .eq. 1) then
       call RESTART_IN_bi(trim(nlst_rt(did)%restart_file), did)
     else
#endif
       call RESTART_IN_nc(trim(nlst_rt(did)%restart_file), did)
#ifdef MPP_LAND
     endif
#endif

!yw  if (trim(nlst_rt(did)%restart_file) /= "") then 
!yw          nlst_rt(did)%restart_file = ""
!yw  endif

   endif
 end subroutine HYDRO_rst_in

! ---------------------------------------------------------
     subroutine HYDRO_time_adv(did)
        implicit none

        integer :: did
        character(len = 19) :: newdate 
 
#ifdef MPP_LAND
   if (IO_id.eq.my_id) then
#endif
        !Qian call geth_newdate(newdate, nlst_rt(did)%olddate, nint( nlst_rt(did)%dt))
        call geth_newdate(newdate, nlst_rt(did)%olddate, nint( nlst_rt(did)%dtrt))
        nlst_rt(did)%olddate = newdate
#ifdef HYDRO_D
        write(6,*) "current time is ",newdate
#endif
#ifdef MPP_LAND
   endif
#endif
     end subroutine HYDRO_time_adv

! ----------------------------------------------------  
   subroutine HYDRO_exe(did)
     implicit none

     integer, intent(in) :: did
     integer             :: rst_out

       call HYDRO_out(did)

! running land surface model
! cpl: 0-- offline run; 
!      1-- coupling with WRF but running offline lsm; 
!      2-- coupling with WRF but do not run offline lsm  
!      3-- coupling with LIS and do not run offline lsm  
!      4-- coupling with CLM

!       if (nlst_rt(did)%SYS_CPL .eq. 0 .or. nlst_rt(did)%SYS_CPL .eq. 1 )then
!          call drive_noahLSF(did,kt)
!       else
!          ! does not run the NOAH LASF model, only read the parameter
!          call read_land_par(did,lsm(did)%ix,lsm(did)%jx)
!       endif

        if (    nlst_rt(did)%GWBASESWCRT      .NE. 0 &   
           .or. nlst_rt(did)%SUBRTSWCRT       .NE. 0 &
           .or. nlst_rt(did)%OVRTSWCRT        .NE. 0 &
           .or. nlst_rt(did)%OVCHRTSWCRT      .EQ. 1 &     ! Qian, 2016.06
           .or. nlst_rt(did)%INTERGWSFW       .EQ. 1 &     ! Qian, 2016.10
           ) THEN

           if ((nlst_rt(did)%OVCHRTSWCRT .EQ. 0) .and. & 
               (nlst_rt(did)%INTERGWSFW  .EQ. 0)) then       ! Qian, 2016.06
              RT_DOMAIN(did)%QSTRMVOLRT_DUM   = RT_DOMAIN(did)%QSTRMVOLRT
              RT_DOMAIN(did)%LAKE_INFLORT_DUM = RT_DOMAIN(did)%LAKE_INFLORT
           endif

           ! step 1) disaggregate specific fields from LSM to Hydro grid
           call disaggregateDomain_drv(did)

           ! step 2) Subsurface Routing  
           call SubsurfaceRouting_drv(did)

           ! step 3) Overland Routing
           call OverlandRouting_drv(did) 

           if ((nlst_rt(did)%OVCHRTSWCRT .EQ. 0) .and. &
               (nlst_rt(did)%INTERGWSFW  .EQ. 0)) then              
              RT_DOMAIN(did)%QSTRMVOLRT_TS   = RT_DOMAIN(did)%QSTRMVOLRT - RT_DOMAIN(did)%QSTRMVOLRT_DUM
              RT_DOMAIN(did)%LAKE_INFLORT_TS = RT_DOMAIN(did)%LAKE_INFLORT - RT_DOMAIN(did)%LAKE_INFLORT_DUM
           endif

           ! step 4) Groudnwater Routing
           call driveGwBaseflow(did)

           ! step 4-2) Interaction of Groundwater and Surfacewater
           call driveInteractionGwSfw(did)
 
           ! step 5) Channel Routing
           call driveChannelRouting(did)

           ! step 5-2) 2D Overland Channel Routing
           call driveOverlandChannelRouting(did)  ! Qian, 2016.06  

           ! step 6) aggregate specific fields from Hydro to LSM grid
           call aggregateDomain(did)

        end if


           ! advance to next time step
           call HYDRO_time_adv(did)

           ! output for history 
           call HYDRO_out(did)

!          write(90 + my_id,*) "finish calling hydro_exe"
!          call flush(90+my_id)
!          call mpp_land_sync()
           
           RT_DOMAIN(did)%SOLDRAIN = 0. 
           if ((nlst_rt(did)%OVCHRTSWCRT .EQ. 0) .and. & 
               (nlst_rt(did)%INTERGWSFW  .EQ. 0)) then  ! Qian, 2016.08
              RT_DOMAIN(did)%QSUBRT = 0
           endif

   end subroutine HYDRO_exe
      
!----------------------------------------------------------      
   subroutine driveGwBaseflow(did)
     implicit none

     integer, intent(in) :: did       
     integer             :: i, j, jj, ii

!------------------------------------------------------------------
!DJG Begin GW/Baseflow Routines
!-------------------------------------------------------------------

  IF (nlst_rt(did)%GWBASESWCRT.GE.1) THEN           ! Switch to activate/specify GW/Baseflow

!  IF (nlst_rt(did)%GWBASESWCRT.GE.1000) THEN       ! Switch to activate/specify GW/Baseflow

    If (     nlst_rt(did)%GWBASESWCRT .EQ. 1 &
       .OR. nlst_rt(did)%GWBASESWCRT  .EQ. 2 ) Then ! Call simple bucket baseflow scheme

#ifdef HYDRO_D
       write(6,*) "*****yw******start simp_gw_buck "
#endif

       call simp_gw_buck(RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx,RT_DOMAIN(did)%ixrt,&
             RT_DOMAIN(did)%jxrt,RT_DOMAIN(did)%numbasns,RT_DOMAIN(did)%basns_area,&
             RT_DOMAIN(did)%gwsubbasmsk, RT_DOMAIN(did)%INFXSRT, &
             RT_DOMAIN(did)%SOLDRAIN, &
             RT_DOMAIN(did)%z_gwsubbas,&
             RT_DOMAIN(did)%qin_gwsubbas,RT_DOMAIN(did)%qout_gwsubbas,&
             RT_DOMAIN(did)%qinflowbase,&
             RT_DOMAIN(did)%gw_strm_msk,RT_DOMAIN(did)%gwbas_pix_ct, &
             RT_DOMAIN(did)%dist,nlst_rt(did)%DT,&
             RT_DOMAIN(did)%gw_buck_coeff,RT_DOMAIN(did)%gw_buck_exp, &
             RT_DOMAIN(did)%z_max,&
             nlst_rt(did)%GWBASESWCRT,nlst_rt(did)%OVRTSWCRT)

#ifdef MPP_LAND
       if(my_id .eq. IO_id) then
#endif
         open (unit=51,file='GW_inflow.txt',form='formatted',&
               status='unknown',position='append')
         open (unit=52,file='GW_outflow.txt',form='formatted',&
               status='unknown',position='append')
         open (unit=53,file='GW_zlev.txt',form='formatted',&
               status='unknown',position='append')
         do i=1,RT_DOMAIN(did)%numbasns
            write (51,951) i,nlst_rt(did)%olddate,rt_domain(did)%qin_gwsubbas(i)
951         FORMAT(I3,1X,A19,1X,F11.3)
            write (52,951) i,nlst_rt(did)%olddate,rt_domain(did)%qout_gwsubbas(i)
            write (53,951) i,nlst_rt(did)%olddate,rt_domain(did)%z_gwsubbas(i)
         end do
         close(51)
         close(52)
         close(53)
#ifdef MPP_LAND
       endif
#endif

#ifdef HYDRO_D 
       write(6,*) "*****yw******end simp_gw_buck "
#endif

    ! For parameter setup runs output the percolation for each basin,
    ! otherwise comment out this output...
    Else If (nlst_rt(did)%gwBaseSwCRT .eq. 3) Then

#ifdef HYDRO_D
       write(6,*) "*** start 2d_gw_model..."
#endif

       ! compute qsgwrt  between lsm and gw with namelist selected coupling method
       ! qsgwrt is defined on the routing grid  and needs to be aggregated for SFLX
       if (nlst_rt(did)%gwsoilcpl .GT. 0) then     
          call gwSoilFlux(did)   
       end if

       ! 2D groundwater routing   
       gw2d(did)%excess = 0.
       call gwstep( gw2d(did)%ix,     gw2d(did)%jx,       gw2d(did)%dx,      &
                    gw2d(did)%ltype,  gw2d(did)%elev,     gw2d(did)%bot,     &
                    gw2d(did)%hycond, gw2d(did)%poros,    gw2d(did)%compres, &
                    gw2d(did)%ho,     gw2d(did)%h,        gw2d(did)%convgw,  &
                    gw2d(did)%excess, gw2d(did)%ebot,     gw2d(did)%eocn,    & 
                    gw2d(did)%aqthick,gw2d(did)%dt,       gw2d(did)%istep,   &
                    rt2d(did)%hend,   rt2d(did)%ac_ocngw, rt2d(did)%ac_fixgw &  ! Qian, 2016.09
                  )
       gw2d(did)%ho = gw2d(did)%h

       ! put surface exceeding groundwater to surface routing inflow
       if (nlst_rt(did)%OVCHRTSWCRT .EQ. 0) &
          RT_DOMAIN(did)%SFCHEADSUBRT = RT_DOMAIN(did)%SFCHEADSUBRT &
	                              + gw2d(did)%excess*1000. ! convert to mm
	  
       ! aggregate qsgw from routing to lsm grid
       call aggregateQsgw(did)

       gw2d(did)%istep =  gw2d(did)%istep + 1
  
#ifdef HYDRO_D
       write(6,*) "*** finish 2d_gw_model "
#endif
      
    End if

  END IF    !DJG (End if for RTE SWC activation)

   end subroutine driveGwBaseflow

!-------------------------------------------
   subroutine driveInteractionGwSfw(did)         ! Qian, 2016.10
     implicit none

     integer, intent(in) :: did

     IF (nlst_rt(did)%INTERGWSFW .EQ. 1) THEN
#ifdef HYDRO_D
#ifdef MPP_LAND
       if (my_id .eq. IO_id) &
#endif
       write(6,*) "*** start Interaction of surface water and groundwater..."
#endif

       if (nlst_rt(did)%OVCHRTSWCRT .EQ. 1) then  ! Qian, 2016.06

          ! fractional water surface area in grid
          call fwater_cal(did)                    ! rt2d(did)%fwater

          ! fracitonal saturated soil in grid
          call fsat_cal(did)                      ! rt2d(did)%fsat

          ! estimate flux between groundwater and stream
          call flux_gw_chan(did, rt2d(did)%dt)    ! rt2d(did)%cgiflux & rt2d(did)%cuiflux & 
                                                  ! rt2d(did)%ac_cgi  & rt2d(did)%ac_cui

          ! estimate flux between groundwater and lake
          call flux_gw_lake(did, rt2d(did)%dt)    ! rt2d(did)%cglflux & rt2d(did)%ac_cgl

          ! compute water into/out of soil     
          call dvolgw_cal(did, rt2d(did)%dt)      ! rt2d(did)%dvolgw  & rt2d(did)%ac_infil_grnd
                                                  ! rt2d(did)%ac_drain_grnd

          ! update groundwater head & soil moisture
          call update_gwhead_soil(did, rt2d(did)%dt) ! rt2d(did)%hend 
                                                     ! rt2d(did)%rsat & rt2d(did)%ac_rsat_grnd

       end if

#ifdef HYDRO_D
#ifdef MPP_LAND
       if (my_id .eq. IO_id) &
#endif
     write(6,*) "*** finish Interaction"
#endif

     ENDIF

   end subroutine driveInteractionGwSfw

!-------------------------------------------
   subroutine driveOverlandChannelRouting(did)   ! Qian, 2016.06
     implicit none

     integer, intent(in) :: did
     integer             :: i, j
     real                :: zwei, timeavnew !!C.Jiang
     real,SAVE           :: timeav

     IF (       nlst_rt(did)%GWBASESWCRT      .EQ. 3 &
           .or. nlst_rt(did)%OVCHRTSWCRT      .EQ. 1 &     ! Qian, 2016.06
           .or. nlst_rt(did)%INTERGWSFW       .EQ. 1 &     ! Qian, 2016.10
        ) THEN

        do j = 1, rt2d(did)%jx
          do i = 1, rt2d(did)%ix

             ! accumulated precipitation, m
             rt2d(did)%ac_precip(i,j) = rt2d(did)%ac_precip(i,j) + rt2d(did)%precip(i,j) * rt2d(did)%dt

             ! accumulated evaporation, m
             rt2d(did)%ac_evap(i,j)   = rt2d(did)%ac_evap(i,j) + rt2d(did)%evap(i,j) * rt2d(did)%dt


             ! accumulated net infiltration, m 
!            rt2d(did)%ac_infil_grnd(i,j) = rt2d(did)%ac_infil_grnd(i,j) + rt2d(did)%dt * &
!                                           (1.0 - rt2d(did)%fbed(i,j)) * rt2d(did)%tinfil(i,j)

             rt2d(did)%ac_infil_grnd(i,j) = rt2d(did)%ac_infil_grnd(i,j) + rt2d(did)%tinfil(i,j) * rt2d(did)%dt


             ! accumulated soil drainage, m
!            rt2d(did)%ac_drain_grnd(i,j) = rt2d(did)%ac_drain_grnd(i,j) + rt2d(did)%dt * &
!                                           (1.0 - rt2d(did)%fbed(i,j)) * rt2d(did)%sinfil(i,j)

            rt2d(did)%ac_drain_grnd(i,j) = rt2d(did)%ac_drain_grnd(i,j) + rt2d(did)%sinfil(i,j) * rt2d(did)%dt


            ! accumulated runoff, m 
!           ac_roff_grnd(i,j) = ac_roff_grnd(i,j) + (1. - fbed(i,j)) * runoff(i,j) * rt2d(did)%dt
            rt2d(did)%ac_roff_grnd(i,j) = rt2d(did)%ac_roff_grnd(i,j) + rt2d(did)%runoff(i,j) * rt2d(did)%dt

          end do
        end do

     ENDIF

     IF (nlst_rt(did)%OVCHRTSWCRT .EQ. 1) THEN

#ifdef HYDRO_D
#ifdef MPP_LAND
       if (my_id .eq. IO_id) &
#endif
       write(6,*) "*** start OverlandChannelRouting..."
#endif

        call drive_OverlandChannelRouting (   &
          rt2d(did)%ix,                   rt2d(did)%jx,                rt2d(did)%dx,  &
          rt2d(did)%ltype,                rt2d(did)%ielev,             rt2d(did)%elev, &
          rt2d(did)%bank,                 rt2d(did)%rough,             rt2d(did)%width, &
          rt2d(did)%fbed,                 rt2d(did)%fwater,            rt2d(did)%fsat, &
          rt2d(did)%hlake,                rt2d(did)%cgiflux,           rt2d(did)%cuiflux, &
          rt2d(did)%cglflux,              rt2d(did)%rsat,              rt2d(did)%runoff,  &
          rt2d(did)%precip,               rt2d(did)%evap,              rt2d(did)%streamflux, &
          rt2d(did)%istreamflux,          rt2d(did)%dt,                rt2d(did)%istep, &
          nlst_rt(did)%DT_2DTHRT,         nlst_rt(did)%OVERCHANRT_METHOD, &
          nlst_rt(did)%ntran,             rt2d(did)%ncross,            rt2d(did)%itran, &
          rt2d(did)%jtran,                rt2d(did)%ktran,             rt2d(did)%fluxtran, &
          rt2d(did)%fluxtranav,           nlst_rt(did)%dt_flowtran,    rt2d(did)%count_flowout,&
          rt2d(did)%ac_ocnsf,             rt2d(did)%ac_fixsf,                              &
          rt2d(did)%ac_pmewsf,            rt2d(did)%ac_roff_grnd,      rt2d(did)%w_sf, &
          rt2d(did)%sdep,                 did,                         rt2d(did)%veloc,&
          rt2d(did)%qsKR,                 nlst_rt(did)%indexeros,      nlst_rt(did)%NumFracts,&
          rt2d(did)%SusSed,               rt2d(did)%DepSed,            rt2d(did)%ErodSed,&
          rt2d(did)%SedFlux,              rt2d(did)%t_SedFlux,         rt2d(did)%qsEH,&
          rt2d(did)%SedConc,              rt2d(did)%T_SedConc,         rt2d(did)%sedfluxtranav, &
          rt2d(did)%T_SedConctranav)

 ! C.Jiang
 ! Update exponential-decay running-mean discharge (disch)
 ! for current values
 !         zwei=exp(-rt2d(did)%dt/(180*86400))
 !            do j = 1, rt2d(did)%jx
 !              do i = 1, rt2d(did)%ix
! 
!               if (rt2d(did)%ltype(i,j).eq.1) then
! 
!              rt2d(did)%disch(i,j)=zwei *rt2d(did)%disch(i,j)+(1.-zwei)*rt2d(did)%streamflux(i,j)
!              rt2d(did)%bank(i,j) =  0.6*(rt2d(did)%disch(i,j)**0.4)
!              rt2d(did)%bank(i,j) = min(30., max(0.5, rt2d(did)%bank(i,j)))
!              rt2d(did)%width(i,j) = 10.*(rt2d(did)%disch(i,j)**0.5)
!              rt2d(did)%width(i,j) = min(rt2d(did)%width(i,j), rt2d(did)%dx)
!             endif
! 
!            enddo
!          enddo

!    
     if (rt2d(did)%istep.eq.0) then
      timeav=0.
     endif
     
       
!       Update (annual) time averages for later use in gw spinup mode

      timeavnew = timeav + rt2d(did)%dt
      do j=1,rt2d(did)%jx
        do i=1,rt2d(did)%ix
          if (rt2d(did)%ltype(i,j).eq.1) then
           rt2d(did)%hlakeav(i,j) = (timeav*rt2d(did)%hlakeav (i,j) + rt2d(did)%dt*rt2d(did)%hlake (i,j)) &
                           / timeavnew
           rt2d(did)%sinfilav(i,j) = (timeav*rt2d(did)%sinfilav(i,j) + rt2d(did)%dt*rt2d(did)%sinfil(i,j)) &
                           / timeavnew
          endif
        enddo
      enddo
      timeav = timeavnew

!--------------------------------------------------------------------------------------------------

        if (nlst_rt(did)%INTERGWSFW .EQ. 1) then

           ! After hlake changes, reset groundwater head in lakes to be water level 
           do j = 1, rt2d(did)%jx
              do i = 1, rt2d(did)%ix
                 if (rt2d(did)%ltype(i,j) .eq. 1) then
                    if (rt2d(did)%hend(i,j) .ge. rt2d(did)%elev(i,j)) then
                       rt2d(did)%hend(i,j) = max(rt2d(did)%hlake(i,j), rt2d(did)%elev(i,j))
                    endif 
                 endif
              end do
           end do

           ! update water storage in unsaturated and saturated soil zone
!           call water_storage(did, rt2d(did)%ix, rt2d(did)%jx, rt2d(did)%hend, rt2d(did)%w_us, rt2d(did)%w_gw)
           rt2d(did)%hbeg = rt2d(did)%hend   ! gwhead_begin(t+1) = gwhead_end(t)

           ! finally update groundwater head
           if (nlst_rt(did)%gwBaseSwCRT .eq. 3) then 
              gw2d(did)%h  = rt2d(did)%hend  
              gw2d(did)%ho = gw2d(did)%h     ! gwhead_begin(t+1) = gwhead_end(t)
           endif

        endif

        do j = 1, rt2d(did)%jx
          do i = 1, rt2d(did)%ix
             ! surface water above land ground is considered as surface water head (mm)
             ! RT_DOMAIN(did)%SFCHEADSUBRT -> RT_DOMAIN(did)%SFCHEADRT -> grid(did)%SFCHEADRT in WRF
             RT_DOMAIN(did)%SFCHEADSUBRT(i,j) = max(0., rt2d(did)%hlake(i,j) - rt2d(did)%elev(i,j)) * 1000.

          end do
        end do
 
        rt2d(did)%istep = rt2d(did)%istep + 1

#ifdef HYDRO_D
#ifdef MPP_LAND
       if (my_id .eq. IO_id) &
#endif
     write(6,*) "*** finish OverlandChannelRouting"
#endif

     ENDIF

   end subroutine driveOverlandChannelRouting
      
! ---------------------------------------------------------      
      subroutine driveChannelRouting(did)
       
       implicit none
       integer, intent(in) :: did
       
  !DJG,DNY  Begin Channel and Lake Routing Routines
  IF ( nlst_rt(did)%CHANRTSWCRT.EQ.1 .or. nlst_rt(did)%CHANRTSWCRT.EQ.2 ) THEN
 
    call drive_CHANNEL(RT_DOMAIN(did)%latval,RT_DOMAIN(did)%lonval, &
       RT_DOMAIN(did)%timestep_flag,RT_DOMAIN(did)%IXRT,RT_DOMAIN(did)%JXRT, &
       nlst_rt(did)%SUBRTSWCRT, RT_DOMAIN(did)%QSUBRT, &
       RT_DOMAIN(did)%LAKE_INFLORT_TS, RT_DOMAIN(did)%QSTRMVOLRT_TS,&
       RT_DOMAIN(did)%TO_NODE, RT_DOMAIN(did)%FROM_NODE, RT_DOMAIN(did)%TYPEL,&
       RT_DOMAIN(did)%ORDER, RT_DOMAIN(did)%MAXORDER, RT_DOMAIN(did)%NLINKS,&
       RT_DOMAIN(did)%CH_NETLNK, RT_DOMAIN(did)%CH_NETRT, &
       RT_DOMAIN(did)%LAKE_MSKRT, nlst_rt(did)%DT, nlst_rt(did)%DTCT,nlst_rt(did)%DTRT,&
       RT_DOMAIN(did)%MUSK, RT_DOMAIN(did)%MUSX,  RT_DOMAIN(did)%QLINK, &
       RT_DOMAIN(did)%HLINK, RT_DOMAIN(did)%ELRT,RT_DOMAIN(did)%CHANLEN,&
       RT_DOMAIN(did)%MannN,RT_DOMAIN(did)%So, RT_DOMAIN(did)%ChSSlp, &
       RT_DOMAIN(did)%Bw,&
       RT_DOMAIN(did)%RESHT, RT_DOMAIN(did)%HRZAREA, RT_DOMAIN(did)%LAKEMAXH,&
       RT_DOMAIN(did)%WEIRC, RT_DOMAIN(did)%WEIRL, RT_DOMAIN(did)%ORIFICEC, &
       RT_DOMAIN(did)%ORIFICEA, &
       RT_DOMAIN(did)%ORIFICEE, RT_DOMAIN(did)%ZELEV, RT_DOMAIN(did)%CVOL, &
       RT_DOMAIN(did)%NLAKES, RT_DOMAIN(did)%QLAKEI, RT_DOMAIN(did)%QLAKEO,&
       RT_DOMAIN(did)%LAKENODE, RT_DOMAIN(did)%dist, &
       RT_DOMAIN(did)%QINFLOWBASE, RT_DOMAIN(did)%CHANXI, &
       RT_DOMAIN(did)%CHANYJ, nlst_rt(did)%channel_option, &
       RT_DOMAIN(did)%RETDEP_CHAN &
       , RT_DOMAIN(did)%node_area &
#ifdef MPP_LAND
       ,RT_DOMAIN(did)%lake_index,RT_DOMAIN(did)%link_location,&
       RT_DOMAIN(did)%mpp_nlinks,RT_DOMAIN(did)%nlinks_index, &
       RT_DOMAIN(did)%yw_mpp_nlinks  &
#endif
       ,nlst_rt(did)%GwBaseSwCRT, gw2d(did)%ho, gw2d(did)%qgw_chanrt, &
       nlst_rt(did)%gwChanCondSw, nlst_rt(did)%gwChanCondConstIn, &
       nlst_rt(did)%gwChanCondConstOut &
       )
    
    if ((nlst_rt(did)%gwBaseSwCRT .eq. 3) .and. (nlst_rt(did)%gwChanCondSw .eq. 1)) then
      
      ! add/rm channel-aquifer exchange contribution    
      gw2d(did)%ho  = gw2d(did)%ho  &
                    + (((gw2d(did)%qgw_chanrt*(-1)) * gw2d(did)%dt / gw2d(did)%dx**2) &
                      /  gw2d(did)%poros)
      
    endif

#ifdef HYDRO_D
    write(6,*) "*****yw******end drive_CHANNEL "
#endif

  ENDIF
      
      end subroutine driveChannelRouting
 
!---------------------------------------------------------- 
      subroutine aggregateDomain(did)
       
       implicit none
       integer, intent(in) :: did

       integer :: i, j, krt, ixxrt, jyyrt, &
                  AGGFACYRT, AGGFACXRT

       real    :: SICEAGGRT(nlst_rt(did)%NSOIL)  ! Qian, 2016.12
       real    :: IRRIGAGGRT                     ! Cong, 2022.08
       real    :: peta, neta  ! for deviation

       neta = -1.E-6
       peta = 1.E-6
#ifdef HYDRO_D
       write(6,*), "*** Beginning Aggregation..."
#endif

        do J=1,RT_DOMAIN(did)%JX
          do I=1,RT_DOMAIN(did)%IX

             RT_DOMAIN(did)%SFCHEADAGGRT = 0.
!DJG Subgrid weighting edit...
             RT_DOMAIN(did)%LSMVOL=0.
             do KRT=1,nlst_rt(did)%NSOIL
!                SMCAGGRT(KRT) = 0.
               SICEAGGRT(KRT) = 0.
               RT_DOMAIN(did)%SH2OAGGRT(KRT) = 0.
             end do
             IRRIGAGGRT = 0.

             do AGGFACYRT=nlst_rt(did)%AGGFACTRT-1,0,-1
              do AGGFACXRT=nlst_rt(did)%AGGFACTRT-1,0,-1


                IXXRT=I*nlst_rt(did)%AGGFACTRT-AGGFACXRT
                JYYRT=J*nlst_rt(did)%AGGFACTRT-AGGFACYRT
#ifdef MPP_LAND
                if(left_id.ge.0) IXXRT=IXXRT+1
                if(down_id.ge.0) JYYRT=JYYRT+1
#else
!yw ????
!       IXXRT=IXXRT+1
!       JYYRT=JYYRT+1
#endif

!State Variables
                RT_DOMAIN(did)%SFCHEADAGGRT = RT_DOMAIN(did)%SFCHEADAGGRT &
                                            + RT_DOMAIN(did)%SFCHEADSUBRT(IXXRT,JYYRT)

!DJG Subgrid weighting edit...
                RT_DOMAIN(did)%LSMVOL = RT_DOMAIN(did)%LSMVOL &
                                      + RT_DOMAIN(did)%SFCHEADSUBRT(IXXRT,JYYRT) &
                                      * RT_DOMAIN(did)%dist(IXXRT,JYYRT,9)

                do KRT=1,nlst_rt(did)%NSOIL
!DJG              SMCAGGRT(KRT)=SMCAGGRT(KRT)+SMCRT(IXXRT,JYYRT,KRT)
                  SICEAGGRT(KRT) = SICEAGGRT(KRT) + max(0., RT_DOMAIN(did)%SICERT(IXXRT,JYYRT,KRT))
                  RT_DOMAIN(did)%SH2OAGGRT(KRT) = RT_DOMAIN(did)%SH2OAGGRT(KRT) &
                                                + max(0., RT_DOMAIN(did)%SMCRT(IXXRT,JYYRT,KRT))
                end do
                 
! aggregate irrigation amount

                IRRIGAGGRT = IRRIGAGGRT + rt2d(did)%actirg(IXXRT,JYYRT)

              end do
             end do
! update constrained irrigation amount

            RT_DOMAIN(did)%IRRIG(I,J) = IRRIGAGGRT/(nlst_rt(did)%AGGFACTRT**2) * (1000.*nlst_rt(did)%dtrt)  ! m/s --> mm

            RT_DOMAIN(did)%SFCHEADRT(I,J) = RT_DOMAIN(did)%SFCHEADAGGRT &
                                          / (nlst_rt(did)%AGGFACTRT**2)

            do KRT=1,nlst_rt(did)%NSOIL
!              SMC(I,J,KRT)=SMCAGGRT(KRT)/(AGGFACTRT**2)
               RT_DOMAIN(did)%SH2OX(I,J,KRT) = RT_DOMAIN(did)%SH2OAGGRT(KRT) &
                                             / (nlst_rt(did)%AGGFACTRT**2)

               RT_DOMAIN(did)%SICE(I,J,KRT)  = SICEAGGRT(KRT)/(nlst_rt(did)%AGGFACTRT**2)

               RT_DOMAIN(did)%SMC(I,J,KRT)   = RT_DOMAIN(did)%SH2OX(I,J,KRT) + &  ! Qian, 2016.12
                                               RT_DOMAIN(did)%SICE(I,J, KRT)
            end do

! Qian, add to debug SMC, SICE, SH2OX
            do KRT=1,nlst_rt(did)%NSOIL
               IF (RT_DOMAIN(did)%SH2OX(I,J,KRT) .LE. neta) THEN
                  print *, "Error negative SH2O...", &
                            RT_DOMAIN(did)%SH2OX(I,J,KRT),I,J,KRT
                  call hydro_stop("aggregateDomain")
               END IF

               IF (RT_DOMAIN(did)%SICE(I,J,KRT) .LE. neta) THEN
                  print *, "Error negative SICE...", &
                           RT_DOMAIN(did)%SICE(I,J,KRT),I,J,KRT
                  call hydro_stop("aggregateDomain")
               END IF

               IF (RT_DOMAIN(did)%SICE(I,J, KRT) .GT. &
                   RT_DOMAIN(did)%SMCMAX1(I,J) + peta) THEN
                  print*, "SICE exceeds porosity on aggregation...", &
                           RT_DOMAIN(did)%SICE(I,J, KRT), &
                           RT_DOMAIN(did)%SMCMAX1(I,J), I, J, KRT
                  call hydro_stop("aggregateDomain")
               ENDIF

               IF (RT_DOMAIN(did)%SMC(I,J, KRT) .GT. &
                   RT_DOMAIN(did)%SMCMAX1(I,J) + peta) THEN
                  print*, "SMC exceeds porosity on aggregation...", &
                           RT_DOMAIN(did)%SMC(I,J, KRT), &
                           RT_DOMAIN(did)%SH2OX(I,J,KRT), &
                           RT_DOMAIN(did)%SICE(I,J,KRT), &
                           RT_DOMAIN(did)%SMCMAX1(I,J), I, J, KRT
                  call hydro_stop("aggregateDomain")
               ENDIF
            end do
! Qian, end       

!DJG Calculate subgrid weighting array...

              do AGGFACYRT=nlst_rt(did)%AGGFACTRT-1,0,-1
                do AGGFACXRT=nlst_rt(did)%AGGFACTRT-1,0,-1
                  IXXRT=I*nlst_rt(did)%AGGFACTRT-AGGFACXRT
                  JYYRT=J*nlst_rt(did)%AGGFACTRT-AGGFACYRT
#ifdef MPP_LAND
       if(left_id.ge.0) IXXRT=IXXRT+1
       if(down_id.ge.0) JYYRT=JYYRT+1
#else
!yw ???
!       IXXRT=IXXRT+1
!       JYYRT=JYYRT+1
#endif
                  if (RT_DOMAIN(did)%LSMVOL .gt. 0.) then
                    RT_DOMAIN(did)%INFXSWGT(IXXRT,JYYRT) &
                                          = RT_DOMAIN(did)%SFCHEADSUBRT(IXXRT,JYYRT) &
                                          * RT_DOMAIN(did)%dist(IXXRT,JYYRT,9) &
					  / RT_DOMAIN(did)%LSMVOL
                  else
                    RT_DOMAIN(did)%INFXSWGT(IXXRT,JYYRT) &
                                          = 1./FLOAT(nlst_rt(did)%AGGFACTRT**2)
                  end if

                  do KRT=1,nlst_rt(did)%NSOIL

!!!yw added for debug
                   if(RT_DOMAIN(did)%SMCRT(IXXRT,JYYRT,KRT) .lt. neta) then
                      print*, "Error negative SMCRT", RT_DOMAIN(did)%SH2OWGT(IXXRT,JYYRT,KRT), RT_DOMAIN(did)%SMCRT(IXXRT,JYYRT,KRT),RT_DOMAIN(did)%SH2OX(I,J,KRT)
                      call hydro_stop("aggregateDomain")
                   endif
                   if(RT_DOMAIN(did)%SH2OWGT(IXXRT,JYYRT,KRT) .lt. neta) then
                      print *, "Error negative SH2OWGT", RT_DOMAIN(did)%SH2OWGT(IXXRT,JYYRT,KRT), RT_DOMAIN(did)%SMCRT(IXXRT,JYYRT,KRT),RT_DOMAIN(did)%SH2OX(I,J,KRT)
                      call hydro_stop("aggregateDomain")
                   endif
!end

                   ! Qian for debug SICE
                   if(RT_DOMAIN(did)%SICERT(IXXRT,JYYRT,KRT) .gt. RT_DOMAIN(did)%SMCMAX1(I,J) + peta) then
                      print*, "SICERT exceeds porosity on aggregation...", &
                               RT_DOMAIN(did)%SICERT(IXXRT,JYYRT,KRT), & 
                               RT_DOMAIN(did)%SMCMAX1(I,J), IXXRT, JYYRT, KRT
                      call hydro_stop("aggregateDomain")
                   endif
                   ! end


                    IF (RT_DOMAIN(did)%SMCRT(IXXRT,JYYRT,KRT) .GT. &
                        RT_DOMAIN(did)%SMCMAX1(I,J) + peta) THEN
!                        RT_DOMAIN(did)%SMCMAXRT(IXXRT,JYYRT,KRT)) THEN
                      print *, "SMCRT exceeded porosity upon aggregation...", &
                           RT_DOMAIN(did)%SMCRT(IXXRT,JYYRT,KRT),  &
                           RT_DOMAIN(did)%SMCMAX1(I,J), IXXRT, JYYRT, KRT
                      call hydro_stop("aggregateDomain")
                    END IF

                    IF (RT_DOMAIN(did)%SH2OX(I,J,KRT) .LE. neta) THEN
                        print *, "Erroneous value of SH2O...", &
                                 RT_DOMAIN(did)%SH2OX(I,J,KRT),I,J,KRT
                        print *, "Error negative SH2OX", RT_DOMAIN(did)%SH2OWGT(IXXRT,JYYRT,KRT), RT_DOMAIN(did)%SMCRT(IXXRT,JYYRT,KRT),RT_DOMAIN(did)%SH2OX(I,J,KRT)
                        call hydro_stop("aggregateDomain")
                    END IF

                    RT_DOMAIN(did)%SH2OWGT(IXXRT,JYYRT,KRT) &
                                 = RT_DOMAIN(did)%SMCRT(IXXRT,JYYRT,KRT) &
                                 / RT_DOMAIN(did)%SH2OX(I,J,KRT)

!?yw
                    RT_DOMAIN(did)%SH2OWGT(IXXRT,JYYRT,KRT) = max(1.0E-30, RT_DOMAIN(did)%SH2OWGT(IXXRT,JYYRT,KRT))
                  end do

                end do
              end do

         end do
        end do

        
#ifdef MPP_LAND
        call MPP_LAND_COM_REAL(RT_DOMAIN(did)%INFXSWGT, &
                               RT_DOMAIN(did)%IXRT,    &
                               RT_DOMAIN(did)%JXRT, 99)
	
        do i = 1, nlst_rt(did)%NSOIL
           call MPP_LAND_COM_REAL(RT_DOMAIN(did)%SH2OWGT(:,:,i), &
                                  RT_DOMAIN(did)%IXRT, &
				  RT_DOMAIN(did)%JXRT, 99)
        end do
#endif

!DJG Update SMC with SICE (unchanged) and new value of SH2O from routing...
!        RT_DOMAIN(did)%SMC = RT_DOMAIN(did)%SH2OX + RT_DOMAIN(did)%SICE

#ifdef HYDRO_D
        write(6,*), "*** Finish Aggregation"
#endif

	
      end subroutine aggregateDomain

! ---------------------------------------------------------
      subroutine HYDRO_ini(ntime, did, ix0, jx0, vegtyp,soltyp)
        implicit none
        integer :: ntime,   did
        integer :: rst_out, ix, jx
!        integer, OPTIONAL:: ix0,jx0
        integer :: ix0, jx0, i, j
        integer, dimension(ix0,jx0), OPTIONAL :: vegtyp, soltyp
        integer, allocatable :: gltype(:,:)

#ifdef MPP_LAND
     call  MPP_LAND_INIT()
#endif

     ! read the namelist
     ! the lsm namelist will be read by rtland sequentially again.
     call read_rt_nlst(nlst_rt(did) )


     IF (      nlst_rt(did)%GWBASESWCRT .eq. 0 &
         .and. nlst_rt(did)%SUBRTSWCRT  .eq. 0 &
         .and. nlst_rt(did)%OVRTSWCRT   .eq. 0 &
         .and. nlst_rt(did)%OVCHRTSWCRT .eq. 0 &  ! Qian, 2016.06
         .and. nlst_rt(did)%INTERGWSFW  .eq. 0 &  ! Qian, 2016.10
        ) return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! get the dimension 
     call get_file_dimension(trim(nlst_rt(did)%geo_static_flnm), ix, jx)
 
#ifdef MPP_LAND

     if (nlst_rt(did)%sys_cpl .eq. 1 .or. nlst_rt(did)%sys_cpl .eq. 4) then
     ! sys_cpl: 1-- coupling with HRLDAS but running offline lsm; 
     !          2-- coupling with WRF but do not run offline lsm  
     !          3-- coupling with LIS and do not run offline lsm  
     !          4-- coupling with CLM

     ! create 2 dimensiaon logical mapping of the CPUs for coupling with CLM or HRLDAS.
         call log_map2d()

         global_nx = ix  ! get from land model
         global_ny = jx  ! get from land model

         call mpp_land_bcast_int1(global_nx)
         call mpp_land_bcast_int1(global_ny)

!!! temp set global_nx to ix 
         rt_domain(did)%ix = global_nx
         rt_domain(did)%jx = global_ny

         ! over write the ix and jx
         call MPP_LAND_PAR_INI(1,rt_domain(did)%ix,rt_domain(did)%jx,&
              nlst_rt(did)%AGGFACTRT)

     else  
         !  coupled with WRF, LIS
         numprocs = node_info(1,1)

         call wrf_LAND_set_INIT(node_info,numprocs,nlst_rt(did)%AGGFACTRT)


         rt_domain(did)%ix = local_nx
         rt_domain(did)%jx = local_ny
     endif

     rt_domain(did)%g_IXRT = global_rt_nx
     rt_domain(did)%g_JXRT = global_rt_ny
     rt_domain(did)%ixrt   = local_rt_nx
     rt_domain(did)%jxrt   = local_rt_ny

#ifdef HYDRO_D
     write(6,*) "rt_domain(did)%g_IXRT, rt_domain(did)%g_JXRT, rt_domain(did)%ixrt, rt_domain(did)%jxrt"
     write(6,*)  rt_domain(did)%g_IXRT, rt_domain(did)%g_JXRT, rt_domain(did)%ixrt, rt_domain(did)%jxrt
     write(6,*) "rt_domain(did)%ix, rt_domain(did)%jx "
     write(6,*)  rt_domain(did)%ix, rt_domain(did)%jx 
     write(6,*) "global_nx, global_ny, local_nx, local_ny"
     write(6,*)  global_nx, global_ny, local_nx, local_ny
#endif

#else
     ! sequential
     rt_domain(did)%ix   = ix
     rt_domain(did)%jx   = jx
     rt_domain(did)%ixrt = ix*nlst_rt(did)%AGGFACTRT 
     rt_domain(did)%jxrt = jx*nlst_rt(did)%AGGFACTRT
#ifdef HYDRO_D
     write(6,*) "rt_domain(did)%ix, rt_domain(did)%jx "
     write(6,*)  rt_domain(did)%ix, rt_domain(did)%jx
#endif
#endif

     if ((nlst_rt(did)%OVCHRTSWCRT .eq. 1) .or. &
         (nlst_rt(did)%INTERGWSFW  .eq. 1)    ) then
        ! allocate RT_DOMAIN & rt2d arrays
        call rt2d_allocate(did,&
                           rt_domain(did)%ix,&
                           rt_domain(did)%jx,&
                           rt_domain(did)%ixrt,&
                           rt_domain(did)%jxrt,&
                           nlst_rt(did)%nsoil, &
                           nlst_rt(did)%NumFracts) ! Qian, 2016.06
#ifdef HYDRO_D
        write(6,*) "finish rt2d_allocate"
#endif

     else
        ! allocate RT_DOMAIN arrays
        call getChanDim(did)
#ifdef HYDRO_D
        write(6,*) "finish getChanDim "
#endif
     endif

     ! allocate gw2d arrays
     if (nlst_rt(did)%GWBASESWCRT .eq. 3 ) then
        call gw2d_allocate(did,&
                           rt_domain(did)%ixrt,&
                           rt_domain(did)%jxrt,&
                           nlst_rt(did)%nsoil)
#ifdef HYDRO_D
        write(6,*) "finish gw2d_allocate"
#endif
     endif

     ! calculate the distance between grids for routing.
     ! decompose the land parameter/data 

!    ix0= rt_domain(did)%ix
!    jx0= rt_domain(did)%jx
    
     if(present(vegtyp) .and. present(soltyp)) then
          call lsm_input(did,ix0=ix0,jx0=jx0,vegtyp0=vegtyp,soltyp0=soltyp)
     else
          call lsm_input(did,ix0=ix0,jx0=jx0)
     endif

#ifdef HYDRO_D
     write(6,*) "finish decomposion"
#endif

     call get_dist_lsm(did) 
     call get_dist_lrt(did)

     ! soil properity disaggregate
     call properity_disaggregate(did)  ! Qian, 2017.02

     ! rt model initilization
     call LandRT_ini(did)
#ifdef HYDRO_D
     write(6,*) "finish LandRT_ini"    
#endif


     ! if nlst_rt(did)%rst_bi_in == 1 
     ! replace some initial condition using restart file 
     if (nlst_rt(did)%rst_bi_in .eq. 1)then 
        call HYDRO_rst_in(did)
        print*, 'finish readin restart file'
     endif


     ! initilization for 2D groundwater routing  
     if (nlst_rt(did)%GWBASESWCRT .eq. 3) then
        call gw2d_ini(did,&
                      nlst_rt(did)%dtrt,&
                      nlst_rt(did)%dxrt0)
#ifdef HYDRO_D                        
        write(6,*) "finish gw2d_ini"      
#endif
     endif

     ! initialization for 2D Overland Channel Routing
     if ((nlst_rt(did)%OVCHRTSWCRT .eq. 1) .or. &
         (nlst_rt(did)%INTERGWSFW  .eq. 1)     )then
        call rt2d_ini(did,&               ! Qian, 2016.06
                      nlst_rt(did)%dtrt,&
                      nlst_rt(did)%dxrt0)
#ifdef HYDRO_D                        
        write(6,*) "finish rt2d_ini"
#endif

        ! cross-section initialization
        if (nlst_rt(did)%ntran .GT. 0) then
#ifdef MPP_LAND
           if (my_id == IO_id) then
              allocate(gltype(global_rt_nx,global_rt_ny))
           else
              allocate(gltype(1,1))
           endif
           call write_IO_RT_int(rt2d(did)%ltype,gltype)

           if (my_id == IO_id) then
              call transect_ini(did, global_rt_nx, global_rt_ny, gltype)
#ifdef HYDRO_D                        
              write(6,*) "finish transect_ini"
#endif
           endif
           deallocate(gltype)  
#else
           call transect_ini(did, rt2d(did)%ix, rt2d(did)%jx, rt2d(did)%ltype)
#ifdef HYDRO_D                        
              write(6,*) "finish transect_ini"
#endif
#endif
        endif
     endif
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! Perform ter rain adjustment of incoming solar
     IF (nlst_rt(did)%TERADJ_SOLAR.EQ.1 .and. nlst_rt(did)%CHANRTSWCRT.NE.2) THEN
#ifdef MPP_LAND
        call MPP_seq_land_SO8(rt_domain(did)%SO8LD_D, rt_domain(did)%SO8LD_Vmax,&
                              rt_domain(did)%TERRAIN, rt_domain(did)%dist_lsm,&
                              rt_domain(did)%ix,      rt_domain(did)%jx,&
                              global_nx,              global_ny)
#else
        call seq_land_SO8(rt_domain(did)%SO8LD_D,rt_domain(did)%SO8LD_Vmax,&
                          rt_domain(did)%TERRAIN,rt_domain(did)%dist_lsm,&
                          rt_domain(did)%ix,     rt_domain(did)%jx)
#endif
     ENDIF


     if (     (nlst_rt(did)%CHANRTSWCRT .eq. 1  &
         .or.  nlst_rt(did)%CHANRTSWCRT .eq. 2  &
         .or.  nlst_rt(did)%GWBASESWCRT .gt. 0) &
         .and. nlst_rt(did)%OVCHRTSWCRT .eq. 0  &      ! Qian, 2016.08
         .and. nlst_rt(did)%INTERGWSFW  .eq. 0  & 
        ) call get_basn_area(did)

     if (nlst_rt(did)%CHANRTSWCRT.EQ.1 .or. nlst_rt(did)%CHANRTSWCRT .eq. 2 ) then
!       call get_basn_area(did)
        call get_node_area(did)
     endif
     
!    if (trim(nlst_rt(did)%restart_file) == "") then
!         ! output at the initial time
!        call HYDRO_out(did)
!        return
!    endif

     ! restart the file

     ! jummp the initial time output
!    rt_domain(did)%out_counts = rt_domain(did)%out_counts + 1
!    rt_domain(did)%his_out_counts = rt_domain(did)%his_out_counts + 1

!     because gw2d(did)%ho, rt2d(did)%hlake, rt2d(did)%hbeg are needed 
!     in gw2d_ini, rt2d_ini
!     call HYDRO_rst_in(did)  ! comment by Qian, 2016.11

!    initialize these variables based on initialized hydro data 
!    RT_DOMAIN(did)%SMC(IX,JX,NSOIL), RT_DOMAIN(did)%SH2OX(IX,JX,NSOIL)
!    RT_DOMAIN(did)%ZWT(IX,JX), RT_DOMAIN(did)%SFCHEADRT(IX,JX) 
     call aggregateDomain(did)

     call HYDRO_out(did)

      end subroutine HYDRO_ini

! ---------------------------------------------------------
      subroutine lsm_input(did,ix0,jx0,vegtyp0,soltyp0)
         implicit none

         integer :: did
         integer :: ix0,jx0
         integer, dimension(ix0,jx0),OPTIONAL :: vegtyp0, soltyp0

         integer, parameter    :: leng = 100
         real, dimension(leng) :: xdum1, MAXSMC,refsmc,wltsmc
         real, dimension(leng) :: SATSUC, indexb                     ! Qian, 2016.06
         integer               :: i, j, nn
         integer, allocatable, dimension(:,:) :: soltyp

         allocate(soltyp(RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx) )
         soltyp   = 0

         call get2d_lsm_soltyp(soltyp,RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx,trim(nlst_rt(did)%geo_static_flnm))
         call get2d_lsm_real("HGT",RT_DOMAIN(did)%TERRAIN,RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx,trim(nlst_rt(did)%geo_static_flnm))
         call get2d_lsm_real("XLAT",RT_DOMAIN(did)%lat_lsm,RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx,trim(nlst_rt(did)%geo_static_flnm))
         call get2d_lsm_real("XLONG",RT_DOMAIN(did)%lon_lsm,RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx,trim(nlst_rt(did)%geo_static_flnm))
         call get2d_lsm_vegtyp(RT_DOMAIN(did)%VEGTYP,RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx,trim(nlst_rt(did)%geo_static_flnm))

!         ! C.Jiang 2020
!         call get2d_lsm_real("SLOPERT",RT_DOMAIN(did)%slope,RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx,trim(nlst_rt(did)%HRLDAS_CONSTANTS_FILE))
!         call get2d_lsm_real("LAI",RT_DOMAIN(did)%LAI,RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx,trim(nlst_rt(did)%HRLDAS_CONSTANTS_FILE))
!         call get2d_lsm_real("K1_FACTOR",RT_DOMAIN(did)%K1_FACTOR,RT_DOMAIN(did)%ix,RT_DOMAIN(did)%jx,trim(nlst_rt(did)%HRLDAS_CONSTANTS_FILE))

         ! coupling with WRF
!         if(nlst_rt(did)%sys_cpl .eq. 2 ) then
           if(present(soltyp0) .and. present(vegtyp0)) then

              where(soltyp0 == 14)  VEGTYP0 = 16
!qxia              where(VEGTYP0 == 16 ) soltyp0 = 14

              soltyp = soltyp0
              RT_DOMAIN(did)%VEGTYP = VEGTYP0

              ! coupled with CLM4, soil type does not connect with vegetable type
              ! because soil moisture in CLM4 depends on soil type
!              where(soltyp == 14) RT_DOMAIN(did)%VEGTYP = 16
!              where(RT_DOMAIN(did)%VEGTYP == 16 ) soltyp = 14

           endif
!         endif

         rt_domain(did)%soiltyp(1:RT_DOMAIN(did)%ix,1:RT_DOMAIN(did)%jx) = &
         soltyp(1:RT_DOMAIN(did)%ix,1:RT_DOMAIN(did)%jx)
#ifdef HYDRO_D
         write(6,*) 'finish reading soil type and vegetation type'
#endif

!        temporary set
!        RT_DOMAIN(did)%SMCRTCHK = 0   ! Qian, 2016.06
         RT_DOMAIN(did)%SMCAGGRT = 0
         RT_DOMAIN(did)%STCAGGRT = 0
         RT_DOMAIN(did)%SH2OAGGRT = 0

         RT_DOMAIN(did)%zsoil(1:nlst_rt(did)%nsoil) = nlst_rt(did)%zsoil8(1:nlst_rt(did)%nsoil)
!Qian         RT_DOMAIN(did)%sldpth(1)   = abs( RT_DOMAIN(did)%zsoil(1))
!Qian         do i = 2, nlst_rt(did)%nsoil
!Qian           RT_DOMAIN(did)%sldpth(i) = RT_DOMAIN(did)%zsoil(i-1)-RT_DOMAIN(did)%zsoil(i)
!Qian         enddo
!Qian         RT_DOMAIN(did)%SOLDEPRT = -1.0*RT_DOMAIN(did)%ZSOIL(nlst_rt(did)%NSOIL)
         ! Qian
         RT_DOMAIN(did)%sldpth(1:nlst_rt(did)%nsoil) = nlst_rt(did)%dzsoil(1:nlst_rt(did)%nsoil)
         RT_DOMAIN(did)%SOLDEPRT = -1.0*RT_DOMAIN(did)%ZSOIL(nlst_rt(did)%NSOIL) + &
                                   nlst_rt(did)%dzsoil(nlst_rt(did)%nsoil)

#ifdef HYDRO_D
         write(6,*) 'finish assignment for RT_DOMAIN%ZSOIL'
#endif

         ! input OV_ROUGH from OVROUGH.TBL
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
             read(71,*) xdum1(i),  MAXSMC(i), refsmc(i), wltsmc(i), &
                        SATSUC(i), indexb(i)  !Qian,2016.06
           end do 
         close(71)

#ifdef MPP_LAND
         endif
         call mpp_land_bcast_real(leng,RT_DOMAIN(did)%OV_ROUGH)
         call mpp_land_bcast_real(leng,xdum1)
         call mpp_land_bcast_real(leng,MAXSMC)
         call mpp_land_bcast_real(leng,refsmc)
         call mpp_land_bcast_real(leng,wltsmc)
         call mpp_land_bcast_real(leng,SATSUC)  ! Qian,2016.06
         call mpp_land_bcast_real(leng,indexb)  ! Qian,2016.06
#endif

         rt_domain(did)%lksat = 0.0
         do j = 1, RT_DOMAIN(did)%jx
            do i = 1, RT_DOMAIN(did)%ix
!yw            rt_domain(did)%lksat(i,j) = xdum1(soltyp(i,j) ) * 1000.0
               rt_domain(did)%lksat(i,j)   = xdum1(soltyp(i,j))
               rt_domain(did)%SUCSAT1(i,j) = SATSUC(soltyp(i,j))   ! Qian,2016.06
               rt_domain(did)%BSW1(i,j)    = indexb(soltyp(i,j))   ! Qian,2016.06
               IF(rt_domain(did)%VEGTYP(i,j) == 1 ) THEN   ! urban
                 rt_domain(did)%SMCMAX1(i,j) = 0.45
                 rt_domain(did)%SMCREF1(i,j) = 0.42
                 rt_domain(did)%SMCWLT1(i,j) = 0.40
               ELSE
                 rt_domain(did)%SMCMAX1(i,j) = MAXSMC(soltyp(I,J))
                 rt_domain(did)%SMCREF1(i,j) = refsmc(soltyp(I,J))
                 rt_domain(did)%SMCWLT1(i,j) = wltsmc(soltyp(I,J))
               ENDIF
            end do
         end do
#ifdef HYDRO_D
         write(6,*) 'finish reading soil properties '
#endif
         deallocate(soltyp)

      end subroutine lsm_input

      subroutine properity_disaggregate(did)
        implicit none
     
        ! this subroutine is to initially disaggregate some properity
        ! properity(1:ix, 1:jx) to properity(1:ixrt, 1:jxrt) 
 
        INTEGER, INTENT(IN) :: did
        INTEGER :: i, j
        INTEGER :: AGGFACYRT, AGGFACXRT
        INTEGER :: IXXRT, JYYRT

#ifdef HYDRO_D
#ifdef MPP_LAND
        if (my_id .eq. IO_id) &
#endif
        print *, "*** Beginning soil properity Disaggregation..."
#endif

         ! soil properity
         do j = 1, RT_DOMAIN(did)%jx
           do i = 1, RT_DOMAIN(did)%ix

             do AGGFACYRT = nlst_rt(did)%AGGFACTRT-1, 0, -1
               do AGGFACXRT = nlst_rt(did)%AGGFACTRT-1, 0, -1

                 IXXRT = I*nlst_rt(did)%AGGFACTRT - AGGFACXRT
                 JYYRT = J*nlst_rt(did)%AGGFACTRT - AGGFACYRT

#ifdef MPP_LAND
                 if(left_id .ge. 0) IXXRT = IXXRT+1
                 if(down_id .ge. 0) JYYRT = JYYRT+1
#endif
                 RT_DOMAIN(did)%SMCSATRT(IXXRT,JYYRT) = rt_domain(did)%SMCMAX1(I,J)  ! porosity, m3/m3 
                 RT_DOMAIN(did)%SUCSATRT(IXXRT,JYYRT) = rt_domain(did)%SUCSAT1(I,J)  ! saturated suction, m
                 RT_DOMAIN(did)%BSWRT(IXXRT,JYYRT)    = rt_domain(did)%BSW1(I,J)     ! parameter index B
                 RT_DOMAIN(did)%HYCONDRT(IXXRT,JYYRT) = rt_domain(did)%LKSAT(I,J)    ! saturated hydraulic conductivity, m/s

               end do
             end do

           end do
         end do
 
#ifdef MPP_LAND
        call MPP_LAND_COM_REAL(RT_DOMAIN(did)%SMCSATRT,RT_DOMAIN(did)%IXRT,RT_DOMAIN(did)%JXRT,99)
        call MPP_LAND_COM_REAL(RT_DOMAIN(did)%SUCSATRT,RT_DOMAIN(did)%IXRT,RT_DOMAIN(did)%JXRT,99)
        call MPP_LAND_COM_REAL(RT_DOMAIN(did)%BSWRT,   RT_DOMAIN(did)%IXRT,RT_DOMAIN(did)%JXRT,99)
        call MPP_LAND_COM_REAL(RT_DOMAIN(did)%HYCONDRT,RT_DOMAIN(did)%IXRT,RT_DOMAIN(did)%JXRT,99)   
#endif

      end subroutine properity_disaggregate

end module module_HYDRO_drv

! ---------------------------------------------------------
      ! stop the job due to the fatal error.
      subroutine HYDRO_stop(msg)
#ifdef MPP_LAND
        use module_mpp_land
#endif
        character(len=*) :: msg
        integer :: ierr
#ifdef HYDRO_D
        write(6,*) "The job is stoped due to the fatal error. ", trim(msg)
        call flush(6)
#endif
#ifdef MPP_LAND
#ifndef HYDRO_D
        print*, "---"
        print*, "ERROR! Program stopped. Recompile with environment &
                 variable HYDRO_D set to 1 for enhanced debug information."
        print*, " "
#endif

!        call mpp_land_sync()
!        write(my_id+90,*) msg
!        call flush(my_id+90)

         call mpp_land_abort()
         call MPI_finalize(ierr)
#else
         stop "Fatal Error"
#endif

     return
     end  subroutine HYDRO_stop  

! ---------------------------------------------------------
     ! stop the job due to the fatal error.
     subroutine HYDRO_finish()
#ifdef MPP_LAND
        USE module_mpp_land
#endif
        integer :: ierr

        print*, "The model finished successfully !!!"
#ifdef MPP_LAND
!         call mpp_land_abort()
         call flush(6)
         call mpp_land_sync()
         call MPI_finalize(ierr)
         stop 
#else
         stop 
#endif

     return
     end  subroutine HYDRO_finish
