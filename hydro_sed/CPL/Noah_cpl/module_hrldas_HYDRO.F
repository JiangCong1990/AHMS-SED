module module_HRLDAS_HYDRO

! NDHMS  module   
#ifdef MPP_LAND
    use module_mpp_land, only: global_nx, global_ny, decompose_data_real, &
                 write_io_real, my_id, mpp_land_bcast_real1, IO_id, &
                 mpp_land_bcast_int1
#endif
    use module_HYDRO_drv, only: HYDRO_ini, HYDRO_exe, HYDRO_rst_out
    use module_rt_data, only:  rt_domain
    use module_namelist, only: nlst_rt
    use module_gw_gw2d_data, only: gw2d

    implicit none
    integer begg, endg
    integer :: numg, numl, numc, nump
    INTEGER, PARAMETER :: double=8
    real(kind=double), pointer :: r2p(:,:) , r1p(:)
     
    integer ::  begl, endl, begc, endc, begp, endp

    real, allocatable, dimension(:,:) :: vg_test
    integer :: nn
    integer :: open_unit_status



CONTAINS

    subroutine hrldas_cpl_HYDRO(STC,SMC,SH2OX,infxsrt,sfcheadrt,soldrain,ii,jj,kk, qsgw)


       implicit none

       integer ii,jj,kk
       integer k
       real,dimension(ii,jj,kk) :: STC,SMC,SH2OX
       real,dimension(ii,jj) ::infxsrt,sfcheadrt, soldrain, qsgw

        integer ::  did

        integer ntime

        integer :: i,j

        real*8 :: t1, t2, dact
        integer :: clock_count_1, clock_count_2, clock_rate 
        save dact


!output flux and state variable

        did = 1

        if(  nlst_rt(did)%SUBRTSWCRT  .eq. 0   .and. &
             nlst_rt(did)%OVRTSWCRT   .eq. 0   .and. &
             nlst_rt(did)%GWBASESWCRT .eq. 0   .and. &
             nlst_rt(did)%OVCHRTSWCRT .eq. 0   .and. &   ! Qian, 2016.08
             nlst_rt(did)%INTERGWSFW  .eq. 0         &
          )  return


! decompose the hrldas 1-d data into routing domain
#ifdef MPP_LAND
        do k = 1, kk
           call decompose_data_real(STC(:,:,k),RT_DOMAIN(did)%STC(:,:,k))
           call decompose_data_real(SMC(:,:,k),RT_DOMAIN(did)%SMC(:,:,k))
           call decompose_data_real(SH2OX(:,:,k),RT_DOMAIN(did)%SH2OX(:,:,k))
        end do
        call decompose_data_real(infxsrt,RT_DOMAIN(did)%infxsrt)
        call decompose_data_real(soldrain,RT_DOMAIN(did)%soldrain)
	
	if(nlst_rt(did)%GWBASESWCRT == 3)&
        call decompose_data_real(qsgw,gw2d(did)%qsgw)
#else
        RT_DOMAIN(did)%STC = STC
        RT_DOMAIN(did)%SMC = SMC
        RT_DOMAIN(did)%SH2OX = SH2OX
        RT_DOMAIN(did)%infxsrt = infxsrt
        RT_DOMAIN(did)%soldrain = soldrain
        if(nlst_rt(did)%GWBASESWCRT == 3) &
	  gw2d(did)%qsgw = qsgw
#endif
  
#ifdef MPP_LAND
#ifdef HYDRO_D
        if(my_id .eq. IO_id) then
           call system_clock(count=clock_count_1, count_rate=clock_rate)
        endif
#endif
#endif

        ntime = 1

        call HYDRO_exe(did)
#ifdef MPP_LAND
#ifdef HYDRO_D
           call mpp_land_sync()
        if(my_id .eq. IO_id) then
           call system_clock(count=clock_count_2, count_rate=clock_rate)
           dact = dact + float(clock_count_2-clock_count_1)/float(clock_rate)
           write(6,*) "timing: ", float(clock_count_2-clock_count_1)/float(clock_rate), &
                 "accumulated time (s): ",dact
        endif
#endif
#endif

! add for update the HRLDAS state variable.
#ifdef MPP_LAND
        do k = 1, kk
           call write_io_real(rt_domain(did)%STC(:,:,k),STC(:,:,k))
           call write_io_real(rt_domain(did)%SMC(:,:,k),SMC(:,:,k))
           call write_io_real(rt_domain(did)%SH2OX(:,:,k),SH2OX(:,:,k))
        end do
         call write_io_real(rt_domain(did)%sfcheadrt,sfcheadrt)
	if(nlst_rt(did)%GWBASESWCRT == 3)&
         call write_io_real(gw2d(did)%qsgw,qsgw)
#else
        STC = rt_domain(did)%STC
        SMC = rt_domain(did)%SMC
        SH2OX = rt_domain(did)%SH2OX
        sfcheadrt = rt_domain(did)%sfcheadrt
	if(nlst_rt(did)%GWBASESWCRT == 3)&
        qsgw = gw2d(did)%qsgw
#endif

!? not sure for the following
!           grid%xice(its:ite,jts:jte) = rt_domain(did)%sice


     return 
     end subroutine hrldas_cpl_HYDRO

    subroutine hrldas_cpl_HYDRO_ini(STC,SMC,SH2OX,infxsrt,sfcheadrt,soldrain,ii,jj,kk,kt,dt, olddate,zsoil)


       implicit none

       integer ii,jj,kk
       integer k, kt
       real :: dt
       real,dimension(ii,jj,kk) :: STC,SMC,SH2OX
       real,dimension(ii,jj) ::infxsrt,sfcheadrt, soldrain
       real, dimension(kk) :: zsoil
       character(len = *) :: olddate

        integer ::  did

        integer ntime

        integer :: i,j

!output flux and state variable

        did = 1


        if(.not. RT_DOMAIN(did)%initialized) then  
           nlst_rt(did)%dt = dt
           nlst_rt(did)%olddate(1:19) = olddate(1:19)
           nlst_rt(did)%startdate(1:19) = olddate(1:19)

           nlst_rt(did)%nsoil = kk
#ifdef MPP_LAND
           call mpp_land_bcast_int1(nlst_rt(did)%nsoil)
#endif
           allocate(nlst_rt(did)%zsoil8(nlst_rt(did)%nsoil))
           nlst_rt(did)%zsoil8(1:nlst_rt(did)%nsoil) = zsoil(1:nlst_rt(did)%nsoil)
           
           
           call HYDRO_ini(ntime,did,ix0=1,jx0=1)

            if(nlst_rt(did)%sys_cpl .ne. 1) then
               write(6,*) "Error: sys_cpl should be 1."
               call hydro_stop("hrldas_cpl_HYDRO_ini")
            endif

           RT_DOMAIN(did)%initialized = .true.
   
          if(  nlst_rt(did)%SUBRTSWCRT  .eq. 0   .and. &
               nlst_rt(did)%OVRTSWCRT   .eq. 0   .and. &
               nlst_rt(did)%GWBASESWCRT .eq. 0   .and. &
               nlst_rt(did)%OVCHRTSWCRT .eq. 0   .and. &   ! Qian, 2016.08
               nlst_rt(did)%INTERGWSFW  .eq. 0         &
            )  return


#ifdef MPP_LAND
           call mpp_land_bcast_real1(nlst_rt(did)%dt)
           call write_io_real(rt_domain(did)%sfcheadrt,sfcheadrt)
           call write_io_real(rt_domain(did)%infxsrt, infxsrt)
!          call write_io_real(rt_domain(did)%soldrain,soldrain)
           if(nlst_rt(did)%rst_typ .eq. 1) then
              do k = 1, kk
                 call write_io_real(rt_domain(did)%STC(:,:,k),STC(:,:,k))
                 call write_io_real(rt_domain(did)%SMC(:,:,k),SMC(:,:,k))
                 call write_io_real(rt_domain(did)%SH2OX(:,:,k),SH2OX(:,:,k))
              end do
           else
              if(nlst_rt(did)%sys_cpl .eq. 1) then
                 where( abs(STC)   .gt. 500) stc = 282
                 where( abs(SMC)   .gt. 500) SMC = 0.25
                 where( abs(SH2OX) .gt. 500) SH2OX = 0.25

                 where( abs(SH2OX) .eq. 1)   SH2OX = 0.25
                 where( abs(SMC)   .eq. 1)   SMC = 0.25
                 where( SH2OX      .le. 0.0001) SH2OX = 0.25
                 where( SMC        .le. 0.0001) SMC = 0.25
              endif
           endif 
#else
           sfcheadrt = rt_domain(did)%sfcheadrt
           infxsrt = rt_domain(did)%infxsrt
           if(nlst_rt(did)%rst_typ .eq. 1) then
              STC = rt_domain(did)%STC
              SMC = rt_domain(did)%SMC
              SH2OX = rt_domain(did)%SH2OX
           else
              if(nlst_rt(did)%sys_cpl .eq. 1) then
                 where( abs(STC) .gt. 500) stc = 282
                 where( abs(SMC) .gt. 500) SMC = 0.25
                 where( abs(SH2OX) .gt. 500) SH2OX = 0.25

                 where( abs(SH2OX) .eq. 1) SH2OX = 0.25
                 where( abs(SMC) .eq. 1) SMC = 0.25
              endif
           endif
#endif
        endif

     return 
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

     subroutine getNameList(varName,pVar)
          implicit none
          integer :: pVar
          integer :: did 
          character(len=*) :: varName
          did = 1
          if(varName .eq. "GWSPINCYCLES") then
               pVar = nlst_rt(did)%gwSpinCycles
          elseif(varName .eq. "GWBASESWCRT") then
               pVar = nlst_rt(did)%GWBASESWCRT
          elseif(varName .eq. "GWSOILCPL") then
               pVar = nlst_rt(did)%gwsoilcpl
          else
               pVar = 0
          endif
     end subroutine getNameList
end module module_HRLDAS_HYDRO
