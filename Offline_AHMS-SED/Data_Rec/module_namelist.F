Module module_namelist

#ifdef MPP_LAND
          USE module_mpp_land
#endif

    IMPLICIT NONE
    INTEGER, PARAMETER :: max_domain=5

#include "namelist.inc"
    TYPE(namelist_rt_field) , dimension(max_domain) :: nlst_rt
    save nlst_rt 

CONTAINS 

    subroutine read_rt_nlst(nlst)     
          implicit none

          TYPE(namelist_rt_field) nlst

          integer:: ierr
          integer:: RT_OPTION = 0
          integer:: CHANRTSWCRT = 0
          integer:: channel_option = 0
          integer:: SUBRTSWCRT = 0
          integer:: OVRTSWCRT = 0
          integer:: AGGFACTRT = 0
          integer:: GWBASESWCRT = 0
          integer:: GW_RESTART = 0                   
          integer:: RSTRT_SWC = 0
          integer:: TERADJ_SOLAR = 0 
          integer:: sys_cpl = 0
          integer:: rst_typ = 0
          integer:: rst_bi_in = 0
          integer:: rst_bi_out = 0
          integer:: gwChanCondSw = 0
          integer:: GwPreDiagInterval = 0
          integer:: gwsoilcpl = 0
          integer:: GwSpinCycles = 0
          integer:: GwPreCycles = 0 
          integer:: INTERGWSFW = 0            ! if activated interaction of groundwater and surface water, Qian, 2016.10  
          integer:: OVCHRTSWCRT = 0           ! if activated 2D Overland_Channel_Routing, Qian, 2016.06
          integer:: OVERCHANRT_METHOD = 0     ! numerical method for 2D Overland_Channel_Routing, Qian, 2016.06
          integer:: num_transect = 0          ! number of transects specified, Qian, 2016.06
          real   :: DT_2DTHRT = 0.            ! timestep interval of terrain hydrologic routing, Qian
          real   :: dt_flowtran = 0.          ! interval for streamflow output, Qian, 2016.10
          integer:: IRRIG_OPTION = 0          ! Irrigation module, C.Jiang, 2021.08
          integer:: indexeros = 0             ! Erosion and Sediment Routing, C.Jiang, 2022.01
          integer:: NumFracts = 0             ! Maximum number of size fractions of soil types, C.Jiang, 2022.01
          real   :: e1 = 0.                     !
          real   :: e2 = 0.                     ! 
          real   :: e3 = 0.                     ! water erosion
          real   :: alpha = 0.                  ! gully erosion 
          real   :: notillf = 0.                ! cropland erosion, C.Jiang,2023.01
          real   :: DTRT = 0.
          real   :: dxrt = 0.
          real   :: gwChanCondConstIn = 0.
          real   :: gwChanCondConstOut = 0.
          real   :: gwIhShift = 0.
          character(len=256) :: rt2d_flnm=""  ! input data for 2D Overland_Channel_Routing, Qian, 2016.06   
          character(len=256) :: route_topo_f=""
          character(len=256) :: route_chan_f=""
          character(len=256) :: route_link_f=""
          character(len=256) :: route_lake_f=""
          character(len=256) :: route_direction_f=""
          character(len=256) :: route_order_f=""
          character(len=256) :: gwbasmskfil =""
          character(len=256) :: gwstrmfil =""
          character(len=256) :: geo_finegrid_flnm =""
          character(len=256) :: RESTART_FILE = ""
          character(len=256) :: geo_static_flnm = ""
          character(len=256) :: RESTART_OUTNAME = "" ! Qian, 2016.11
          character(len=256) :: GW2D_OUTNAME = ""    ! Qian, 2016.11
          character(len=256) :: RT2D_OUTNAME = ""    ! Qian, 2016.11
          integer :: SOLVEG_INITSWC = 0
          real    :: out_dt = 0.
          real    :: rst_dt = 0.
          logical :: GwPreDiag = .False.
          logical :: GwSpinUp = .False.
          integer :: split_output_count = 0
          integer :: order_to_write = 0
          integer :: igrid = 0
          integer :: DEEPGWSPIN = 0
       
          integer :: CHRTOUT_DOMAIN = 0       ! Netcdf point timeseries output at all channel points
          integer :: CHRTOUT_GRID = 0         ! Netcdf grid of channel streamflow values
          integer :: LSMOUT_DOMAN = 0         ! Netcdf grid of variables passed between LSM and routing components
          integer :: RTOUT_DOMAIN = 0         ! Netcdf grid of terrain routing variables on routing grid
          integer :: output_gw = 0            ! switch to output gw variables
          integer :: outlake = 0              ! switch to output lake variables
          integer :: output_rt2d = 0          ! switch to output 2D-routing variables, Qian, 2016.10

          integer :: i

!!! add the following two dummy variables 

       namelist /HYDRO_nlist/ &
            RESTART_FILE,SPLIT_OUTPUT_COUNT,IGRID,&
            geo_static_flnm, &
            out_dt, rst_dt, &
            DEEPGWSPIN, SOLVEG_INITSWC, &
            RT_OPTION, CHANRTSWCRT, channel_option, &
            SUBRTSWCRT,OVRTSWCRT,AGGFACTRT, dtrt,dxrt,&
            GwSpinCycles, GwPreCycles, GwSpinUp, GwPreDiag, GwPreDiagInterval, gwIhShift, &
            GWBASESWCRT, gwChanCondSw, gwChanCondConstIn, gwChanCondConstOut , &
            route_topo_f,route_chan_f,route_link_f,route_lake_f, &
            route_direction_f,route_order_f,gwbasmskfil, geo_finegrid_flnm,&
            gwstrmfil,GW_RESTART,RSTRT_SWC,TERADJ_SOLAR, sys_cpl, &
            order_to_write , rst_typ, rst_bi_in, rst_bi_out, gwsoilcpl, &
            CHRTOUT_DOMAIN,CHRTOUT_GRID,LSMOUT_DOMAN,RTOUT_DOMAIN, output_gw, outlake, &
            output_rt2d, INTERGWSFW, OVCHRTSWCRT, OVERCHANRT_METHOD, num_transect, &
            DT_2DTHRT, &
            dt_flowtran, rt2d_flnm, RESTART_OUTNAME, GW2D_OUTNAME, RT2D_OUTNAME, IRRIG_OPTION, indexeros, NumFracts,e1,e2,e3,alpha,notillf  ! Qian, 2016.06
                    
#ifdef MPP_LAND
       if(IO_id .eq. my_id) then
#endif
          open(30, file="hydro.namelist", form="FORMATTED")
          read(30, HYDRO_nlist, iostat=ierr)
          close(30)
#ifdef MPP_LAND
       endif
#endif

!Qian   rst_bi_in = 0
!Qian   rst_bi_out = 0

#ifdef MPP_LAND
!  call mpp_land_bcast_real1(DT)
  call mpp_land_bcast_int1(SPLIT_OUTPUT_COUNT)
  call mpp_land_bcast_int1(IGRID)
  call mpp_land_bcast_real1(out_dt)
  call mpp_land_bcast_real1(rst_dt)
  call mpp_land_bcast_int1(DEEPGWSPIN)
  call mpp_land_bcast_int1(SOLVEG_INITSWC)
#endif


#ifdef MPP_LAND
      call mpp_land_bcast_int1(nlst%NSOIL)
      do i = 1, nlst%NSOIL
        call mpp_land_bcast_real1(nlst%ZSOIL8(i))
      end do
#ifdef HYDRO_D
      write(6,*) "nlst%NSOIL = ", nlst%NSOIL
      write(6,*) "nlst%ZSOIL8 = ",nlst%ZSOIL8
#endif
#endif

!  nlst%DT = DT
  nlst%RESTART_FILE       = RESTART_FILE
  nlst%SPLIT_OUTPUT_COUNT = SPLIT_OUTPUT_COUNT
  nlst%IGRID              = IGRID
  nlst%geo_static_flnm    = geo_static_flnm
  nlst%out_dt             = out_dt
  nlst%rst_dt             = rst_dt
  nlst%DEEPGWSPIN         = DEEPGWSPIN
  nlst%SOLVEG_INITSWC     = SOLVEG_INITSWC

! modification for DT_2DTHRT
  if (OVCHRTSWCRT .EQ. 1) then
     if (DT_2DTHRT .LE. 0.) then
        DT_2DTHRT = DTRT
     elseif (DT_2DTHRT .GT. DTRT) then
        DT_2DTHRT = DTRT
     endif
  endif

#ifdef MPP_LAND
  call mpp_land_bcast_char(256,nlst%RESTART_FILE)
#endif

  write(nlst%hgrid,'(I1)') igrid


  if(RESTART_FILE .eq. "") rst_typ = 0

#ifdef MPP_LAND
  !bcast namelist variable.
  call mpp_land_bcast_int1(rt_option)
  call mpp_land_bcast_int1(CHANRTSWCRT)
  call mpp_land_bcast_int1(channel_option)
  call mpp_land_bcast_int1(SUBRTSWCRT)
  call mpp_land_bcast_int1(OVRTSWCRT)
  call mpp_land_bcast_int1(INTERGWSFW)        ! Qian
  call mpp_land_bcast_int1(OVCHRTSWCRT)       ! Qian
  call mpp_land_bcast_int1(OVERCHANRT_METHOD) ! Qian
  call mpp_land_bcast_int1(num_transect)      ! Qian
  call mpp_land_bcast_int1(AGGFACTRT)
  call mpp_land_bcast_real1(DTRT)
  call mpp_land_bcast_real1(DXRT)
  call mpp_land_bcast_real1(DT_2DTHRT)        ! Qian
  call mpp_land_bcast_real1(dt_flowtran)      ! Qian
  call mpp_land_bcast_real1(gwChanCondConstIn)
  call mpp_land_bcast_real1(gwChanCondConstOut)
  call mpp_land_bcast_real1(gwIhShift)
  call mpp_land_bcast_int1(GWBASESWCRT)
  call mpp_land_bcast_int1(GWSOILCPL)
  call mpp_land_bcast_int1(gwChanCondSw)
  call mpp_land_bcast_int1(GwSpinCycles)
  call mpp_land_bcast_int1(GwPreCycles)
  call mpp_land_bcast_log1(GwPreDiag)
  call mpp_land_bcast_log1(GwSpinUp)
  call mpp_land_bcast_int1(GwPreDiagInterval)
  call mpp_land_bcast_int1(GW_RESTART)
  call mpp_land_bcast_int1(RSTRT_SWC  )
  call mpp_land_bcast_int1(TERADJ_SOLAR)
  call mpp_land_bcast_int1(sys_cpl)
  call mpp_land_bcast_int1(rst_typ)
  call mpp_land_bcast_int1(rst_bi_in)
  call mpp_land_bcast_int1(rst_bi_out)
  call mpp_land_bcast_int1(order_to_write)
  call mpp_land_bcast_int1(CHRTOUT_DOMAIN)
  call mpp_land_bcast_int1(output_gw)
  call mpp_land_bcast_int1(outlake)
  call mpp_land_bcast_int1(output_rt2d)       ! Qian
  call mpp_land_bcast_int1(CHRTOUT_GRID)
  call mpp_land_bcast_int1(LSMOUT_DOMAN)
  call mpp_land_bcast_int1(RTOUT_DOMAIN)
#endif
    nlst%CHRTOUT_DOMAIN = CHRTOUT_DOMAIN
    nlst%output_gw      = output_gw
    nlst%outlake      = outlake
    nlst%output_rt2d  = output_rt2d
    nlst%CHRTOUT_GRID = CHRTOUT_GRID
    nlst%LSMOUT_DOMAN = LSMOUT_DOMAN
    nlst%RTOUT_DOMAIN = RTOUT_DOMAIN
    nlst%RT_OPTION = RT_OPTION
    nlst%CHANRTSWCRT = CHANRTSWCRT
    nlst%GW_RESTART  = GW_RESTART 
    nlst%RSTRT_SWC   = RSTRT_SWC  
    nlst%channel_option = channel_option
    nlst%DTRT   = DTRT
    nlst%DTCT   = DTRT
    nlst%SUBRTSWCRT = SUBRTSWCRT
    nlst%OVRTSWCRT = OVRTSWCRT
    nlst%dxrt0 = dxrt
    nlst%AGGFACTRT = AGGFACTRT
    nlst%GWBASESWCRT = GWBASESWCRT
    nlst%GWSOILCPL= GWSOILCPL
    nlst%gwChanCondSw = gwChanCondSw
    nlst%gwChanCondConstIn = gwChanCondConstIn
    nlst%gwChanCondConstOut = gwChanCondConstOut
    nlst%gwIhShift = gwIhShift
    nlst%GwSpinCycles = GwSpinCycles
    nlst%GwPreCycles = GwPreCycles
    nlst%GwPreDiag = GwPreDiag
    nlst%GwSpinUp = GwSpinUp
    nlst%GwPreDiagInterval = GwPreDiagInterval
    nlst%TERADJ_SOLAR = TERADJ_SOLAR
    nlst%sys_cpl = sys_cpl
    nlst%rst_typ = rst_typ
    nlst%rst_bi_in = rst_bi_in
    nlst%rst_bi_out = rst_bi_out
    nlst%order_to_write = order_to_write

    nlst%INTERGWSFW        = INTERGWSFW         ! Qian, 2016.10
    nlst%OVCHRTSWCRT       = OVCHRTSWCRT        ! Qian, 2016.06
    nlst%OVERCHANRT_METHOD = OVERCHANRT_METHOD  ! Qian, 2016.06
    nlst%ntran             = num_transect       ! Qian, 2016.06
    nlst%DT_2DTHRT         = DT_2DTHRT          ! Qian
    nlst%dt_flowtran       = dt_flowtran        ! Qian, 2016.10
    nlst%IRRIG_OPTION      = IRRIG_OPTION       ! C.Jiang,2022.01
    nlst%indexeros         = indexeros          ! C.Jiang,2022.01
    nlst%NumFracts         = NumFracts          ! C.Jiang,2022.01
    nlst%e1                = e1
    nlst%e2                = e2
    nlst%e3                = e3
    nlst%alpha             = alpha
    nlst%notillf           = notillf  
 
! files
    nlst%rt2d_flnm       = rt2d_flnm            ! Qian, 2016.06
    nlst%RESTART_OUTNAME = RESTART_OUTNAME      ! Qian, 2016.11
    nlst%GW2D_OUTNAME    = GW2D_OUTNAME         ! Qian, 2016.11
    nlst%RT2D_OUTNAME    = RT2D_OUTNAME         ! Qian, 2016.11
    nlst%route_topo_f =  route_topo_f
    nlst%route_chan_f = route_chan_f 
    nlst%route_link_f = route_link_f
    nlst%route_lake_f =route_lake_f
    nlst%route_direction_f =  route_direction_f
    nlst%route_order_f =  route_order_f
    nlst%gwbasmskfil =  gwbasmskfil
    nlst%gwstrmfil =  gwstrmfil
    nlst%geo_finegrid_flnm =  geo_finegrid_flnm

#ifdef MPP_LAND
  if(my_id .eq. IO_id) then
#endif
#ifdef HYDRO_D
    write(6,*) " ---------------------------------------" 
    write(6,*) " output of the namelist file "

    write(6,*) " nlst%RT_OPTION ", RT_OPTION
    write(6,*) " nlst%CHANRTSWCRT ", CHANRTSWCRT
    write(6,*) " nlst%GW_RESTART  ", GW_RESTART 
    write(6,*) " nlst%RSTRT_SWC   ", RSTRT_SWC  
    write(6,*) " nlst%channel_option ", channel_option
    write(6,*) " nlst%DTRT   ", DTRT
    write(6,*) " nlst%SUBRTSWCRT ", SUBRTSWCRT
    write(6,*) " nlst%OVRTSWCRT ", OVRTSWCRT
    write(6,*) " nlst%INTERGWSFW ", INTERGWSFW                ! Qian, 2016.10
    write(6,*) " nlst%OVERCHANRT_SWCRT ", OVCHRTSWCRT         ! Qian, 2016.06
    write(6,*) " nlst%OVERCHANRT_METHOD ", OVERCHANRT_METHOD  ! Qian, 2016.06
    write(6,*) " nlst%ntran ", num_transect                   ! Qian, 2016.06
    write(6,*) " nlst%DT_2DTHRT ", DT_2DTHRT                  ! Qian
    write(6,*) " nlst%dt_flowtran ", dt_flowtran              ! Qian, 2016.10
    write(6,*) " nlst%IRRIG_OPTION ", IRRIG_OPTION            ! cjiang
    write(6,*) " nlst%indexeros ", indexeros
    write(6,*) " nlst%NumFracts ", NumFracts
    write(6,*) " nlst%e1",e1
    write(6,*) " nlst%e2",e2
    write(6,*) " nlst%e3",e3
    write(6,*) " nlst%alpha", alpha
    write(6,*) " nlst%notillf", notillf
    write(6,*) " nlst%dxrt0 ", dxrt
    write(6,*) " nlst%AGGFACTRT ", AGGFACTRT
    write(6,*) " nlst%GWBASESWCRT ", GWBASESWCRT
    write(6,*) " nlst%GWSOILCPL ", GWSOILCPL
    write(6,*) " nlst%gwChanCondSw ", gwChanCondSw
    write(6,*) " nlst%gwChanCondConstIn ", gwChanCondConstIn
    write(6,*) " nlst%gwChanCondConstOut ", gwChanCondConstOut
    write(6,*) " nlst%gwIhShift ", gwIhShift
    write(6,*) " nlst%GwSpinCycles ", GwSpinCycles
    write(6,*) " nlst%GwPreDiag ", GwPreDiag
    write(6,*) " nlst%GwPreDiagInterval ", GwPreDiagInterval
    write(6,*) " nlst%TERADJ_SOLAR ", TERADJ_SOLAR
    write(6,*) " nlst%sys_cpl ", sys_cpl
    write(6,*) " nlst%rst_typ ", rst_typ
    write(6,*) " nlst%order_to_write ", order_to_write
    write(6,*) " nlst%route_topo_f   ",  route_topo_f
    write(6,*) " nlst%route_chan_f ", route_chan_f 
    write(6,*) " nlst%route_link_f ", route_link_f
    write(6,*) " nlst%route_lake_f ",route_lake_f
    write(6,*) " nlst%route_direction_f ",  route_direction_f
    write(6,*) " nlst%route_order_f ",  route_order_f
    write(6,*) " nlst%gwbasmskfil ",  gwbasmskfil
    write(6,*) " nlst%gwstrmfil ",  gwstrmfil
    write(6,*) " nlst%geo_finegrid_flnm ",  geo_finegrid_flnm
    write(6,*) " nlst%rt2d_flnm ", rt2d_flnm                  ! Qian, 2016.06
    write(6,*) " nlst%RESTART_OUTNAME ", RESTART_OUTNAME      ! Qian, 2016.11
    write(6,*) " nlst%GW2D_OUTNAME ", GW2D_OUTNAME            ! Qian, 2016.11
    write(6,*) " nlst%RT2D_OUTNAME ", RT2D_OUTNAME            ! Qian, 2016.11

#endif
#ifdef MPP_LAND
  endif
#endif

#ifdef MPP_LAND
  !bcast other  variable.
      call mpp_land_bcast_real1(nlst%dt)
#endif
      return
    end subroutine read_rt_nlst


end module module_namelist
