   TYPE namelist_rt_field  
      
       integer                        :: nsoil, SOLVEG_INITSWC
       real, allocatable,dimension(:) :: ZSOIL8
       real, allocatable,dimension(:) :: DZSOIL
       real                           :: out_dt, rst_dt, dt
       integer                        :: START_YEAR
       integer                        :: START_MONTH
       integer                        :: START_DAY
       integer                        :: START_HOUR
       integer                        :: START_MIN
       character(len=256)             :: restart_file = ""
       integer                        :: split_output_count
       integer                        :: igrid
       integer                        :: rst_bi_in    ! used for parallel io with large restart file.
       integer                        :: rst_bi_out   ! used for parallel io with large restart file.
                                                      ! each process will output the restart tile.
       character(len=256)             :: geo_static_flnm = ""
!       character(len=256)             :: HRLDAS_CONSTANTS_FILE = ""
       integer                        :: DEEPGWSPIN
       integer                        :: order_to_write, rst_typ
       
!      additional character
       character         :: hgrid
       character(len=19) :: olddate="123456"
       character(len=19) :: startdate="123456"
       character(len=19) :: sincedate="123456"

          integer:: RT_OPTION,         CHANRTSWCRT,  channel_option,            &
                    SUBRTSWCRT,        OVRTSWCRT,    AGGFACTRT,                 &
                    GWBASESWCRT,       GW_RESTART,   RSTRT_SWC,   TERADJ_SOLAR, &
                    sys_cpl,           gwChanCondSw, GwPreCycles, GwSpinCycles, &
                    GwPreDiagInterval, gwsoilcpl

          integer :: INTERGWSFW          ! if activated interaction of groundwater and surface water, Qian, 2016.10
          integer :: OVCHRTSWCRT         ! if activated 2D Overland_Channel_Routing, Qian, 2016.06      
          integer :: OVERCHANRT_METHOD   ! numerical method for 2D Overland_Channel_Routing, Qian. 2016.06
          real    :: DT_2DTHRT           ! timestep interval of 2D Terrain Hydrologic Routing
          integer :: ntran               ! number of transect, Qian, 2016.06
          real    :: dt_flowtran         ! interval for averaging-streamflow of specified transects (sec.)
          integer :: IRRIG_OPTION        ! if activated irrigation module = 1, else = 0
          integer :: indexeros           ! Erosion option index (1: erosion is simulated; not equal 1: otherwise)
          integer :: NumFracts           ! Maximum number of size fractions found any of the soil types defined in the soil type input grid
          real    :: e1,e2,e3            ! water erosion coefficient
          real    :: alpha, notillf      ! gully and crop coefficients

          logical :: GwPreDiag, GwSpinUp
          real    :: DTRT, dxrt0, DTCT, gwChanCondConstIn, gwChanCondConstOut, gwIhShift
          character(len=256) :: route_topo_f=""
          character(len=256) :: route_chan_f=""
          character(len=256) :: route_link_f=""
          character(len=256) :: route_lake_f=""
          character(len=256) :: route_direction_f=""
          character(len=256) :: route_order_f=""
          character(len=256) :: gwbasmskfil =""
          character(len=256) :: gwstrmfil =""
          character(len=256) :: geo_finegrid_flnm =""
          character(len=256) :: rt2d_flnm =""  ! Qian, 2016.06
          character(len=256) :: RESTART_OUTNAME = "" ! Qian, 2016.11
          character(len=256) :: GW2D_OUTNAME = ""    ! Qian, 2016.11
          character(len=256) :: RT2D_OUTNAME = ""    ! Qian, 2016.11 

          integer ::frxst_pts_out            ! ASCII point timeseries output at user specified points
          integer ::CHRTOUT_DOMAIN           ! Netcdf point timeseries output at all channel points
          integer ::CHRTOUT_GRID             ! Netcdf grid of channel streamflow values
          integer ::LSMOUT_DOMAN             ! Netcdf grid of variables passed between LSM and routing components
          integer ::RTOUT_DOMAIN             ! Netcdf grid of terrain routing variables on routing grid
          integer ::output_gw                ! Netcdf grid of GW
          integer ::outlake                  ! Netcdf grid of lake
	  integer ::output_rt2d		     ! Netcdf grid of rt2d

   END TYPE namelist_rt_field 

