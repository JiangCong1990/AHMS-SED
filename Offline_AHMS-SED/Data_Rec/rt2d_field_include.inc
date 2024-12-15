
 TYPE rt2d_field

! Qian, add for 2D Routing

        integer :: allo_status = -99
        real    :: dx             ! resoultion, m, dx = dy
        real    :: dt             ! timestep interval for routing, from hydro.namelist 
        integer :: ix             ! grid number in west-east
        integer :: jx             ! grid number in south-north
        integer :: method         ! numerical mothod for Overland-Channel routing
	integer :: count_flowout=0! counter for averaging-streamflow output 
        real    :: its, ite, jts, jte
        integer :: istep = 0
        integer :: NumFracts      ! Maximum number of size fractions found any of the
                                  ! soil types defined in the soil type input grid


        integer, allocatable, dimension(:,:) :: &
                ltype,          & ! land surface type
                ielev,          & ! sub-basin index
                istreamflux       ! quadrant for streamflux (1-8)

        real,    allocatable, dimension(:,:) :: &
                elev,           & ! elevation of hydrologic topography [m] 
                topo,           & ! elevation of topography [m]
		hlake, 		& ! elevation of river/lake surface [m]
		bank,		& ! depth of streambed below local elevation [m]
                sdep,           & ! Storage depth of water [m] !C.Jiang,2019.10
	        width, 		& ! stream width [m]
                width_obs,      & ! stream width from obs [m]  !C.Jiang,2020.01
		fbed, 		& ! fractional area of streambed in grid [m2/m2]
                rough,          & ! Manning roughness coefficient
                hycond,         & ! saturated hydraulic conductivity [m/s]
                hycond_g,        & !Cong,saturated hydraulic conductivity cgi [m/s]
                hycond_u,        & !Cong,saturated hydraulic conductivity cui [m/s]
                poros,          & ! porosity [m3/m3]
		aqthick,        & ! unconfined aquifer thickness [m]
                bot,            & ! lowest groundwater head [m]
		slope,          & ! slope 
                GINDEX,         & ! Gully index
                IRFRACT,        & !
		sdelev,         & ! standard deviation of elevation
                sbedthick,      & ! thickness of streambed [m]
                sabs,           & ! slope abs
                LAI,            & ! leaf area index
                kfactor           ! k factor soil erodibility


        real(kind=8), allocatable, dimension(:,:) :: &
		volaq             ! water storage of aquifer (below bottom of soil layers) [m]

        real,    allocatable, dimension(:,:) :: &
                irg,            & !
                actirg,         & ! actual irrigation with constrained by available water in river
                irgswfract,     & !
                irgrw,          & !
                irggw             !

        real,    allocatable, dimension(:,:) :: &
                fwater,         & ! fractional area of water in grid [m2/m2]
                fsat,           & ! fractional saturated area in grid [m2/m2]
		dvolgw,         & ! change of groundwater amount [m]
                hbeg,           & ! groundwater head at the beginning of timestep [m]
                hend,           & ! groundwater head at the beginning of timestep [m]
                zhead,          & ! depth of groundwater head below ground surface [m]
                precip,         & ! precipitation rate [m/s]
                evap,           & ! evaporation rate [m/s]
                runoff,         & ! surface runoff rate [m/s]
                                  ! when coupling with complete HYDRO, this runoff means infiltration-exess runoff rate 
                                  ! when only coupling with surface routing, this runoff means surface runoff rate 
                sinfil,         & ! soil drainage rate [m/s] (+ to aquifer)
                tinfil,         & ! total soil infiltration rate [m/s] (+ to soil)
                rsat,           & ! saturation excess runoff rate [m/s]
                cgiflux,        & ! cell-mean flux between stream and aquifer (+ to aquifer) [m/s]
                cuiflux,        & ! cell-mean flux,stream -> unsaturated -> aquiferi (+ to aquifer) [m/s]
                cglflux,        & ! cell-mean flux between lake and aquifer (+ to aquifer) [m/s]
                disch,          & ! running mean streamflux throungh cell, [m3/s]              
                hlakeav,        & ! annual average hlake for gw spinup mode
                sinfilav,       & ! annual average sinfil for gw spinup mode  	
	        streamflux,  	& ! stream flux in grid [m3/s]
                veloc,          & !
                qsKR,           & ! Transport Capacity using the Kilinc-Richardson equation
                t_SedFlux,      & ! 
                T_SedConc

         real,    allocatable, dimension(:,:,:) :: &
                qsEH,           &
                SusSed,         &
                DepSed,         &
                ErodSed,        &
                SedFlux,        &
                SedConc


        real,    allocatable, dimension(:,:) :: &
                w_sf,          &  ! grid-mean surface water storage [m]
                w_us,          &  ! grid-mean water storage in unsaturated zone [m]
                w_gw,          &  ! grid-mean water storage in saturated zone [m]
                ac_precip,     &  ! accumulated precipitation [m]
                ac_evap,       &  ! accumulated evaporation [m]
                ac_irgrw,      &  ! accumulated irrigation from river water [m]
                ac_irggw,      &  ! accumulated irrigation from groundwater [m]
                ac_pmewsf,     &  ! accumulated P-E over water surface [m]
                ac_infil_grnd, &  ! accumulated infiltration at soil surface [m]
		ac_drain_grnd, &  ! accumulated soil drainage for soil column [m]  
                ac_roff_grnd,  &  ! accumulated surface runoff at soil surface [m]
                ac_rsat_grnd,  &  ! accumulated saturation-excess runoff at soil surface [m]
                ac_cgi,        &  ! accumulated cgiflux*dt [m]
                ac_cgl,        &  ! accumulated cglflux*dt [m]
                ac_cui,        &  ! accumulated cuiglux*dt [m]
                ac_precevap,   &  ! accumulated precipitation minus evaporation on the water surface river/lake [m] 
                ac_fixsf,      &  ! accumulated water to fix surface water level [m3]
                ac_fixgw,      &  ! accumulated water to fix groundwater [m3]
                ac_ocnsf,      &  ! accumulated water to ocean through rivers [m3]
                ac_ocngw          ! accumulated water to ocean through groundwater [m3]

        real,    allocatable, dimension(:,:,:)  :: &
                atran             ! coordinates (lon,lat) of two ending points in transects [deg]

        real,    allocatable, dimension(:,:)    :: &
                xtran,         &  ! x-grid of two points in transects
                ytran             ! y-grid of two points in transects

        integer, allocatable, dimension(:,:)    :: &
                itran,         &
                jtran,         &
                ktran      
     
        real,    allocatable, dimension(:)      :: &
                fluxtran,      &  ! flux of specified transects 
                sedfluxtranav, &  ! averaging flux of specified transects
                T_SedConctranav, &
		fluxtranav        ! averaging flux of specified transects

        integer, allocatable, dimension(:)      :: &
                ncross            ! number of vectors crossing correspoing transect 


 END TYPE rt2d_field
