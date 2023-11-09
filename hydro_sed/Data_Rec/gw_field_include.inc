
 type gw_field
      integer :: ix        ! grid number in west-east
      integer :: jx        ! grid number in south-north
      integer :: allo_status = -99

      real :: dx           ! resolution, m
      real :: dt           ! timestep interval, sec.
      real :: tarea_basin  ! basin area, m2

      integer, allocatable, dimension(:,:) ::  &
        ltype              ! land-sfc type

      real,    allocatable, dimension(:,:) ::  &
        elev,           &  ! elev/bathymetry of sfc rel to sl (m)
        bot,            &  ! lowest groundwater head (m) 
        aqthick,        &  ! aquifer thickness (m)
        hycond,         &  ! hydraulic conductivity (m/s per m/m)
        poros,          &  ! porosity (m3/m3)
        compres,        &  ! compressibility (1/Pa)
        ho                 ! head at start of timestep (m)

      real,    allocatable, dimension(:,:) ::  &
        h,              &  ! head, after ghmcompute (m)
        convgw,         &  ! convergence due to gw flow (m/s)
        excess             ! surface exceeding groundwater (mm)

      real,    allocatable, dimension(:,:) ::  &
	qdarcyRT,       &  ! approximated flux between soil and groundwater for coupled simulations on routing grid
	qsgwrt,         &  ! flux between soil and groundwater for coupled simulations on routing grid
	qsgw,           &  ! flux between soil and groundwater for coupled simulations on lsm grid
	qgw_chanrt         ! flux between groundwater and channel

      real    :: ebot      ! total fix for groundwater head, m/s
      real    :: eocn      ! total flux into ocean from groundwater, m/s
      integer :: istep = 0 ! gwstep timestep counter
      
      real    :: its, ite, jts, jte

  end type gw_field

