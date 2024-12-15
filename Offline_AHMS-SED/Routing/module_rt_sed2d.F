!
!                   *** Purpose ***
!This subroutine computes splash erosion by raindrop impact and flow erosion by overland flow
!-----------------------------------------------------------------------------------------------
!   Splash erosion model is adopted from EUROSEM model developed by Morgan (2001).
!   Rill/interill erosion model is adoped from Modified ANSWERS model by Park et al.(1982)
!   C.Jiang,2019.08, Modified from J. Jeong and N. Kannan, BRC.(SWAT)
!------------------------------------------------------------------------------------------------

subroutine rt_sed
IMPLICIT NONE
!----------------------------INPUT VARIABLES---------------------------------
INTEGER, INTENT(IN) :: did
REAL   , INTENT(IN) :: dx
REAL   ,   parameter  :: botinc   = 1.0
REAL,   DIMENSION(:,:), ALLOCATABLE :: roughness
    cht(:)      !m             !canopy height
    fimp(:)     !fraction      !fraction of HRU area that is
                               !impervious (both directly and
                               !indirectly connected)
    hhqday(:)   !mm H2O        !surface runoff generated each timestep 
                               !of day in HRU
    hru_km(:)   !km2           !area of HRU in square kilometers
    idt         !minutes       !length of time step used to report
 
REAL   , INTENT(IN) :: dt    ! time interval for HYDRO, = nlst_rt(did)%dtrt



    inum1       !none          !subbasin number
    laiday(:)   !m2/m2         !leaf area index
    rainsub(:,:)!mm H2O        !precipitation for the time step during the
                               !day in HRU
    PREC
    eros_spl    !none        !coefficient of splash erosion varing 0.9-3.1
    urblu(:)    !none          !urban land type identification number from
                               !urban.dat
    usle_k(:)   !              !USLE equation soil erodibility (K) factor

!----------------------------OUTPUT VARIABLES---------------------------------
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hhsedy(:,:)|tons           |sediment yield from HRU drung a time step
!!                               |applied to HRU


!----------------------------LOCAL VARIABLES---------------------------------
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
REAL                                    :: FFF    !runoff decay factor (m-1)
!                               UNITS               Definition
REAL          :: bed_shear      !N/m2		   |shear stress b/w stream bed and flow	
REAL          :: erod_k         !g/J		   |soil detachability value	
INTEGER       :: jj             !none              |HRU number
INTEGER       ::  kk              !none              |time step of the day
!!	  ke_direct	::	!J/m2/mm	   |rainfall kinetic energy of direct throughfall
!!	  ke_leaf	::	!J/m2/mm	   |rainfall kinetic energy of leaf drainage
!!	  ke_total	::	!J/m2   	   |total kinetic energy of rainfall
!!	  percent_clay::	!percent	   |percent clay
!!	  percent_sand::	!percent	   |percent sand
!!	  percent_silt::	!percent	   |percent silt
!!	  pheff     ::	!m			   |effective plant height
!!	  rdepth_direct::	!mm			   |rainfall depth of direct throughfall
!!	  rdepth_leaf::	!mm			   |rainfall depth of leaf drainage
!!	  rdepth_tot::	!mm			   |total rainfall depth 
!!    rintnsty	  ::  !mm/hr         |rainfall intensity
!!	  sedspl::		|tons		   |sediment yield by rainfall impact during time step
!!	  sedov ::		|tons		   |sediment yield by overland flow during time step
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


      use parm

	integer :: k, j
	real :: erod_k
	real :: ke_direct, ke_leaf, ke_total,pheff, c
	real :: rdepth_direct, rdepth_leaf, rdepth_tot, canopy_cover
	real :: bed_shear, sedov, sedspl, rain_d50, rintnsty

	j = ihru
	
!! Fraction of sand
      percent_clay = sol_clay(1,j)
	percent_silt = sol_silt(1,j) 
	percent_sand = 100. - percent_clay - percent_silt

!! Soil detachability values adopted from EUROSEM User Guide (Table 1)
	if (ISLTYP=1)  then
	  erod_k = 1.9 !Sand
      elseif (ISLTYP=2) then
	  erod_k = 3.0 !Loamy Sand 
      elseif (ISLTYP=3) then
	  erod_k = 2.6 !Sandy Loam
      elseif (ISLTYP=4) then
	  erod_k = 1.5 !Silt Loam
      elseif (ISLTYP=5) then
	  erod_k = 1.2 !Silt
      elseif (ISLTYP=6) then
	  erod_k = 2.0 !Loam
      elseif (ISLTYP=7) then
	  erod_k = 2.1 !Sandy Clay Loamy
      elseif (ISLTYP=8) then
	  erod_k = 1.6 !Silty Clay Loam
      elseif (ISLTYP=9) then
	  erod_k = 1.7 !Clay Loam
      elseif (ISLTYP=10) then
	  erod_k = 1.9 !Sandy Clay
      elseif (ISLTYP=11) then
	  erod_k = 1.6 !Silty Clay
      elseif (ISLTYP=12) then
	  erod_k = 2.0 !Clay
      else
     erod_k = 0.!not soil 
      end if
!------------------------------------
!!canopy cover based on leaf area index
!!canopy cover is assumed to be 100% if LAT>=1

	  if(LAI(j)>=1.) then
	    canopy_cover = 1.
	  else
	    canopy_cover = LAI(j)
	  end if

!! Rainfall kinetic energy generated by direct throughfall (J/m^2/mm)
        
        ke_direct = 8.95 + 8.44 * log10(rintnsty)

        if(ke_direct<0.) ke_direct = 0.

!! Rainfall kinetic energy generated by leaf drainage (J/m^2)

	    pheff = 0.5 * cht(j)
	    ke_leaf = 15.8 * pheff ** 0.5 - 5.87
	    if (ke_leaf<0) ke_leaf = 0.
	  !! Depth of rainfall
	    rdepth_tot = rainsub(j,k) / (idt * 60.)
	    rdepth_leaf = rdepth_tot * canopy_cover 
	    rdepth_direct = rdepth_tot - rdepth_leaf 
	  else
	    ke_direct = 0.
	    ke_leaf = 0.
		rdepth_tot = 0.
	    rdepth_leaf = 0.
	    rdepth_direct = 0.
	  endif







t
