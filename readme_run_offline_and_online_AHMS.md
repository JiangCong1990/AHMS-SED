! Qian Xia, 2018.01
! Cong Jiang, 2023.01
! ********************************************************************************************
! To configure, compile and run offline AHMS (NoahMP-HMS) or online AHMS
! ********************************************************************************************
0-1. configure and compile WRF, and get real.exe & wrf.exe in WRFV3/main
	steps:	do steps by steps as WRF User Guide shows.
	Note:	to coupled WRF and HYDRO, the environmental variables should be set as 
		WRF_EM_CORE = 1, WRF_HYDRO = 1, HYDRO_D = 1 

0-2. configure and compile WPS, and get geogrid.exe, ungrib.exe, & metgrid.exe in WPS
	steps:	do steps by steps as WRF User Guide shows.

0-3. configure and compile NoahMP-HMS, and get wrf_hydro_NoahMP.exe  
	steps:  compile NoahMP-HMS in WRFV3/hydro by ./compile_offline_NoahMP.csh
	Note:	to couple Noah-MP and HYDRO, the environmental variables must be WRF_HYDRO = 1, HYDRO_D = 1 

! ********************************************************************************************
! To run offline AHMS (NoahMP-HMS) or online AHMS, the following datasets are required.
! For NoahMP-HMS: geo_em.d01.nc, wrfinput_d01, hydro_basin.nc, YYYYMMDDHH.LDASIN_DOMAIN1
! For AHMS:       geo_em.d01.nc, wrfinput_d01, hydro_basin.nc, wrfbdy_d01 
! ********************************************************************************************
1-1. geo_em.d01.nc, providing geographyical data
	steps:  setup in WPS/namelist.wps (&geogrid)
		./geogrid.exe  

1-2. wrfinput_d01, providing land surface data 
	steps:  setup in WPS/namelist.wps (&ungrib)
		./ungrib.exe

		setup in WPS/namelist.wps (metgrid)
		./metgrid.exe

		setup in WRFV3/wrf_run/namelist.input
		ln -sf ../main/real.exe .
		ln -sf ../main/wrf.exe .
		./real.exe

	Note: 	As only wrfinput_d01 is required for running NoahMP-HMS, the simulation period does not matter. 
		So, you can specify the running time as short as possible.
		For AHMS, however, you need to specify the starting and ending time exactly as you want to simulate. 

1-3. wrfbdy_d01, providing meteorological data for AHMS (only used when running AHMS)
	steps:	obtain wrfinput_d01 and wrfbdy_d01 after successfully executing real.exe

1-4. YYYYMMDDHH.LDASIN_DOMAIN1, providing hourly meteorological focring data for offline NoahMP-HMS (only used when running NoahMP-HMS)
	a. steps: setting in NCEP_6HOUR_1HOUR_obs/include/*.h
		  compile and run global_NCEP_obs_domain_ATM.f90

        b. steps: ncl 'interp_opt="bilinear"' 'srcGridName="input_files/*****"' 'dstGridName="geo_em.d01.nc"' ERA52WRFHydro_generate_weights.ncl
                  ncl 'srcFileName="*****"' 'dstGridName="geo_em.d01.nc"' ERA52WRFHydro_regrid.ncl	            

1-5. hydro_basin.nc, providing distributed hydrological parameters
	steps:	specify input filenames in create_netCDF_for_hydro_input/create_hydro_netcdf.f90
		compile and run create_hydro_netcdf.f90 


! ********************************************************************************************
! To run offline AHMS (NoahMP-HMS) 
! ******************************************************************************************** 
2-1. WRFV3/hydro/Run/hydro_parameter.input, specifying values of hydrological parameters

2-2. WRFV3/hydro/Run/hydro_transects.input, specifying longitude & latitude for cross-sections  

2-3. WRFV3/hydro/Run/namelist.hrldas,       simulation setup for land surface model Noah-MP

2-4. WRFV3/hydro/Run/hydro.namelist,        simulation setup for hydrological model

2-5. run NoahMP-HMS in WRFV3/hydro/Run by ./wrf_hydro_NoahMP.exe


! ********************************************************************************************
! To run online AHMS
! ********************************************************************************************
3-1. WRFV3/wrf_run/hydro_parameter.input, specifying values of hydrological parameters

3-2. WRFV3/wrf_run/hydro_transects.input, specifying longitude & latitude for cross-sections

3-3. WRFV3/wrf_run/namelist.input,        simulation setup for WRF

3-4. WRFV3/wrf_run/hydro.namelist,        simulation setup for hydrological model

3-5. run AHMS in WRFV3/wrf_run by ./wrf.exe

