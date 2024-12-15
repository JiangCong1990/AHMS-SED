#!/bin/csh -f

if(! $?WRF_HYDRO) setenv WRF_HYDRO 1

# update code
rm -rf  LandModel LandModel_cpl
if(! -e Land_models/NoahMP/MPP) then
  cd Land_models/NoahMP
  ln -sf ../../MPP .
  cd ../..
endif
ln -sf Land_models/NoahMP LandModel
ln -sf CPL/NoahMP_cpl LandModel_cpl

# make clean  
make clean; rm -f Run/wrf_hydro_NoahMP.exe

cat macros LandModel/user_build_options.bak > LandModel/user_build_options

# environment variables
export HYDRO_D=0
export WRF_HYDRO=1

# make to compile 
make

# update executable file, *.exe 
cd Run
mv  wrf_hydro.exe wrf_hydro_NoahMP.exe; ln -sf wrf_hydro_NoahMP.exe wrf_hydro.exe
