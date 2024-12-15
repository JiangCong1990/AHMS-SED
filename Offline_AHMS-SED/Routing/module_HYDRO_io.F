module module_HYDRO_io
#ifdef MPP_LAND
   use module_mpp_land
#endif
   use module_HYDRO_utils, only: get_dist_ll
   use module_namelist, only: nlst_rt
   use module_RT_data, only: rt_domain
   use module_gw_gw2d_data, only: gw2d   
   USE module_rt_rt2d_data, ONLY: rt2d    !Qian, 2016.06  
 
   implicit none
#include <netcdf.inc>

     contains
        integer function get2d_real(var_name,out_buff,ix,jx,fileName)
          implicit none
          integer :: ivar, iret,varid,ncid,ix,jx
          real out_buff(ix,jx)
          character(len=*), intent(in) :: var_name
          character(len=*), intent(in) :: fileName
          get2d_real = -1

#ifdef PARALLELIO
          iret = nf90_open_par(trim(fileName), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#else
          iret = nf_open(trim(fileName), NF_NOWRITE, ncid)
#endif
          if (iret .ne. 0) then
#ifdef HYDRO_D 
            print*,"failed to open the netcdf file: ",trim(fileName)
#endif
            out_buff = -9999.
            return 
          endif
          ivar = nf_inq_varid(ncid,trim(var_name),  varid)
          if(ivar .ne. 0) then
            ivar = nf_inq_varid(ncid,trim(var_name//"_M"),  varid)
            if(ivar .ne. 0) then
#ifdef HYDRO_D
               write(6,*) "Read Variable Error file: ",trim(fileName)
               write(6,*) "Read Error: could not find ",trim(var_name)
#endif
                 return 
            endif
          end if
          iret = nf_get_var_real(ncid, varid, out_buff)
          iret = nf_close(ncid)
          get2d_real =  ivar
      end function get2d_real
     
     subroutine get2d_lsm_real(var_name,out_buff,ix,jx,fileName)
         implicit none
         integer ix,jx, status
         character (len=*),intent(in) :: var_name, fileName
         real,dimension(ix,jx):: out_buff

#ifdef MPP_LAND
         real,allocatable, dimension(:,:) :: buff_g
#ifndef PARALLELIO
         write(6,*) "start to read variable ", var_name
         if(my_id .eq. IO_id) then
            allocate(buff_g(global_nx,global_ny) )
            status = get2d_real(var_name,buff_g,global_nx,global_ny,fileName)
         else
            allocate(buff_g(1,1) )
         end if

         call decompose_data_real(buff_g,out_buff)
         deallocate(buff_g)
#else
         status = get2d_real(var_name,out_buff,ix,jx,fileName)
#endif
#else         
         status = get2d_real(var_name,out_buff,ix,jx,fileName)
#endif
     end subroutine get2d_lsm_real

     subroutine get2d_lsm_vegtyp(out_buff,ix,jx,fileName)
         implicit none
         integer ix,jx, status,land_cat, iret, dimid,ncid
         character (len=*),intent(in) :: fileName
         character (len=256) units 
         integer,dimension(ix,jx):: out_buff
         real, dimension(ix,jx) :: xdum
#ifdef MPP_LAND
         real,allocatable, dimension(:,:) :: buff_g


#ifndef PARALLELIO
         if(my_id .eq. IO_id) then
            allocate(buff_g (global_nx,global_ny) )
         else
            allocate(buff_g (1,1) )
         endif
         if(my_id .eq. IO_id) then
#endif
#endif
                ! Open the NetCDF file.
#ifdef PARALLELIO
              iret = nf90_open_par(trim(fileName), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#else
              iret = nf_open(fileName, NF_NOWRITE, ncid)
#endif
              if (iret /= 0) then
                 write(*,'("Problem opening geo_static file: ''", A, "''")') &
                      trim(fileName)
                 call hydro_stop("get2d_lsm_vegtyp")
              endif

            iret = nf_inq_dimid(ncid, "land_cat", dimid)
            if (iret /= 0) then
              iret = nf_inq_dimid(ncid, "land_cat_stag", dimid)  ! Qian
              if (iret /= 0) then
                 print*, "nf_inq_dimid:  land_cat"
                 call hydro_stop("get2d_lsm_vegtyp")
              endif
             endif

            iret = nf_inq_dimlen(ncid, dimid, land_cat)
            if (iret /= 0) then
               print*, "nf_inq_dimlen:  land_cat"
               call hydro_stop("get2d_lsm_vegtyp")
            endif

#ifdef MPP_LAND
#ifndef PARALLELIO
            call get_landuse_netcdf(ncid, buff_g, units, global_nx ,global_ny, land_cat)
         end if
         call decompose_data_real(buff_g,xdum)     
         deallocate(buff_g)
#else
          call get_landuse_netcdf(ncid, xdum,   units, ix, jx, land_cat)
#endif
          iret = nf_close(ncid)
          
#else         
          call get_landuse_netcdf(ncid, xdum,   units, ix, jx, land_cat)
          iret = nf_close(ncid)
#endif
         out_buff = nint(xdum)
     end subroutine get2d_lsm_vegtyp



     subroutine get_file_dimension(fileName, ix,jx)
            implicit none
            character(len=*) fileName
            integer ncid , iret, ix,jx, dimid
#ifdef MPP_LAND
#ifndef PARALLELIO
            if(my_id .eq. IO_id) then
#endif
#endif
#ifdef PARALLELIO
            iret = nf90_open_par(trim(fileName), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#else
            iret = nf_open(fileName, NF_NOWRITE, ncid)
#endif
            if (iret /= 0) then
               write(*,'("Problem opening geo_static file: ''", A, "''")') &
                    trim(fileName)
               call hydro_stop("get_file_dimension")
            endif
        
            iret = nf_inq_dimid(ncid, "west_east", dimid)
        
            if (iret /= 0) then
               print*, "nf_inq_dimid:  west_east"
               call hydro_stop("get_file_dimension")
            endif
        
            iret = nf_inq_dimlen(ncid, dimid, ix)
            if (iret /= 0) then
               print*, "nf_inq_dimlen:  west_east"
               call hydro_stop("get_file_dimension")
            endif
        
            iret = nf_inq_dimid(ncid, "south_north", dimid)
            if (iret /= 0) then
                       print*, "nf_inq_dimid:  south_north"
                       call hydro_stop("get_file_dimension")
            endif
        
            iret = nf_inq_dimlen(ncid, dimid, jx)
            if (iret /= 0) then
               print*, "nf_inq_dimlen:  south_north"
               call hydro_stop("get_file_dimension")
            endif
            iret = nf_close(ncid)
#ifdef MPP_LAND
#ifndef PARALLELIO
            endif
            call mpp_land_bcast_int1(ix)
            call mpp_land_bcast_int1(jx)
#endif
#endif

     end subroutine get_file_dimension

     subroutine get2d_lsm_soltyp(out_buff,ix,jx,fileName)
         implicit none
         integer ix,jx, status,land_cat, iret, dimid,ncid
         character (len=*),intent(in) :: fileName
         character (len=256) units 
         integer,dimension(ix,jx):: out_buff
         real, dimension(ix,jx) :: xdum
#ifdef MPP_LAND
#ifndef PARALLELIO
         real,allocatable, dimension(:,:) :: buff_g

         allocate(buff_g (global_nx,global_ny) )

         if(my_id .eq. IO_id) then
#endif
#endif
                ! Open the NetCDF file.
#ifdef PARALLELIO
            iret = nf90_open_par(trim(fileName), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#else
            iret = nf_open(fileName, NF_NOWRITE, ncid)
#endif
              if (iret /= 0) then
                 write(*,'("Problem opening geo_static file: ''", A, "''")') &
                      trim(fileName)
                 call hydro_stop("get2d_lsm_soltyp")
              endif

            iret = nf_inq_dimid(ncid, "soil_cat", dimid)
            if (iret /= 0) then
                iret = nf_inq_dimid(ncid, "soil_cat_stag", dimid)  ! Qian
                if (iret /= 0) then
                   print*, "nf_inq_dimid:  soil_cat"
                   call hydro_stop("get2d_lsm_soltyp")
                endif
            endif

            iret = nf_inq_dimlen(ncid, dimid, land_cat)
            if (iret /= 0) then
               print*, "nf_inq_dimlen:  soil_cat"
               call hydro_stop("get2d_lsm_soltyp")
            endif

#ifdef MPP_LAND
#ifndef PARALLELIO
            call get_soilcat_netcdf(ncid, buff_g, units, global_nx ,global_ny, land_cat)
         end if
         call decompose_data_real(buff_g,xdum)     
         deallocate(buff_g)
#else
          call get_soilcat_netcdf(ncid, xdum,   units, ix, jx, land_cat)
#endif
          iret = nf_close(ncid)
#else         
          call get_soilcat_netcdf(ncid, xdum,   units, ix, jx, land_cat)
          iret = nf_close(ncid)
#endif
          out_buff = nint(xdum)
     end subroutine get2d_lsm_soltyp


  subroutine get_landuse_netcdf(ncid, array, units, idim, jdim, ldim)
    implicit none
#include <netcdf.inc>
    integer, intent(in) :: ncid
    integer, intent(in) :: idim, jdim, ldim
    real, dimension(idim,jdim), intent(out) :: array
    character(len=256), intent(out) :: units
    integer :: iret, varid
    real, dimension(idim,jdim,ldim) :: xtmp
    integer, dimension(1) :: mp
    integer :: i, j, l
    character(len=24), parameter :: name = "LANDUSEF"

    units = ""

    iret = nf_inq_varid(ncid,  trim(name),  varid)
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       print*, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_landuse_netcdf:  nf_inq_varid"
       call hydro_stop("get_landuse_netcdf")
    endif

    iret = nf_get_var_real(ncid, varid, xtmp)
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       print*, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_landuse_netcdf:  nf_get_var_real"
       call hydro_stop("get_landuse_netcdf")
    endif

    do i = 1, idim
       do j = 1, jdim
          mp = maxloc(xtmp(i,j,:))
          array(i,j) = mp(1)
          do l = 1,ldim
            if(xtmp(i,j,l).lt.0) array(i,j) = -9999.0
          enddo
       enddo
    enddo

  end subroutine get_landuse_netcdf


  subroutine get_soilcat_netcdf(ncid, array, units, idim, jdim, ldim)
    implicit none
#include <netcdf.inc>

    integer, intent(in) :: ncid
    integer, intent(in) :: idim, jdim, ldim
    real, dimension(idim,jdim), intent(out) :: array
    character(len=256), intent(out) :: units
    integer :: iret, varid
    real, dimension(idim,jdim,ldim) :: xtmp
    integer, dimension(1) :: mp
    integer :: i, j, did
    character(len=24), parameter :: name = "SOILCTOP"

    did = 1
    units = ""

    iret = nf_inq_varid(ncid,  trim(name),  varid)
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       print*, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_soilcat_netcdf:  nf_inq_varid"
       call hydro_stop("get_soilcat_netcdf")
    endif

    iret = nf_get_var_real(ncid, varid, xtmp)
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       print*, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_soilcat_netcdf:  nf_get_var_real"
       call hydro_stop("get_soilcat_netcdf")
    endif

    do i = 1, idim
       do j = 1, jdim
          mp = maxloc(xtmp(i,j,:))
          array(i,j) = mp(1)
       enddo
    enddo

!     if(nlst_rt(did)%GWBASESWCRT .ne. 3) then
!        where (array == 14) array = 1   ! DJG remove all 'water' soils...
!     endif

  end subroutine get_soilcat_netcdf


subroutine get_greenfrac_netcdf(ncid, array3, units, idim, jdim, ldim,mm,dd)
    implicit none
#include <netcdf.inc>
    integer, intent(in) :: ncid,mm,dd
    integer, intent(in) :: idim, jdim, ldim
    real, dimension(idim,jdim) :: array
    real, dimension(idim,jdim) :: array2
    real, dimension(idim,jdim) :: diff
    real, dimension(idim,jdim), intent(out) :: array3
    character(len=256), intent(out) :: units
    integer :: iret, varid
    real, dimension(idim,jdim,ldim) :: xtmp
    integer, dimension(1) :: mp
    integer :: i, j, mm2,daytot
    real :: ddfrac
    character(len=24), parameter :: name = "GREENFRAC"

    units = "fraction"

    iret = nf_inq_varid(ncid,  trim(name),  varid)
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       print*, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_greenfrac_netcdf:  nf_inq_varid"
       call hydro_stop("get_greenfrac_netcdf")
    endif

    iret = nf_get_var_real(ncid, varid, xtmp)
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       print*, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_greenfrac_netcdf:  nf_get_var_real"
       call hydro_stop("get_greenfrac_netcdf")
    endif


    if (mm.lt.12) then 
      mm2 = mm+1
    else
      mm2 = 1
    end if

!DJG_DES Set up dates for daily interpolation...
          if (mm.eq.1.OR.mm.eq.3.OR.mm.eq.5.OR.mm.eq.7.OR.mm.eq.8.OR.mm.eq.10.OR.mm.eq.12) then
             daytot = 31
          else if (mm.eq.4.OR.mm.eq.6.OR.mm.eq.9.OR.mm.eq.11) then 
             daytot = 30
          else if (mm.eq.2) then
             daytot = 28
          end if
          ddfrac = float(dd)/float(daytot)
          if (ddfrac.gt.1.0) ddfrac = 1.0   ! Assumes Feb. 29th change is same as Feb 28th

    print *,"DJG_DES Made it past netcdf read...month = ",mm,mm2,dd,daytot,ddfrac

    do i = 1, idim
       do j = 1, jdim
          array(i,j) = xtmp(i,j,mm)   !GREENFRAC in geogrid in units of fraction from month 1
          array2(i,j) = xtmp(i,j,mm2)   !GREENFRAC in geogrid in units of fraction from month 1
          diff(i,j) = array2(i,j) - array(i,j)
          array3(i,j) = array(i,j) + ddfrac * diff(i,j) 
       enddo
    enddo

end subroutine get_greenfrac_netcdf



subroutine get_albedo12m_netcdf(ncid, array3, units, idim, jdim, ldim,mm,dd)
    implicit none
#include <netcdf.inc>
    integer, intent(in) :: ncid,mm,dd
    integer, intent(in) :: idim, jdim, ldim
    real, dimension(idim,jdim) :: array
    real, dimension(idim,jdim) :: array2
    real, dimension(idim,jdim) :: diff
    real, dimension(idim,jdim), intent(out) :: array3
    character(len=256), intent(out) :: units
    integer :: iret, varid
    real, dimension(idim,jdim,ldim) :: xtmp
    integer, dimension(1) :: mp
    integer :: i, j, mm2,daytot
    real :: ddfrac
    character(len=24), parameter :: name = "ALBEDO12M"


    units = "fraction"

    iret = nf_inq_varid(ncid,  trim(name),  varid)
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       print*, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_albedo12m_netcdf:  nf_inq_varid"
       call hydro_stop("get_albedo12m_netcdf")
    endif

    iret = nf_get_var_real(ncid, varid, xtmp)
    if (iret /= 0) then
       print*, 'name = "', trim(name)//'"'
       print*, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_albedo12m_netcdf:  nf_get_var_real"
       call hydro_stop("get_albedo12m_netcdf")
    endif

    if (mm.lt.12) then 
      mm2 = mm+1
    else
      mm2 = 1
    end if

!DJG_DES Set up dates for daily interpolation...
          if (mm.eq.1.OR.mm.eq.3.OR.mm.eq.5.OR.mm.eq.7.OR.mm.eq.8.OR.mm.eq.10.OR.mm.eq.12) then
             daytot = 31
          else if (mm.eq.4.OR.mm.eq.6.OR.mm.eq.9.OR.mm.eq.11) then 
             daytot = 30
          else if (mm.eq.2) then
             daytot = 28
          end if
          ddfrac = float(dd)/float(daytot)
          if (ddfrac.gt.1.0) ddfrac = 1.0   ! Assumes Feb. 29th change is same as Feb 28th

#ifdef HYDRO_D
    print *,"DJG_DES Made it past netcdf read...month = ",mm,mm2,dd,daytot,ddfrac
#endif

    do i = 1, idim
       do j = 1, jdim
          array(i,j) = xtmp(i,j,mm) / 100.0   !Convert ALBEDO12M from % to fraction...month 1
          array2(i,j) = xtmp(i,j,mm2) / 100.0   !Convert ALBEDO12M from % to fraction... month 2
          diff(i,j) = array2(i,j) - array(i,j)
          array3(i,j) = array(i,j) + ddfrac * diff(i,j) 
       enddo
    enddo

end subroutine get_albedo12m_netcdf




  subroutine get_2d_netcdf(name, ncid, array, units, idim, jdim, &
       fatal_if_error, ierr)
    implicit none
#include <netcdf.inc>
    character(len=*), intent(in) :: name
    integer, intent(in) :: ncid
    integer, intent(in) :: idim, jdim
    real, dimension(idim,jdim), intent(out) :: array
    character(len=256), intent(out) :: units
    integer :: iret, varid
    ! .TRUE._IF_ERROR:  an input code value:
    !      .TRUE. if an error in reading the data should stop the program.
    !      Otherwise the, IERR error flag is set, but the program continues.
    logical, intent(in) :: fatal_if_error 
    integer, intent(out) :: ierr

    units = ""

    iret = nf_inq_varid(ncid,  name,  varid)

    if (iret /= 0) then
       if (fatal_IF_ERROR) then
          print*, 'name = "', trim(name)//'"'
          print*, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_2d_netcdf:  nf_inq_varid"
          call hydro_stop("get_2d_netcdf")
       else
          ierr = iret
          return
       endif
    endif


    iret = nf_get_var_real(ncid, varid, array)
    if (iret /= 0) then
       if (fatal_IF_ERROR) then
          print*, 'name = "', trim(name)//'"'
          print*, "MODULE_NOAHLSM_HRLDAS_INPUT:  get_2d_netcdf:  nf_get_var_real"
          call hydro_stop("get_2d_netcdf")
       else
          ierr = iret
          return
       endif
    endif
! ERA5 land ocean mask is inconsistent with WPS geo data
    WHERE(array > 1e19) array=0.

    ierr = 0;
  end subroutine get_2d_netcdf


      subroutine get_2d_netcdf_cows(var_name,ncid,var, &
            ix,jx,tlevel,fatal_if_error,ierr)
#include <netcdf.inc>
          character(len=*), intent(in) :: var_name
          integer,intent(in) ::  ncid,ix,jx,tlevel
          real, intent(out):: var(ix,jx)
          logical, intent(in) :: fatal_if_error
          integer ierr, iret
          integer varid
          integer start(4),count(4)
          data count /1,1,1,1/
          data start /1,1,1,1/
          count(1) = ix
          count(2) = jx
          start(4) = tlevel
      iret = nf_inq_varid(ncid,  var_name,  varid)

      if (iret /= 0) then
        if (fatal_IF_ERROR) then
           print*, "MODULE_NOAHLSM_HRLDAS_INPUT: get_2d_netcdf_cows:nf_inq_varid"
           call hydro_stop("get_2d_netcdf_cows")
        else
          ierr = iret
          return
        endif
      endif
      iret = nf_get_vara_real(ncid, varid, start,count,var)

      return
      end subroutine get_2d_netcdf_cows

!---------------------------------------------------------
!DJG Subroutinesfor inputting routing fields...
!DNY   first reads the files to get the size of the 
!DNY   LINKS arrays
!DJG   - Currently only hi-res topo is read 
!DJG   - At a future time, use this routine to input
!DJG     subgrid land-use classification or routing
!DJG     parameters 'overland roughness' and 'retention
!DJG     depth'
!
!DJG,DNY - Update this subroutine to read in channel and lake
!           parameters if activated       11.20.2005
!---------------------------------------------------------

       SUBROUTINE READ_ROUTEDIM(IXRT,JXRT,route_chan_f,route_link_f, &
            route_direction_f, route_lake_f, NLINKS, NLAKES, &
            CH_NETLNK, channel_option, geo_finegrid_flnm)

         implicit none
#include <netcdf.inc>
        INTEGER                                      :: I,J,channel_option,iret,jj
        INTEGER, INTENT(INOUT)                       :: NLINKS, NLAKES
        INTEGER, INTENT(IN)                          :: IXRT,JXRT
        INTEGER                                      :: CHNID,cnt
        INTEGER, DIMENSION(IXRT,JXRT)                :: CH_NETRT   !- binary channel mask
        INTEGER, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: CH_NETLNK  !- each node gets unique id
        INTEGER, DIMENSION(IXRT,JXRT)                :: DIRECTION  !- flow direction
        INTEGER, DIMENSION(IXRT,JXRT)                :: LAKE_MSKRT
        REAL, DIMENSION(IXRT,JXRT)                   :: LAT, LON

!!Dummy read in grids for inverted y-axis


        CHARACTER(len=256)       :: route_chan_f, route_link_f,route_direction_f,route_lake_f
        CHARACTER(len=256)       :: InputLine
        CHARACTER(len=256)       :: geo_finegrid_flnm
        CHARACTER(len=256)       :: var_name
     
        NLINKS = 0
        NLAKES = 0
        CH_NETRT = -9999
        CH_NETLNK = -9999


        cnt = 0 
#ifdef HYDRO_D
       print *, "Channel Option in Routedim is ", channel_option
#endif

       IF(channel_option.eq.3) then  !get maxnodes and links from grid

         var_name = "CHANNELGRID"
         call readRT2d_int(var_name,CH_NETRT,ixrt,jxrt,&
                   trim(geo_finegrid_flnm))
         
         var_name = "FLOWDIRECTION"
         call readRT2d_int(var_name,DIRECTION,ixrt,jxrt,&
                   trim(geo_finegrid_flnm))

         var_name = "LAKEGRID"
         call readRT2d_int(var_name,LAKE_MSKRT,ixrt,jxrt,&
                   trim(geo_finegrid_flnm))


        var_name = "LATITUDE"
        call readRT2d_real(var_name,LAT,ixrt,jxrt,&
                     trim(geo_finegrid_flnm))
        var_name = "LONGITUDE"
        call readRT2d_real(var_name,LON,ixrt,jxrt,&
                     trim(geo_finegrid_flnm))
          
! temp fix for buggy Arc export...
        do j=1,jxrt
          do i=1,ixrt
            if(DIRECTION(i,j).eq.-128) DIRECTION(i,j)=128
          end do
        end do

!DJG inv         do j=jxrt,1,-1
         do j=1,jxrt
             do i = 1, ixrt
               if (CH_NETRT(i,j) .ge.0.AND.CH_NETRT(i,j).lt.100) then 
                 NLINKS = NLINKS + 1
               endif
            end do 
         end do 
#ifdef HYDRO_D
         print *, "NLINKS IS ", NLINKS 
#endif


!DJG inv         DO j = JXRT,1,-1  !rows
         DO j = 1,JXRT  !rows
          DO i = 1 ,IXRT   !colsumns
           If (CH_NETRT(i, j) .ge. 0) then !get its direction
            If ((DIRECTION(i, j) .EQ. 64) .AND. (j+1 .LE. JXRT).AND.(CH_NETRT(i,j+1) .ge.0) ) then !North
             cnt = cnt + 1
             CH_NETLNK(i,j) = cnt
            else if ((DIRECTION(i, j) .EQ. 128) .AND. (i + 1 .LE. IXRT) &
               .AND. (j + 1 .LE. JXRT) .AND. (CH_NETRT(i+1,j+1) .ge.0)) then !North East
             cnt = cnt + 1
             CH_NETLNK(i,j) = cnt
            else if ((DIRECTION(i, j) .EQ. 1) .AND. (i + 1 .LE. IXRT).AND.(CH_NETRT(i+1,j) .ge. 0)) then !East
             cnt = cnt + 1
             CH_NETLNK(i,j) = cnt 
            else if ((DIRECTION(i, j) .EQ. 2) .AND. (i + 1 .LE. IXRT) &
                    .AND. (j - 1 .NE. 0).AND.(CH_NETRT(i+1,j-1).ge.0)) then !south east
             cnt = cnt + 1
             CH_NETLNK(i,j) = cnt
            else if ((DIRECTION(i, j) .EQ. 4).AND.(j - 1 .NE. 0).AND.(CH_NETRT(i,j-1).ge.0)) then !due south
             cnt = cnt + 1
             CH_NETLNK(i,j) = cnt
            else if ((DIRECTION(i, j) .EQ. 8) .AND. (i - 1 .GT. 0) &
                    .AND. (j - 1 .NE. 0).AND. (CH_NETRT(i-1,j-1).ge.0) ) then !south west
             cnt = cnt + 1
             CH_NETLNK(i,j) = cnt
            else if ((DIRECTION(i, j) .EQ. 16) .AND. (i - 1 .GT. 0).AND.(CH_NETRT(i-1,j).ge.0)) then !West
             cnt = cnt + 1
             CH_NETLNK(i,j) = cnt
            else if ((DIRECTION(i, j) .EQ. 32) .AND. (i - 1 .GT. 0) &
                    .AND. (j + 1 .LE. JXRT) .AND. (CH_NETRT(i-1,j+1).ge.0)) then !North West
             cnt = cnt + 1
             CH_NETLNK(i,j) = cnt 
           else 
#ifdef HYDRO_D
             write(*,135) "PrPt/LkIn", CH_NETRT(i,j), DIRECTION(i,j), LON(i,j), LAT(i,j),i,j 
135             FORMAT(A9,1X,I3,1X,I3,1X,F10.5,1X,F9.5,1X,I4,1X,I4)
#endif
             if (DIRECTION(i,j) .eq. 0) then
#ifdef HYDRO_D
               print *, "Direction i,j ",i,j," of point ", cnt, "is invalid"
#endif
             endif

           End If
         End If !CH_NETRT check for this node
        END DO
       END DO 
#ifdef HYDRO_D
       print *, "found type 0 nodes", cnt
#endif
      
!Find out if the boundaries are on an edge or flow into a lake
!DJG inv       DO j = JXRT,1,-1
       DO j = 1,JXRT
         DO i = 1 ,IXRT
          If (CH_NETRT(i, j) .ge. 0) then !get its direction

           If ( ((DIRECTION(i, j).EQ. 64) .AND. (j + 1 .GT. JXRT))         & !-- 64's can only flow north
               .OR. ((DIRECTION(i, j) .EQ. 64).and. (j<jxrt) .AND. (CH_NETRT(i,j+1) .lt. 0))) then !North
              cnt = cnt + 1
              CH_NETLNK(i,j) = cnt
#ifdef HYDRO_D
              print *, "Boundary Pour Point N", cnt,CH_NETRT(i,j), i,j
#endif
           else if ( ((DIRECTION(i, j) .EQ. 128) .AND. (i + 1 .GT. IXRT))  & !-- 128's can flow out of the North or East edge
               .OR.  ((DIRECTION(i, j) .EQ. 128) .AND. (j + 1 .GT. JXRT))  & !   this is due north edge     
               .OR.  ((DIRECTION(i, j) .EQ. 128) .AND. (i<ixrt .and. j<jxrt) .AND.(CH_NETRT(i + 1, j + 1).lt.0))) then !North East
              cnt = cnt + 1
              CH_NETLNK(i,j) = cnt
#ifdef HYDRO_D
              print *, "Boundary Pour Point NE", cnt, CH_NETRT(i,j),i,j
#endif
           else if (((DIRECTION(i, j) .EQ. 1) .AND. (i + 1 .GT. IXRT))     & !-- 1's can only flow due east
               .OR. ((DIRECTION(i, j) .EQ. 1) .and. (i<ixrt) .AND. (CH_NETRT(i + 1, j) .lt. 0))) then !East
              cnt = cnt + 1
              CH_NETLNK(i,j) = cnt
#ifdef HYDRO_D
              print *, "Boundary Pour Point E", cnt,CH_NETRT(i,j), i,j
#endif
           else if ( ((DIRECTION(i, j) .EQ. 2) .AND. (i + 1 .GT. IXRT))    &      !-- 2's can flow out of east or south edge
                .OR. ((DIRECTION(i, j) .EQ. 2) .AND. (j - 1 .EQ. 0))       &      !-- this is the south edge
                .OR. ((DIRECTION(i, j) .EQ. 2) .and. (i<ixrt .and. j>1) .AND.(CH_NETRT(i + 1, j - 1) .lt.0))) then !south east
              cnt = cnt + 1
              CH_NETLNK(i,j) = cnt
#ifdef HYDRO_D
              print *, "Boundary Pour Point SE", cnt,CH_NETRT(i,j), i,j
#endif
           else if ( ((DIRECTION(i, j) .EQ. 4) .AND. (j - 1 .EQ. 0))       &      !-- 4's can only flow due south
                .OR. ((DIRECTION(i, j) .EQ. 4) .and. (j>1) .AND.(CH_NETRT(i, j - 1) .lt. 0))) then !due south
              cnt = cnt + 1
              CH_NETLNK(i,j) = cnt
#ifdef HYDRO_D
              print *, "Boundary Pour Point S", cnt,CH_NETRT(i,j), i,j
#endif
           else if ( ((DIRECTION(i, j) .EQ. 8) .AND. (i - 1 .LE. 0))      &      !-- 8's can flow south or west
               .OR.  ((DIRECTION(i, j) .EQ. 8) .AND. (j - 1 .EQ. 0))      &      !-- this is the south edge
               .OR.  ((DIRECTION(i, j).EQ.8).and. (i>1 .and. j>1) .AND.(CH_NETRT(i - 1, j - 1).lt.0))) then !south west
              cnt = cnt + 1
              CH_NETLNK(i,j) = cnt
#ifdef HYDRO_D
              print *, "Boundary Pour Point SW", cnt,CH_NETRT(i,j), i,j
#endif
           else if ( ((DIRECTION(i, j) .EQ. 16) .AND. (i - 1 .LE. 0))       &      !-- 16's can only flow due west 
               .OR.  ((DIRECTION(i, j).EQ.16) .and. (i>1) .AND.(CH_NETRT(i - 1, j).lt.0))) then !West
              cnt = cnt + 1
              CH_NETLNK(i,j) = cnt              
#ifdef HYDRO_D
              print *, "Boundary Pour Point W", cnt,CH_NETRT(i,j), i,j
#endif
           else if ( ((DIRECTION(i, j) .EQ. 32) .AND. (i - 1 .LE. 0))      &      !-- 32's can flow either west or north
               .OR.  ((DIRECTION(i, j) .EQ. 32) .AND. (j + 1 .GT. JXRT))   &      !-- this is the north edge
               .OR.  ((DIRECTION(i, j).EQ.32) .and. (i>1 .and. j<jxrt) .AND.(CH_NETRT(i - 1, j + 1).lt.0))) then !North West
              cnt = cnt + 1
              CH_NETLNK(i,j) = cnt
#ifdef HYDRO_D
              print *, "Boundary Pour Point NW", cnt,CH_NETRT(i,j), i,j
#endif
           endif
          endif !CH_NETRT check for this node
         END DO
       END DO 

#ifdef HYDRO_D
       print *, "total number of channel elements", cnt
       print *, "total number of NLINKS          ", NLINKS
#endif



      !-- get the number of lakes
       if (cnt .ne. NLINKS) then 
#ifdef HYDRO_D
         print *, "Apparent error in network topology", cnt, NLINKS
         print* , "ixrt =", ixrt, "jxrt =", jxrt
#endif
         call hydro_stop("READ_ROUTEDIM")
       endif
!DJG inv       do j=jxrt,1,-1
       do j=1,jxrt
          do i = 1,ixrt
           if (LAKE_MSKRT(i,j) .gt. NLAKES) then 
             NLAKES = LAKE_MSKRT(i,j)
           endif
        end do
       end do
#ifdef HYDRO_D
       write(6,*) "finish read_red ..  nlakes = ", nlakes
#endif

       return

     else  ! get nlinks from the ascii file of links
        open(unit=79,file=trim(route_link_f),          & !link
                   form='formatted',status='old')
#ifdef HYDRO_D
                   write(6,*) "read file ",route_link_f
#endif
1011    read(79,*,end= 1999) InputLine
        NLINKS = NLINKS + 1
        goto 1011
1999   continue
        NLINKS = NLINKS - 1 !-- first line is a comment 
        close(79)
#ifdef HYDRO_D
        print *, "Number of Links on sparse network", NLINKS
#endif
        return
     end if

     END SUBROUTINE READ_ROUTEDIM

     SUBROUTINE nreadRT2d_real(var_name, inv, ixrt, jxrt, fileName) 
         implicit none
         INTEGER :: iret
         INTEGER, INTENT(IN) :: ixrt,jxrt
         INTEGER :: i, j, ii,jj
         CHARACTER(len=*):: var_name,fileName
         real, INTENT(OUT), dimension(ixrt,jxrt) :: inv
         real, dimension(ixrt,jxrt) :: inv_tmp
#ifdef MPP_LAND
         real, allocatable,dimension(:,:) :: g_inv_tmp, g_inv
         if(my_id .eq. io_id) then
              allocate(g_inv_tmp(global_rt_nx,global_rt_ny))
              allocate(g_inv(global_rt_nx,global_rt_ny))
              g_inv_tmp = -9999.9
              iret =  get2d_real(var_name,g_inv_tmp,global_rt_nx,global_rt_ny,&
                     trim(fileName))
              do i=1,global_rt_nx
                 jj=global_rt_ny
                 do j=1,global_rt_ny
                   g_inv(i,j)=g_inv_tmp(i,jj)
                   jj=global_rt_ny-j
                 end do
              end do
         else
              allocate(g_inv_tmp(1,1))
              allocate(g_inv(1,1))
         endif 
         call decompose_RT_real(g_inv,inv,global_rt_nx,global_rt_ny,IXRT,JXRT)
         deallocate(g_inv_tmp, g_inv)
#else
         inv_tmp = -9999.9
         iret =  get2d_real(var_name,inv_tmp,ixrt,jxrt,&
                     trim(fileName))
         do i=1,ixrt
            jj=jxrt
         do j=1,jxrt
           inv(i,j)=inv_tmp(i,jj)
           jj=jxrt-j
         end do
        end do
#endif

        
     end SUBROUTINE nreadRT2d_real

     SUBROUTINE nreadRT2d_int(var_name, inv, ixrt, jxrt, fileName) 
         implicit none
         INTEGER, INTENT(IN) :: ixrt,jxrt
         INTEGER :: i, j, ii,jj, iret
         CHARACTER(len=*):: var_name,fileName
         integer, INTENT(OUT), dimension(ixrt,jxrt) :: inv
         integer, dimension(ixrt,jxrt) :: inv_tmp
#ifdef MPP_LAND
         integer, allocatable,dimension(:,:) :: g_inv_tmp, g_inv
         if(my_id .eq. io_id) then
              allocate(g_inv_tmp(global_rt_nx,global_rt_ny))
              allocate(g_inv(global_rt_nx,global_rt_ny))
              g_inv_tmp = -9999.9
              call  get2d_int(var_name,g_inv_tmp,global_rt_nx,global_rt_ny,&
                     trim(fileName))
              do i=1,global_rt_nx
                 jj=global_rt_ny
                do j=1,global_rt_ny
                  g_inv(i,j)=g_inv_tmp(i,jj)
                  jj=global_rt_ny-j
                end do
              end do
         else
              allocate(g_inv_tmp(1,1))
              allocate(g_inv(1,1))
         endif
         call decompose_RT_int(g_inv,inv,global_rt_nx,global_rt_ny,IXRT,JXRT)
         deallocate(g_inv_tmp, g_inv)
#else
         call  get2d_int(var_name,inv_tmp,ixrt,jxrt,&
                     trim(fileName))
         do i=1,ixrt
            jj=jxrt
         do j=1,jxrt
           inv(i,j)=inv_tmp(i,jj)
           jj=jxrt-j
         end do
        end do
#endif
     end SUBROUTINE nreadRT2d_int
!---------------------------------------------------------
!DJG -----------------------------------------------------

     SUBROUTINE readRT2d_real(var_name, inv, ixrt, jxrt, fileName) 
         implicit none
         INTEGER :: iret
         INTEGER, INTENT(IN) :: ixrt,jxrt
         INTEGER :: i, j, ii,jj
         CHARACTER(len=*):: var_name,fileName
         real, INTENT(OUT), dimension(ixrt,jxrt) :: inv
         real, dimension(ixrt,jxrt) :: inv_tmp

         inv_tmp = -9999.9

         iret =  get2d_real(var_name,inv_tmp,ixrt,jxrt,&
                     trim(fileName))
         do i=1,ixrt
            jj=jxrt
         do j=1,jxrt
           inv(i,j)=inv_tmp(i,jj)
           jj=jxrt-j
         end do
        end do
     end SUBROUTINE readRT2d_real

     SUBROUTINE readRT2d_int(var_name, inv, ixrt, jxrt, fileName) 
         implicit none
         INTEGER, INTENT(IN) :: ixrt,jxrt
         INTEGER :: i, j, ii,jj
         CHARACTER(len=*):: var_name,fileName
         integer, INTENT(OUT), dimension(ixrt,jxrt) :: inv
         integer, dimension(ixrt,jxrt) :: inv_tmp
         call  get2d_int(var_name,inv_tmp,ixrt,jxrt,&
                     trim(fileName))
         do i=1,ixrt
            jj=jxrt
         do j=1,jxrt
           inv(i,j)=inv_tmp(i,jj)
           jj=jxrt-j
         end do
        end do
     end SUBROUTINE readRT2d_int
!---------------------------------------------------------
!DJG -----------------------------------------------------

#ifdef MPP_LAND
  subroutine MPP_READ_SIMP_GW(IX,JX,IXRT,JXRT,GWSUBBASMSK,gwbasmskfil,&
          gw_strm_msk,numbasns,ch_netrt,AGGFACTRT)

   USE module_mpp_land
    
    integer, intent(in)                     :: IX,JX,IXRT,JXRT,AGGFACTRT
    integer, intent(out)                    :: numbasns
    integer, intent(out), dimension(IX,JX)  :: GWSUBBASMSK
    integer, intent(out), dimension(IXRT,JXRT)  :: gw_strm_msk
    integer, intent(in), dimension(IXRT,JXRT)  :: ch_netrt
    character(len=256)                      :: gwbasmskfil
    !integer,dimension(global_nX,global_ny) ::  g_GWSUBBASMSK
    !yw integer,dimension(global_rt_nx, global_rt_ny) ::  g_gw_strm_msk,g_ch_netrt

    integer,allocatable,dimension(:,:) ::  g_GWSUBBASMSK
    integer,allocatable,dimension(:, :) ::  g_gw_strm_msk,g_ch_netrt
    
     if(my_id .eq. IO_id) then
          allocate(g_gw_strm_msk(global_rt_nx, global_rt_ny))
          allocate(g_ch_netrt(global_rt_nx, global_rt_ny))
          allocate(g_GWSUBBASMSK(global_nX,global_ny))
     else
          allocate(g_gw_strm_msk(1,1))
          allocate(g_ch_netrt(1,1))
          allocate(g_GWSUBBASMSK(1,1))
     endif


     call write_IO_rt_int(ch_netrt,g_ch_netrt)

     if(my_id .eq. IO_id) then
       call READ_SIMP_GW(global_nX,global_ny,global_rt_nx,global_rt_ny,&
             g_GWSUBBASMSK,gwbasmskfil,g_gw_strm_msk,numbasns,&
             g_ch_netrt,AGGFACTRT) 
     endif
     call decompose_data_int(g_GWSUBBASMSK,GWSUBBASMSK)
     call decompose_RT_int(g_gw_strm_msk,gw_strm_msk,  &
          global_rt_nx, global_rt_ny,ixrt,jxrt)
     call mpp_land_bcast_int1(numbasns)

     deallocate(g_gw_strm_msk)
     deallocate(g_ch_netrt)
     deallocate(g_GWSUBBASMSK)

  return
  end subroutine MPP_READ_SIMP_GW
#endif

!DJG -----------------------------------------------------
!   SUBROUTINE READ_SIMP_GW
!DJG -----------------------------------------------------

  subroutine READ_SIMP_GW(IX,JX,IXRT,JXRT,GWSUBBASMSK,gwbasmskfil,&
          gw_strm_msk,numbasns,ch_netrt,AGGFACTRT)
    implicit none
#include <netcdf.inc>

    integer, intent(in)                     :: IX,JX,IXRT,JXRT,AGGFACTRT
    integer, intent(in), dimension(IXRT,JXRT)  :: ch_netrt
    integer, intent(out)                    :: numbasns
    integer, intent(out), dimension(IX,JX)  :: GWSUBBASMSK
    integer, intent(out), dimension(IXRT,JXRT)  :: gw_strm_msk
    character(len=256)                      :: gwbasmskfil
    integer                                 :: i,j,aggfacxrt,aggfacyrt,ixxrt,jyyrt

    numbasns = 0
    gw_strm_msk = -9999

!Open files...
    open(unit=91,file=trim(gwbasmskfil),          &
            form='formatted',status='old')

!Read in sub-basin mask...
    do j=jx,1,-1
          read (91,*) (GWSUBBASMSK(i,j),i=1,ix)
    end do
    close(91)


!Loop through to count number of basins and assign basin indices to chan grid
     do J=1,JX
       do I=1,IX

!Determine max number of basins...(assumes basins are numbered
!   sequentially from 1 to max number of basins...)
        if (GWSUBBASMSK(i,j).gt.numbasns) then
          numbasns = GWSUBBASMSK(i,j)   ! get count of basins...
        end if

!Assign gw basin index values to channel grid...
        do AGGFACYRT=AGGFACTRT-1,0,-1
          do AGGFACXRT=AGGFACTRT-1,0,-1

             IXXRT=I*AGGFACTRT-AGGFACXRT
             JYYRT=J*AGGFACTRT-AGGFACYRT
             IF(ch_netrt(IXXRT,JYYRT).ge.0) then  !If channel grid cell
               gw_strm_msk(IXXRT,JYYRT) = GWSUBBASMSK(i,j)  ! assign coarse grid basn indx to chan grid
             END IF

           end do !AGGFACXRT
         end do !AGGFACYRT

      end do   !I-ix
    end do    !J-jx

#ifdef HYDRO_D
      write(6,*) "numbasns = ", numbasns
#endif

    return

!DJG -----------------------------------------------------
   END SUBROUTINE READ_SIMP_GW
!DJG -----------------------------------------------------

! -------------------------------------------------------- 
! read the static input fields and initial data  
! for 2D OverlandChannel Routing scheme, Qian, 2016.06
! --------------------------------------------------------
  subroutine readRT2d(ix,     jx,    ltype,  ielev,  &
                      elev,   bank,  hlake,  rough,  &
                      sdelev, sdep,  width_obs,      &!!Cong,20190715
                      hycond, poros, head,   aqthick,&
                      filename,topo,irgswfract)
    implicit none

#include <netcdf.inc>

    INTEGER, INTENT(IN)                      :: ix
    INTEGER, INTENT(IN)                      :: jx
    INTEGER, INTENT(INOUT), DIMENSION(ix,jx) :: ltype
    INTEGER, INTENT(INOUT), DIMENSION(ix,jx) :: ielev
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: elev ! hydrology topography, m
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: sdelev ! standard deviation of elevation,-
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: bank
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: width_obs ! river/lake width from obs
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: sdep      ! water storage deposition
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: hlake
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: rough
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: hycond
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: poros
    CHARACTER(LEN=256), INTENT(IN)           :: filename
    REAL,    OPTIONAL, INTENT(INOUT), DIMENSION(ix,jx) :: head
    REAL,    OPTIONAL, INTENT(INOUT), DIMENSION(ix,jx) :: aqthick
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: topo
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: irgswfract

#ifdef MPP_LAND
    INTEGER, ALLOCATABLE,   DIMENSION(:,:)   :: gltype
    INTEGER, ALLOCATABLE,   DIMENSION(:,:)   :: gielev
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: gelev
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: gbank
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: gwidth_obs
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: gsdep
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: ghlake
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: grough
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: ghycond
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: gsdelev
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: gporos
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: ghead
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: gaqthick
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: gihead
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: gtopo  ! topography, m
    REAL,    ALLOCATABLE,   DIMENSION(:,:)   :: girgswfract  !
#endif

    REAL,    DIMENSION(ix,jx) :: ihead !, topo 

    INTEGER :: i, j ! Qian, 2016.08

#ifdef MPP_LAND
    if (my_id .eq. IO_id) then
       allocate(gltype (global_rt_nx,global_rt_ny))
       allocate(gielev (global_rt_nx,global_rt_ny))
       allocate(gelev  (global_rt_nx,global_rt_ny))
       allocate(gbank  (global_rt_nx,global_rt_ny))
       allocate(gwidth_obs  (global_rt_nx,global_rt_ny)) !
       allocate(gsdep  (global_rt_nx,global_rt_ny))      !
       allocate(ghlake (global_rt_nx,global_rt_ny))
       allocate(grough (global_rt_nx,global_rt_ny))
       allocate(ghycond(global_rt_nx,global_rt_ny))
       allocate(sdelev(global_rt_nx,global_rt_ny))
       allocate(gporos (global_rt_nx,global_rt_ny))
       if (present(irgswfract)) allocate(girgswfract (global_rt_nx,global_rt_ny))
       if (present(aqthick)) allocate(gaqthick(global_rt_nx,global_rt_ny)) 
       if (present(head)) allocate(ghead  (global_rt_nx,global_rt_ny))
       if (present(ihead)) allocate(gihead (global_rt_nx,global_rt_ny))
       if (present(gtopo)) allocate(gtopo  (global_rt_nx,global_rt_ny))
    else
       allocate(gltype (1,1))
       allocate(gielev (1,1))
       allocate(gelev  (1,1))
       allocate(gbank  (1,1))
       allocate(gwidth (1,1))
       allocate(gsdep  (1,1))
       allocate(ghlake (1,1))
       allocate(grough (1,1))
       allocate(ghycond(1,1))
       allocate(sdelev(1,1))
       allocate(gporos (1,1))
       if (present(aqthick)) allocate(gaqthick(1,1))
       if (present(head)) allocate(ghead  (1,1))
       if (present(ihead)) allocate(gihead (1,1))
       if (present(topo)) allocate(gtopo  (1,1))
       if (present(irgswfract)) allocate(girgswfract  (1,1))
    endif

#ifndef PARALLELIO
    if (my_id .eq. IO_id) then
#endif
#ifdef HYDRO_D
    print*, "2D OverlandChannelRouting scheme selected, readin ..."
#endif

#endif

    ! land type
    call get2d_int("LTYPE",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   gltype, global_rt_nx, global_rt_ny, &
#else
                   ltype,  ix,        jx,        &
#endif

#else
                   ltype,  ix,        jx,        &
#endif
                   trim(filename) )

    ! sub-basin
    call get2d_int("SUBBASIN",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   gielev, global_rt_nx, global_rt_ny, &
#else
                   ielev,  ix,        jx,        &
#endif

#else
                   ielev,  ix,        jx,        &
#endif
                   trim(filename) )

    ! hydrological topography, filled 
    i = get2d_real("ELEV",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   gelev, global_rt_nx, global_rt_ny, &
#else
                   elev,  ix,        jx,        &
#endif

#else
                   elev,  ix,        jx,        &
#endif
                   trim(filename) )

    ! stream/lake depth
    i = get2d_real("BANK",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   gbank, global_rt_nx, global_rt_ny, &
#else
                   bank,  ix,        jx,        &
#endif

#else
                   bank,  ix,        jx,        &
#endif
                   trim(filename) )

    ! stream/lake width from obs
    i = get2d_real("width_obs",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   gwidth_obs, global_rt_nx, global_rt_ny, &
#else
                   width_obs,  ix,        jx,        &
#endif

#else
                   width_obs,  ix,        jx,        &
#endif
                   trim(filename) )

    ! Storage deposition
    i = get2d_real("d_sfwater",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   gsdep, global_rt_nx, global_rt_ny, &
#else
                   sdep,  ix,        jx,        &
#endif

#else
                   sdep,  ix,        jx,        &
#endif
                   trim(filename) )

    ! initial water level in stream / lake
    i = get2d_real("HLAKE",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   ghlake, global_rt_nx, global_rt_ny, &
#else
                   hlake,  ix,        jx,        &
#endif

#else
                   hlake,  ix,        jx,        &
#endif
                   trim(filename) ) 

    ! standard deviation of elevation
    i = get2d_real("sdelev",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   gsdelev, global_rt_nx, global_rt_ny, &
#else
                   sdelev,  ix,        jx,        &
#endif

#else
                   sdelev,  ix,        jx,        &
#endif
                   trim(filename) )

    ! area irrigated with surface water expressed as percentage of total area
    ! equipped for irrigation
    i = get2d_real("IRSWFRACT",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   girgswfract, global_rt_nx, global_rt_ny, &
#else
                   irgswfract,  ix,        jx,        &
#endif

#else
                   irgswfract,  ix,        jx,        &
#endif
                   trim(filename) )
     irgswfract=irgswfract/100.

    ! Manning roughness 
!    i = get2d_real("ROUGH",  &
!#ifdef MPP_LAND

!#ifndef PARALLELIO
!                   grough, global_rt_nx, global_rt_ny, &
!#else
!                   rough,  ix,        jx,        &
!#endif

!#else
!                   rough,  ix,        jx,        &
!#endif
!                   trim(filename) ) 


    ! acquifer hydraulic conductivity 
    i = get2d_real("cond",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   ghycond, global_rt_nx, global_rt_ny, &
#else
                   hycond,  ix,        jx,        &
#endif

#else
                   hycond,  ix,        jx,        &
#endif
                   trim(filename) )


    ! acquifer porosity
    i = get2d_real("poros",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   gporos, global_rt_nx, global_rt_ny, &
#else
                   poros,  ix,        jx,        &
#endif

#else
                   poros,  ix,        jx,        &
#endif
                   trim(filename) ) 

    ! topography
    if (present(head)) then
    i = get2d_real("HGT_M",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   gtopo, global_rt_nx, global_rt_ny, &
#else
                   topo,  ix,        jx,        &
#endif

#else
                   topo,  ix,        jx,        &
#endif
                   trim(filename) )
    endif

    ! initial groundwater head
    if (present(head)) then
    ! intial groundwater head from global hydrologic modeling
    i = get2d_real("HGWRT",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   gihead, global_rt_nx, global_rt_ny, &
#else
                   ihead,  ix,        jx,        &
#endif

#else
                   ihead,  ix,        jx,        &
#endif
                   trim(filename) )
    endif

    ! aquifer thickness
    if (present(aqthick)) then
        i = get2d_real("THICK", &
#ifdef MPP_LAND
#ifndef PARALLELIO 
                       gaqthick, global_rt_nx, global_rt_ny, &
#else
                       aqthick, ix, jx,  &
#endif
#else
                       aqthick, ix, jx,  &
#endif
                       trim(filename))
    endif

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! terrain height, topo (gtopo), in land surface model
! hydrological topography, elev (gelev), in hydrological model
! Perhaps, they are different. 
! The following steps are to guarantee identical groundwater depths in two models 
#ifdef MPP_LAND
#ifndef PARALLELIO
   do j = 1, global_rt_ny
      do i = 1, global_rt_nx
         !ghycond(i,j) = 1.0   ! = RT_DOMAIN(did)%LKSATRT(IXRT,JXRT)
         !gporos(i,j)  = 1.0   ! = RT_DOMAIN(did)%SMCCSATRT(IXRT,JXRT)
         grough(i,j)  = 0.022  ! unique value from CLHMS, Qian

         if (present(head)) then
            ghead(i,j) = gelev(i,j) - max(0., gtopo(i,j) - gihead(i,j)) 
         endif
      end do
   end do
#else
      write(6,*) "uncompleted parallel readin"
      call hydro_stop("readRT2d")
#endif
#else
   do j = 1, jx
      do i = 1, ix
         !hycond(i,j) = 1.0   ! = RT_DOMAIN(did)%LKSATRT(IXRT,JXRT)
         poros(i,j)  = min(1.0, poros(i,j))   ! = RT_DOMAIN(did)%SMCCSATRT(IXRT,JXRT)
         rough(i,j)  = 0.022  ! unique value from CLHMS, Qian

         if (present(head)) then
            head(i,j) = elev(i,j) - max(0., topo(i,j) - ihead(i,j))
!C.Jiang,2020
!Simple initial groundwater head without Global dataset
!            head(i,j) = elev(i,j) - 20.0
            head(i,j) = ihead(i,j)
         endif
      end do
   end do
#endif
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef MPP_LAND
#ifndef PARALLELIO
!    gltype(1, :) = 2			! commmented by Qian
!    gltype(:, 1) = 2
!    gltype(global_rt_nx, :) = 2
!    gltype(:, global_rt_ny) = 2
#else 
     ! to do parallel io for  
#endif
#else
!       ltype(1,:) = 2
!       ltype(:,1) = 2
!       ltype(ix,:)= 2
!       ltype(:,jx)= 2
#endif

#ifdef MPP_LAND
#ifndef PARALLELIO
    endif
    call decompose_rt_int (gltype,  ltype,  global_rt_nx, global_rt_ny, ix, jx)
    call decompose_rt_int (gielev,  ielev,  global_rt_nx, global_rt_ny, ix, jx)
    call decompose_rt_real(gelev,   elev,   global_rt_nx, global_rt_ny, ix, jx)
    call decompose_rt_real(gbank,   bank,   global_rt_nx, global_rt_ny, ix, jx)
    call decompose_rt_real(gwidth_obs, width_obs,   global_rt_nx, global_rt_ny, ix, jx)
    call decompose_rt_real(gsdep,   sdep,   global_rt_nx, global_rt_ny, ix, jx)
    call decompose_rt_real(ghlake,  hlake,  global_rt_nx, global_rt_ny, ix, jx)
    call decompose_rt_real(grough,  rough,  global_rt_nx, global_rt_ny, ix, jx)
    call decompose_rt_real(ghycond, hycond, global_rt_nx, global_rt_ny, ix, jx)
    call decompose_rt_real(gporos,  poros,  global_rt_nx, global_rt_ny, ix, jx)
    call decompose_rt_real(gsdelev,sdelev,  global_rt_nx, global_rt_ny, ix, jx)
    if (present(head)) call decompose_rt_real(ghead, head, global_rt_nx, global_rt_ny, ix, jx)
    if (present(aqthick))  call decompose_rt_real(gaqthick, aqthick, global_rt_nx, global_rt_ny, ix, jx)
    deallocate(gltype, gielev, gbank, gwidth_obs, gsdep, ghlake, grough, ghycond, gporos, gsdelev)
    if (present(head)) deallocate(ghead, gihead, gtopo)
    if (present(aqthick)) deallocate(gaqthick)
#endif
#endif 

    RETURN

  end subroutine readRT2d

! -------------------------------------------------------- 

  ! BF read the static input fields needed for the 2D GW scheme
  subroutine readGW2d(ix,      jx,  hc,    head,    &
                      aqthick, por, ltype, elev,    &
                      ihShift, filename )
  implicit none
#include <netcdf.inc>
  integer, intent(in) :: ix, jx
  real,    intent(in) :: ihShift
  integer, dimension(ix,jx), intent(inout) :: ltype
  real,    dimension(ix,jx), intent(inout) :: hc, head, aqthick, por
  real,    dimension(ix,jx), intent(inout) :: elev                   ! Qian, 2016.08
  CHARACTER(LEN=256), INTENT(IN)           :: filename 

#ifdef MPP_LAND
  integer, dimension(:,:), allocatable ::  gLtype
  real,    dimension(:,:), allocatable ::  gHC, gHEAD, gaqthick, gPOR
  real,    dimension(:,:), allocatable ::  gELEV
  real,    dimension(:,:), allocatable ::  gIHEAD, gTOPO
#endif

  real,    dimension(ix,jx) :: topo, ihead
  integer :: i
  integer :: j  

#ifdef MPP_LAND
  if(my_id .eq. IO_id) then
      allocate(gHC(global_rt_nx, global_rt_ny))
      allocate(gIHEAD(global_rt_nx, global_rt_ny))
      allocate(gAQTHICK(global_rt_nx, global_rt_ny))
      allocate(gPOR(global_rt_nx, global_rt_ny))
      allocate(gLtype(global_rt_nx, global_rt_ny))
      allocate(gELEV(global_rt_nx, global_rt_ny))
      allocate(gHEAD(global_rt_nx, global_rt_ny))
      allocate(gTOPO(global_rt_nx, global_rt_ny))
  else
      allocate(gHC(1, 1))
      allocate(gIHEAD(1, 1))
      allocate(gAQTHICK(1, 1))
      allocate(gPOR(1, 1))
      allocate(gLtype(1, 1))
      allocate(gELEV(1, 1))
      allocate(gHEAD(1, 1))
      allocate(gTOPO(1, 1))
  endif
 
#ifndef PARALLELIO 
  if(my_id .eq. IO_id) then
#endif
#ifdef HYDRO_D
  print*, "2D GW-Scheme selected, retrieving files from gwhires.nc ..."
#endif
#endif

        ! comments by Qian, 2016.08
        ! hydraulic conductivity
!        i = get2d_real("HC", &
!#ifdef MPP_LAND
!#ifndef PARALLELIO 
!                       gHC, global_rt_nx, global_rt_ny,  &
!#else
!                       hc, ix, jx,  &
!#endif
!#else
!                       hc, ix, jx,  &
!#endif
!                       trim("./gwhires.nc"))

        ! initial head
!        i = get2d_real("IHEAD", &
!#ifdef MPP_LAND
!                       gIHEAD, global_rt_nx, global_rt_ny, &
!#else
!                       ihead,  ix, jx, &
!#endif
!                       trim("./gwhires.nc"))

        ! topotgraphy
        i = get2d_real("HGT_M", &
#ifdef MPP_LAND
#ifndef PARALLELIO
                       gTOPO, global_rt_nx, global_rt_ny, &
#else
                       topo,  ix, jx, &
#endif
#else
                       topo,  ix, jx, &
#endif
                       trim(filename))

        ! initial groundwater head from global modelling
        i = get2d_real("HGWRT", &
#ifdef MPP_LAND
#ifndef PARALLELIO
                       gIHEAD, global_rt_nx, global_rt_ny, &
#else
                       ihead,  ix, jx, &
#endif
#else
                       ihead,  ix, jx, &
#endif
                       trim(filename))

        ! aquifer thickness          
        i = get2d_real("THICK", &
#ifdef MPP_LAND
#ifndef PARALLELIO 
                       gAQTHICK, global_rt_nx, global_rt_ny, &
#else
                       aqthick, ix, jx,  &
#endif
#else
                       aqthick, ix, jx,  &
#endif
                       trim(filename))
                       
!	! aquifer porosity
!        i = get2d_real("POR", &
!#ifdef MPP_LAND
!#ifndef PARALLELIO 
!                       gPOR, global_rt_nx, global_rt_ny, &
!#else
!                       por, ix, jx,  &
!#endif
!#else
!                       por, ix, jx,  &
!#endif
!                       trim("./gwhires.nc"))

    ! land type, Qian, 2016.08
    ! ltype(i,j) = 0, ocean point
    ! ltype(i,j) = 1, land point in the watershed
    ! ltype(i,j) = 2, land point out of the watershed
        call get2d_int("LTYPE", &
#ifdef MPP_LAND
#ifndef PARALLELIO 
                       gLtype, global_rt_nx, global_rt_ny, &
#else
                       ltype, ix, jx, &
#endif
#else
                       ltype, ix, jx,  &
#endif
                       trim(filename))


    ! hydrological topography 
    ! C.Jiang,2020.11 topography 
    i = get2d_real("ELEV",  &
#ifdef MPP_LAND

#ifndef PARALLELIO
                   gelev, global_rt_nx, global_rt_ny, &
#else
                   elev,  ix,        jx,        &
#endif

#else
                   elev,  ix,        jx,        &
#endif
                   trim(filename) )

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! terrain height, topo (gtopo), in land surface model
! hydrological topography, elev (gelev), in hydrological model
! Perhaps, they are different. 
! The following steps are to guarantee identical groundwater depths in two models
#ifdef MPP_LAND
#ifndef PARALLELIO
   do j = 1, global_rt_ny
      do i = 1, global_rt_nx
         gHC(i,j) = 1.0
         gPOR(i,j)  = 1.0
         ghead(i,j) = gelev(i,j) - max(0., gtopo(i,j) - gihead(i,j))
      end do
   end do
#else
      write(6,*) "uncompleted parallel readin"
      call hydro_stop("readGW2d")
#endif
#else
   do j = 1, jx
      do i = 1, ix
         hc(i,j) = 1.0
         por(i,j)  = 1.0
         head(i,j) = elev(i,j) - max(0., topo(i,j) - ihead(i,j))

        !CJiang,2020 
        !Simple initial groundwater head without Global dataset
!         head(i,j) = elev(i,j) - 20.0
         head(i,j) = ihead(i,j)
      end do
   end do


#endif
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
#ifdef MPP_LAND
#ifndef PARALLELIO 
!       gLtype(1,:) = 2		! commented by Qian
!       gLtype(:,1) = 2
!       gLtype(global_rt_nx,:) = 2
!       gLtype(:,global_rt_ny) = 2 
#else
! BF TODO parallel io for gw ltype
#endif
#else
! deleted by Qian, 2016.08
!       ltype(1,:) = 2
!       ltype(:,1) = 2
!       ltype(ix,:)= 2
!       ltype(:,jx)= 2
#endif

#ifdef MPP_LAND  
#ifndef PARALLELIO 
  endif
     call decompose_rt_int (gLtype, ltype, global_rt_nx, global_rt_ny, ix, jx)
     call decompose_rt_real(gHC,hc,global_rt_nx, global_rt_ny, ix, jx)
     call decompose_rt_real(gHEAD,head,global_rt_nx, global_rt_ny, ix, jx)
     call decompose_rt_real(gAQTHICK,aqthick,global_rt_nx, global_rt_ny, ix, jx)
     call decompose_rt_real(gPOR,por,global_rt_nx, global_rt_ny, ix, jx)
     deallocate(gHC, gIHEAD, gAQTHICK, gPOR, gHEAD, gTOPO)
#endif
#endif
    
  ihead = ihead + ihShift
  
!  where(ltype .eq. 0) 
!    hc = 0.
!   por = 10.E21    
!    por = 1.0       ! Qian
!  end where
  
  !bftodo: make filename accessible in namelist
  return
  end subroutine readGW2d
  !BF
 
! --------------------------------------------------------
! read crop water requirement data (met by irrigation)
! for 2D OverlandChannel Routing water use scheme, Jiang, 2020.12
! --------------------------------------------------------
  subroutine readIRG2d(ix,        jx,  irgM1,  irgM2,  &
                      irgM3,  irgM4,  irgM5,  irgM6,  &
                      irgM7,  irgM8,  irgM9,  irgM10, &
                      irgM11, irgM12, filename )
    implicit none

#include <netcdf.inc>

    INTEGER, INTENT(IN)                      :: ix
    INTEGER, INTENT(IN)                      :: jx
    CHARACTER(LEN=256), INTENT(IN)           :: filename
    REAL,    INTENT(INOUT), DIMENSION(ix,jx) :: irgM1, irgM2, irgM3, irgM4, irgM5, irgM6, & 
                                                irgM7, irgM8, irgM9, irgM10,irgM11,irgM12 
    integer :: i

#ifdef MPP_LAND
    real,    dimension(:,:), allocatable ::       girgM1, girgM2, girgM3, girgM4, girgM5, girgM6, &
                                                girgM7, girgM8, girgM9, girgM10,girgM11,girgM12
#endif
    
#ifdef MPP_LAND
  if(my_id .eq. IO_id) then
      allocate(girgM1(global_rt_nx, global_rt_ny))
      allocate(girgM2(global_rt_nx, global_rt_ny))
      allocate(girgM3(global_rt_nx, global_rt_ny))
      allocate(girgM4(global_rt_nx, global_rt_ny))
      allocate(girgM5(global_rt_nx, global_rt_ny))
      allocate(girgM6(global_rt_nx, global_rt_ny))
      allocate(girgM7(global_rt_nx, global_rt_ny))
      allocate(girgM8(global_rt_nx, global_rt_ny))
      allocate(girgM9(global_rt_nx, global_rt_ny))
      allocate(girgM10(global_rt_nx, global_rt_ny))
      allocate(girgM11(global_rt_nx, global_rt_ny))
      allocate(girgM12(global_rt_nx, global_rt_ny))
  else
      allocate(girgM1(1, 1))
      allocate(girgM2(1, 1))
      allocate(girgM3(1, 1))
      allocate(girgM4(1, 1))
      allocate(girgM5(1, 1))
      allocate(girgM6(1, 1))
      allocate(girgM7(1, 1))
      allocate(girgM8(1, 1))
      allocate(girgM9(1, 1))
      allocate(girgM10(1, 1))
      allocate(girgM11(1, 1))
      allocate(girgM12(1, 1))
  endif


#ifndef PARALLELIO
    if (my_id .eq. IO_id) then
#endif
#ifdef HYDRO_D
    print*, "2D OverlandChannelRouting Water Use Module, readin ..."
#endif

#endif

    ! Crop water requirement met by irrigation in Janunary
    i = get2d_real("BW_irr_M01",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM1, global_rt_nx, global_rt_ny, &
#else
                   irgM1,  ix,        jx,        &
#endif

#else
                   irgM1,  ix,        jx,        &
#endif
                   trim(filename) )

    ! Crop water requirement met by irrigation in February
    i = get2d_real("BW_irr_M02",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM2, global_rt_nx, global_rt_ny, &
#else
                   irgM2,  ix,        jx,        &
#endif

#else
                   irgM2,  ix,        jx,        &
#endif
                   trim(filename) )

    ! Crop water requirement met by irrigation in March
    i = get2d_real("BW_irr_M03",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM3, global_rt_nx, global_rt_ny, &
#else
                   irgM3,  ix,        jx,        &
#endif

#else
                   irgM3,  ix,        jx,        &
#endif
                   trim(filename) )

    ! Crop water requirement met by irrigation in April
    i = get2d_real("BW_irr_M04",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM4, global_rt_nx, global_rt_ny, &
#else
                   irgM4,  ix,        jx,        &
#endif

#else
                   irgM4,  ix,        jx,        &
#endif
                   trim(filename) )

    ! Crop water requirement met by irrigation in May
    i = get2d_real("BW_irr_M05",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM5, global_rt_nx, global_rt_ny, &
#else
                   irgM5,  ix,        jx,        &
#endif

#else
                   irgM5,  ix,        jx,        &
#endif
                   trim(filename) )

    ! Crop water requirement met by irrigation in June
    i = get2d_real("BW_irr_M06",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM6, global_rt_nx, global_rt_ny, &
#else
                   irgM6,  ix,        jx,        &
#endif

#else
                   irgM6,  ix,        jx,        &
#endif
                   trim(filename) )

    ! Crop water requirement met by irrigation in July
    i = get2d_real("BW_irr_M07",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM7, global_rt_nx, global_rt_ny, &
#else
                   irgM7,  ix,        jx,        &
#endif

#else
                   irgM7,  ix,        jx,        &
#endif
                   trim(filename) )

    ! Crop water requirement met by irrigation in August
    i = get2d_real("BW_irr_M08",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM8, global_rt_nx, global_rt_ny, &
#else
                   irgM8,  ix,        jx,        &
#endif

#else
                   irgM8,  ix,        jx,        &
#endif
                   trim(filename) )


    ! Crop water requirement met by irrigation in September
    i = get2d_real("BW_irr_M09",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM9, global_rt_nx, global_rt_ny, &
#else
                   irgM9,  ix,        jx,        &
#endif

#else
                   irgM9,  ix,        jx,        &
#endif
                   trim(filename) )

    ! Crop water requirement met by irrigation in Octorber
    i = get2d_real("BW_irr_M10",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM10, global_rt_nx, global_rt_ny, &
#else
                   irgM10,  ix,        jx,        &
#endif

#else
                   irgM10,  ix,        jx,        &
#endif
                   trim(filename) )

    ! Crop water requirement met by irrigation in November
    i = get2d_real("BW_irr_M11",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM11, global_rt_nx, global_rt_ny, &
#else
                   irgM11,  ix,        jx,        &
#endif

#else
                   irgM11,  ix,        jx,        &
#endif
                   trim(filename) )


    ! Crop water requirement met by irrigation in December
    i = get2d_real("BW_irr_M12",  &
#ifdef MPP_LAND
#ifndef PARALLELIO
                   girgM12, global_rt_nx, global_rt_ny, &
#else
                   irgM12,  ix,        jx,        &
#endif

#else
                   irgM12,  ix,        jx,        &
#endif
                   trim(filename) )

#ifdef MPP_LAND
#ifndef PARALLELIO
  endif
!     call decompose_rt_int (gLtype, ltype, global_rt_nx, global_rt_ny, ix, jx)
!     call decompose_rt_real(gHC,hc,global_rt_nx, global_rt_ny, ix, jx)
!     deallocate(gHC, gIHEAD, gAQTHICK, gPOR, gHEAD, gTOPO)
#endif
#endif


  return

  end subroutine readIRG2d


  subroutine output_rt(igrid, split_output_count, ixrt, jxrt, nsoil, &
       startdate, date, QSUBRT,ZWATTABLRT,SMCRT,SUB_RESID,       &
       q_sfcflx_x,q_sfcflx_y,soxrt,soyrt,QSTRMVOLRT,SFCHEADSUBRT, &
       geo_finegrid_flnm,dt,sldpth,LATVAL,LONVAL,dist,CHRTOUT_GRID,  &
       QBDRYRT)

!output the routing variables over routing grid.
    implicit none
#include <netcdf.inc>

    integer,                                  intent(in) :: igrid
    integer,                                  intent(in) :: split_output_count
    integer,                                  intent(in) :: ixrt,jxrt
    real,                                     intent(in) :: dt
    real,                                     intent(in) :: dist(ixrt,jxrt,9)
    integer,                                  intent(in) :: nsoil
    integer,                                  intent(in) :: CHRTOUT_GRID
    character(len=*),                         intent(in) :: startdate
    character(len=*),                         intent(in) :: date
    character(len=*),          intent(in)                :: geo_finegrid_flnm
    real,             dimension(nsoil),       intent(in) :: sldpth
    real, allocatable, DIMENSION(:,:)                   :: xdumd  !-- decimated variable
    real*8, allocatable, DIMENSION(:)                   :: xcoord_d
    real*8, allocatable, DIMENSION(:)                   :: ycoord_d, ycoord

    integer, save :: ncid,ncstatic
    integer, save :: output_count
    real,    dimension(nsoil) :: asldpth

    integer :: dimid_ix, dimid_jx, dimid_times, dimid_datelen, varid, n
    integer :: iret, dimid_soil, i,j,ii,jj
    character(len=256) :: output_flnm
    character(len=19)  :: date19
    character(len=32)  :: convention
    character(len=34)  :: sec_since_date

    character(len=30)  :: soilm

    real                                :: long_cm,lat_po,fe,fn, chan_in
    real, dimension(2)                  :: sp

    real, dimension(ixrt,jxrt) :: xdum,QSUBRT,ZWATTABLRT,SUB_RESID
    real, dimension(ixrt,jxrt) :: q_sfcflx_x,q_sfcflx_y
    real, dimension(ixrt,jxrt) :: QSTRMVOLRT
    real, dimension(ixrt,jxrt) :: SFCHEADSUBRT
    real, dimension(ixrt,jxrt) :: soxrt,soyrt
    real, dimension(ixrt,jxrt) :: LATVAL,LONVAL, QBDRYRT
    real, dimension(ixrt,jxrt,nsoil) :: SMCRT

    character(len=2) :: strTmp

    integer :: seconds_since, decimation, ixrtd,jxrtd, hires_flag
    sec_since_date = 'seconds since '//date(1:4)//'-'//date(6:7)//'-'//date(9:10)//' '//date(12:13)//':'//date(15:16)//' UTC'
    seconds_since = int(dt)*output_count

    decimation = 1 !-- decimation factor
#ifdef MPP_LAND
    ixrtd = int(global_rt_nx/decimation)
    jxrtd = int(global_rt_ny/decimation)
#else
    ixrtd = int(ixrt/decimation)
    jxrtd = int(jxrt/decimation)
#endif

#ifdef MPP_LAND
    if(my_id .eq. io_id) then
#endif
       allocate(xdumd(ixrtd,jxrtd))
       allocate(xcoord_d(ixrtd))
       allocate(ycoord_d(jxrtd))
       allocate(ycoord(jxrtd))

       xdumd = -999
       xcoord_d = -999 
       ycoord_d = -999
       ycoord = -999
#ifdef MPP_LAND
    else
       allocate(xdumd(1,1))
       allocate(xcoord_d(1))
       allocate(ycoord_d(1))
       allocate(ycoord(1))
    endif
#endif
    ii = 0

!DJG Dump timeseries for channel inflow accum. for calibration...(8/28/09)
    chan_in = 0.0
    do j=1,jxrt
      do i=1,ixrt
        chan_in=chan_in+QSTRMVOLRT(I,J)/1000.0*(dist(i,j,9))  !(units m^3)
      enddo
    enddo
#ifdef MPP_LAND
      call sum_real1(chan_in)
#endif
#ifdef MPP_LAND
    if(my_id .eq. io_id) then
#endif
       open (unit=46,file='qstrmvolrt_accum.txt',form='formatted',&
             status='unknown',position='append')
        write (46,713) chan_in
       close (46)
#ifdef MPP_LAND
    endif
#endif
713 FORMAT (F20.7)
!    return
!DJG end dump of channel inflow for calibration....

    if (CHRTOUT_GRID.eq.0) return  ! return if hires flag eq 1, if =2 output full grid

    if (output_count == 0) then  

   !-- Open the  finemesh static files to obtain projection information
#ifdef HYDRO_D
      write(*,'("geo_finegrid_flnm: ''", A, "''")') trim(geo_finegrid_flnm)
#endif

#ifdef MPP_LAND
   if(my_id .eq. io_id) then
#endif
      iret = nf_open(trim(geo_finegrid_flnm), NF_NOWRITE, ncstatic)
#ifdef MPP_LAND
   endif
   call mpp_land_bcast_int1(iret)
#endif

      if (iret /= 0) then
         write(*,'("Problem opening geo_finegrid file: ''", A, "''")') &
         trim(geo_finegrid_flnm)
         write(*,*) "HIRES_OUTPUT will not be georeferenced..."
        hires_flag = 0
      else
        hires_flag = 1
      endif

#ifdef MPP_LAND
   if(my_id .eq. io_id) then
#endif

     if(hires_flag.eq.1) then !if/then hires_georef
      ! Get Latitude (X)
      iret = NF_INQ_VARID(ncstatic,'x',varid)
      if(iret .eq. 0) iret = NF_GET_VAR_DOUBLE(ncstatic, varid, xcoord_d)
      ! Get Longitude (Y)
      iret = NF_INQ_VARID(ncstatic,'y',varid)
      if(iret .eq. 0) iret = NF_GET_VAR_DOUBLE(ncstatic, varid, ycoord)
     else
      ycoord_d = 0.
      xcoord_d = 0.
     end if  !endif hires_georef 

     jj = 0
#ifdef MPP_LAND
     do j=global_rt_ny,1,-1*decimation
#else
     do j=jxrt,1,-1*decimation
#endif
        jj = jj+1
        if (jj<= jxrtd) then
         ycoord_d(jj) = ycoord(j)
        endif
     enddo

       

     if(hires_flag.eq.1) then !if/then hires_georef
      ! Get projection information from finegrid netcdf file
      iret = NF_INQ_VARID(ncstatic,'lambert_conformal_conic',varid)
      if(iret .eq. 0) iret = NF_GET_ATT_REAL(ncstatic, varid, 'longitude_of_central_meridian', long_cm)  !-- read it from the static file
      iret = NF_GET_ATT_REAL(ncstatic, varid, 'latitude_of_projection_origin', lat_po)  !-- read it from the static file
      iret = NF_GET_ATT_REAL(ncstatic, varid, 'false_easting', fe)  !-- read it from the static file
      iret = NF_GET_ATT_REAL(ncstatic, varid, 'false_northing', fn)  !-- read it from the static file
      iret = NF_GET_ATT_REAL(ncstatic, varid, 'standard_parallel', sp)  !-- read it from the static file
     end if  !endif hires_georef 
      iret = nf_close(ncstatic)


!-- create the fine grid routing file
       write(output_flnm, '(A12,".RTOUT_DOMAIN",I1)') date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16), igrid
#ifdef HYDRO_D
       print*, 'output_flnm = "'//trim(output_flnm)//'"'
#endif
#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
       write(6,*) "using large netcdf file for RTOUT_DOMAIN"
       iret = nf_create(trim(output_flnm), IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
#else
       write(6,*) "using normal netcdf file for RTOUT_DOMAIN"
#ifdef PARALLELIO
      iret = nf90_open_par(trim(geo_finegrid_flnm), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#else
       iret = nf_create(trim(output_flnm), NF_CLOBBER, ncid)
#endif
#endif
#ifdef HYDRO_D
       if (iret /= 0) then
         print*, "Problem nf_create"
         call hydro_finish("output_rt")
       endif
#endif

       iret = nf_def_dim(ncid, "time", NF_UNLIMITED, dimid_times)
       iret = nf_def_dim(ncid, "x", ixrtd, dimid_ix)  !-- make a decimated grid
       iret = nf_def_dim(ncid, "y", jxrtd, dimid_jx)
       iret = nf_def_dim(ncid, "depth", nsoil, dimid_soil)  !-- 3-d soils

!--- define variables
!     !- time definition, timeObs
         iret = nf_def_var(ncid,"time",NF_INT, 1, (/dimid_times/), varid)
         iret = nf_put_att_text(ncid,varid,'units',34,sec_since_date)

       !- x-coordinate in cartesian system
        iret = nf_def_var(ncid,"x",NF_DOUBLE, 1, (/dimid_ix/), varid)
        iret = nf_put_att_text(ncid,varid,'long_name',26,'x coordinate of projection')
        iret = nf_put_att_text(ncid,varid,'standard_name',23,'projection_x_coordinate')
        iret = nf_put_att_text(ncid,varid,'units',5,'Meter')

       !- y-coordinate in cartesian ssystem
          iret = nf_def_var(ncid,"y",NF_DOUBLE, 1, (/dimid_jx/), varid)
          iret = nf_put_att_text(ncid,varid,'long_name',26,'y coordinate of projection')
          iret = nf_put_att_text(ncid,varid,'standard_name',23,'projection_y_coordinate')
          iret = nf_put_att_text(ncid,varid,'units',5,'Meter')

       !- LATITUDE
        iret = nf_def_var(ncid,"LATITUDE",NF_FLOAT, 2, (/dimid_ix,dimid_jx/), varid)
        iret = nf_put_att_text(ncid,varid,'long_name',8,'LATITUDE')
        iret = nf_put_att_text(ncid,varid,'standard_name',8,'LATITUDE')
        iret = nf_put_att_text(ncid,varid,'units',5,'deg North')

       !- LONGITUDE
          iret = nf_def_var(ncid,"LONGITUDE",NF_FLOAT, 2, (/dimid_ix,dimid_jx/), varid)
          iret = nf_put_att_text(ncid,varid,'long_name',9,'LONGITUDE')
          iret = nf_put_att_text(ncid,varid,'standard_name',9,'LONGITUDE')
          iret = nf_put_att_text(ncid,varid,'units',5,'deg east')

       !-- z-level is soil
        iret = nf_def_var(ncid,"depth", NF_FLOAT, 1, (/dimid_soil/),varid)
        iret = nf_put_att_text(ncid,varid,'units',2,'cm')
        iret = nf_put_att_text(ncid,varid,'long_name',19,'depth of soil layer')

         do n = 1, NSOIL
             write(strTmp,'(I2)') n
             iret = nf_def_var(ncid,  "SOIL_M"//trim(strTmp),  NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
         end do
            iret = nf_put_att_text(ncid,varid,'units',7,'m^3/m^3')
            iret = nf_put_att_text(ncid,varid,'description',16,'moisture content')
            iret = nf_put_att_text(ncid,varid,'long_name',26,soilm)
!           iret = nf_put_att_text(ncid,varid,'coordinates',5,'x y z')
            iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
            iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

!      iret = nf_def_var(ncid,"ESNOW2D",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)

!       iret = nf_def_var(ncid,"QSUBRT",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
!       iret = nf_put_att_text(ncid,varid,'units',6,'m3 s-1')
!       iret = nf_put_att_text(ncid,varid,'long_name',15,'subsurface flow')
!       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
!       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
!       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

          iret = nf_def_var(ncid,"ZWATTABLRT",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
          iret = nf_put_att_text(ncid,varid,'units',1,'m')
          iret = nf_put_att_text(ncid,varid,'long_name',17,'water table depth')
          iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
          iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
          iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

!       iret = nf_def_var(ncid,"Q_SFCFLX_X",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
!       iret = nf_put_att_text(ncid,varid,'units',6,'m3 s-1')
!       iret = nf_put_att_text(ncid,varid,'long_name',14,'surface flux x')
!       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
!       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
!       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

!       iret = nf_def_var(ncid,"Q_SFCFLX_Y",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
!       iret = nf_put_att_text(ncid,varid,'units',6,'m3 s-1')
!       iret = nf_put_att_text(ncid,varid,'long_name',14,'surface flux y')
!       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
!       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
!       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

       iret = nf_def_var(ncid,"QSTRMVOLRT",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'mm')
       iret = nf_put_att_text(ncid,varid,'long_name',14,'channel inflow')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

       iret = nf_def_var(ncid,"SFCHEADSUBRT",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'mm')
       iret = nf_put_att_text(ncid,varid,'long_name',12,'surface head')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

!      iret = nf_def_var(ncid,"SOXRT",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
!      iret = nf_put_att_text(ncid,varid,'units',1,'1')
!      iret = nf_put_att_text(ncid,varid,'long_name',7,'slope x')
!      iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
!      iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
!      iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

!      iret = nf_def_var(ncid,"SOYRT",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
!      iret = nf_put_att_text(ncid,varid,'units',1,'1')
!      iret = nf_put_att_text(ncid,varid,'long_name',7,'slope 7')
!      iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
!      iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
!      iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

!       iret = nf_def_var(ncid,"SUB_RESID",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)

       iret = nf_def_var(ncid,"QBDRYRT",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'mm')
       iret = nf_put_att_text(ncid,varid,'long_name',70, &
          'accumulated value of the boundary flux, + into domain, - out of domain')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

!-- place projection information
     if(hires_flag.eq.1) then !if/then hires_georef
      iret = nf_def_var(ncid,"lambert_conformal_conic",NF_INT,0, 0,varid)
      iret = nf_put_att_text(ncid,varid,'grid_mapping_name',23,'lambert_conformal_conic')
      iret = nf_put_att_real(ncid,varid,'longitude_of_central_meridian',NF_FLOAT,1,long_cm)
      iret = nf_put_att_real(ncid,varid,'latitude_of_projection_origin',NF_FLOAT,1,lat_po)
      iret = nf_put_att_real(ncid,varid,'false_easting',NF_FLOAT,1,fe)
      iret = nf_put_att_real(ncid,varid,'false_northing',NF_FLOAT,1,fn)
      iret = nf_put_att_real(ncid,varid,'standard_parallel',NF_FLOAT,2,sp)
     end if   !endif hires_georef

!      iret = nf_def_var(ncid,"Date",   NF_CHAR,  2, (/dimid_datelen,dimid_times/),     varid)

      date19(1:19) = "0000-00-00_00:00:00"
      date19(1:len_trim(startdate)) = startdate
      convention(1:32) = "CF-1.0"
      iret = nf_put_att_real(ncid, NF_GLOBAL, "missing_value", NF_FLOAT, 1, -9E15)
      iret = nf_put_att_text(ncid, NF_GLOBAL, "Conventions",6, convention)
      iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_start",19, date19)
      iret = nf_put_att_int(ncid,NF_GLOBAL,"output_decimation_factor",NF_INT, 1,decimation)

       ! iret = nf_redef(ncid)
       date19(1:len_trim(date)) = date
       iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_end", 19, date19)
       ! iret = nf_enddef(ncid)

      iret = nf_enddef(ncid)

!!-- write latitude and longitude locations
         iret = nf_inq_varid(ncid,"x", varid)
         iret = nf_put_vara_double(ncid, varid, (/1/), (/ixrtd/), xcoord_d) !-- 1-d array

         iret = nf_inq_varid(ncid,"y", varid)
         iret = nf_put_vara_double(ncid, varid, (/1/), (/jxrtd/), ycoord_d) !-- 1-d array

#ifdef MPP_LAND
    endif

        call write_IO_rt_real(LATVAL,xdumd)
    if( my_id .eq. io_id) then
#else
        xdumd = LATVAL
#endif
        iret = nf_inq_varid(ncid,"LATITUDE", varid)
        iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ixrtd,jxrtd/), xdumd)


#ifdef MPP_LAND
    endif   !!! end if block of my_id .eq. io_id

        call write_IO_rt_real(LONVAL,xdumd)

    if( my_id .eq. io_id) then
#else
        xdumd = LONVAL
#endif
        iret = nf_inq_varid(ncid,"LONGITUDE", varid)
        iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ixrtd,jxrtd/), xdumd)

#ifdef MPP_LAND
    endif

    if( my_id .eq. io_id) then
#endif






       do n = 1,nsoil
        if(n == 1) then
         asldpth(n) = -sldpth(n)
        else
         asldpth(n) = asldpth(n-1) - sldpth(n)
        endif
       enddo

       iret = nf_inq_varid(ncid,"depth", varid)
       iret = nf_put_vara_real(ncid, varid, (/1/), (/nsoil/), asldpth)
!yw       iret = nf_close(ncstatic)
#ifdef MPP_LAND
    endif  ! end of my_id .eq. io_id
#endif
   endif !!! end of if block output_count == 0
    output_count = output_count + 1

!!-- time
#ifdef MPP_LAND
   if(my_id .eq. io_id) then
#endif
        iret = nf_inq_varid(ncid,"time", varid)
        iret = nf_put_vara_int(ncid, varid, (/output_count/), (/1/), seconds_since)
#ifdef MPP_LAND
   endif
#endif

!-- 3-d soils
     do n = 1, nsoil
#ifdef MPP_LAND
          call write_IO_rt_real(smcrt(:,:,n),xdumd)
#else
          xdumd(:,:) = smcrt(:,:,n)
#endif
#ifdef MPP_LAND
    if(my_id .eq. io_id) then
#endif
          write(strTmp,'(I2)') n
          iret = nf_inq_varid(ncid,  "SOIL_M"//trim(strTmp), varid)
          iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
#ifdef MPP_LAND
    endif
#endif
    enddo !-n soils



#ifdef MPP_LAND
          call write_IO_rt_real(ZWATTABLRT,xdumd)
#else
          xdumd(:,:) = ZWATTABLRT(:,:)
#endif
#ifdef MPP_LAND
    if(my_id .eq. io_id) then
#endif
        iret = nf_inq_varid(ncid,  "ZWATTABLRT", varid)
        iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
#ifdef MPP_LAND
     endif
#endif



#ifdef MPP_LAND
          call write_IO_rt_real(QBDRYRT,xdumd)
#else
          xdumd(:,:) = QBDRYRT(:,:)
#endif
#ifdef MPP_LAND
    if(my_id .eq. io_id) then
#endif
     iret = nf_inq_varid(ncid,  "QBDRYRT", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
#ifdef MPP_LAND
     endif
#endif


#ifdef MPP_LAND
          call write_IO_rt_real(QSTRMVOLRT,xdumd)
#else
          xdumd(:,:) = QSTRMVOLRT(:,:)
#endif
#ifdef MPP_LAND
    if(my_id .eq. io_id) then
#endif
     iret = nf_inq_varid(ncid,  "QSTRMVOLRT", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
#ifdef MPP_LAND
     endif
#endif
 

#ifdef MPP_LAND
          call write_IO_rt_real(SFCHEADSUBRT,xdumd)
#else
          xdumd(:,:) = SFCHEADSUBRT(:,:)
#endif
#ifdef MPP_LAND
    if(my_id .eq. io_id) then
#endif
     iret = nf_inq_varid(ncid,  "SFCHEADSUBRT", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
#ifdef MPP_LAND
     endif
#endif
 
#ifdef MPP_LAND
    if(my_id .eq. io_id) then
#endif


!yw      iret = nf_sync(ncid)
      if (output_count == split_output_count) then
        output_count = 0
        iret = nf_close(ncid)
      endif
#ifdef MPP_LAND
     endif
     call mpp_land_bcast_int1(output_count)
#endif

     deallocate(xdumd)
     deallocate(xcoord_d)
     deallocate(ycoord_d)
     deallocate(ycoord)
    
#ifdef HYDRO_D 
     write(6,*) "end of output_rt" 
#endif

  end subroutine output_rt


!BF output section for gw2d model
!bftodo: clean up an customize for GW usage

  subroutine output_gw_spinup(igrid, split_output_count, ixrt, jxrt, &
       startdate, date, HEAD, convgw, excess, &
       geo_finegrid_flnm,dt,LATVAL,LONVAL,dist,output_gw,outname)

#ifdef MPP_LAND
       USE module_mpp_land
#endif
!output the routing variables over routing grid.
    implicit none
#include <netcdf.inc>

    integer,                                  intent(in) :: igrid
    integer,                                  intent(in) :: split_output_count
    integer,                                  intent(in) :: ixrt,jxrt
    real,                                     intent(in) :: dt
    real,                                     intent(in) :: dist(ixrt,jxrt,9)
    integer,                                  intent(in) :: output_gw
    character(len=*),                         intent(in) :: outname
    character(len=*),                         intent(in) :: startdate
    character(len=*),                         intent(in) :: date
    character(len=*),          intent(in)                :: geo_finegrid_flnm
    real, allocatable, DIMENSION(:,:)                   :: xdumd  !-- decimated variable
    real*8, allocatable, DIMENSION(:)                   :: xcoord_d, xcoord
    real*8, allocatable, DIMENSION(:)                   :: ycoord_d, ycoord

    integer, save :: ncid,ncstatic
    integer, save :: output_count

    integer :: dimid_ix, dimid_jx, dimid_times, dimid_datelen, varid, n
    integer :: iret, dimid_soil, i,j,ii,jj, iret_x, iret_y
    character(len=256) :: output_flnm
    character(len=19)  :: date19
    character(len=32)  :: convention
    character(len=34)  :: sec_since_date

    character(len=30)  :: soilm

    real                                :: long_cm,lat_po,fe,fn, chan_in
    real, dimension(2)                  :: sp

    real, OPTIONAL,dimension(ixrt,jxrt) :: head, convgw, excess, &
                                  latval, lonval

    integer :: seconds_since, decimation, ixrtd,jxrtd, hires_flag
    
#ifdef MPP_LAND
    real, dimension(global_rt_nx,global_rt_ny) :: gHead, gConvgw, gExcess                                                  
    real, dimension(global_rt_nx,global_rt_ny) :: gLatval, gLonval
#endif
    
#ifdef MPP_LAND
    call MPP_LAND_COM_REAL(convgw, ixrt, jxrt, 99)
    call write_IO_rt_real(latval,gLatval)
    call write_IO_rt_real(lonval,gLonval)
    call write_IO_rt_real(head,gHead)
    call write_IO_rt_real(convgw,gConvgw)
    call write_IO_rt_real(excess,gExcess)
    

   if(my_id.eq.IO_id) then
     

#endif
    sec_since_date = 'seconds since '//date(1:4)//'-'//date(6:7)//'-'//date(9:10)//' '//date(12:13)//':'//date(15:16)//' UTC'
    seconds_since = int(dt)*output_count

    decimation = 1 !-- decimation factor
#ifdef MPP_LAND
    ixrtd = int(global_rt_nx/decimation)
    jxrtd = int(global_rt_ny/decimation)
#else
    ixrtd = int(ixrt/decimation)
    jxrtd = int(jxrt/decimation)
#endif
    allocate(xdumd(ixrtd,jxrtd))
    allocate(xcoord_d(ixrtd))
    allocate(ycoord_d(jxrtd))
    allocate(xcoord(ixrtd))
    allocate(ycoord(jxrtd))
    ii = 0
    jj = 0

    if (output_gw.eq.0) return  ! return if hires flag eq 0, if =1 output full grid

    if (output_count == 0) then

   !-- Open the  finemesh static files to obtain projection information
#ifdef HYDRO_D
      write(*,'("geo_finegrid_flnm: ''", A, "''")') trim(geo_finegrid_flnm)

#endif
      iret = nf_open(trim(geo_finegrid_flnm), NF_NOWRITE, ncstatic)

      if (iret /= 0) then
#ifdef HYDRO_D
         write(*,'("Problem opening geo_finegrid file: ''", A, "''")') &
         trim(geo_finegrid_flnm)
         write(*,*) "HIRES_OUTPUT will not be georeferenced..."
#endif
        hires_flag = 0
      else
        hires_flag = 1
      endif

     if(hires_flag.eq.1) then !if/then hires_georef
      ! Get Latitude (X)
      iret_x = NF_INQ_VARID(ncstatic,'x',varid)
      if(iret_x .eq. 0) iret = NF_GET_VAR_DOUBLE(ncstatic, varid, xcoord)
      ! Get Longitude (Y)
      iret_y = NF_INQ_VARID(ncstatic,'y',varid)
      if(iret_y .eq. 0) iret = NF_GET_VAR_DOUBLE(ncstatic, varid, ycoord)
     else
      xcoord_d = 0.
      ycoord_d = 0.
     end if  !endif hires_georef 

     do j=jxrtd,1,-1*decimation
        jj = jj+1
        if (jj<= jxrtd) then
         ycoord_d(jj) = ycoord(j)
        endif
     enddo

!yw     do i = 1,ixrt,decimation
!yw        ii = ii + 1
!yw        if (ii <= ixrtd) then 
!yw         xcoord_d(ii) = xcoord(i)
         xcoord_d = xcoord
!yw        endif
!yw     enddo
       

     if(hires_flag.eq.1) then !if/then hires_georef
      ! Get projection information from finegrid netcdf file
      iret = NF_INQ_VARID(ncstatic,'lambert_conformal_conic',varid)
      if(iret .eq. 0) iret = NF_GET_ATT_REAL(ncstatic, varid, 'longitude_of_central_meridian', long_cm)  !-- read it from the static file
      iret = NF_GET_ATT_REAL(ncstatic, varid, 'latitude_of_projection_origin', lat_po)  !-- read it from the static file
      iret = NF_GET_ATT_REAL(ncstatic, varid, 'false_easting', fe)  !-- read it from the static file
      iret = NF_GET_ATT_REAL(ncstatic, varid, 'false_northing', fn)  !-- read it from the static file
      iret = NF_GET_ATT_REAL(ncstatic, varid, 'standard_parallel', sp)  !-- read it from the static file
     end if  !endif hires_georef 
      iret = nf_close(ncstatic)

!-- create the fine grid routing file
     ! outname option is added by Qian, 2016.11
     if (outname == "") then  
        write(output_flnm, '(A12,".GW_SPINUP",I1)') date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16), igrid
     else 
        write(output_flnm, '(A, "_SPINUP_", A12)') trim(outname), date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16)
     endif

#ifdef HYDRO_D
       print*, 'output_flnm = "'//trim(output_flnm)//'"'
#endif


#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
       iret = nf_create(trim(output_flnm), IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
#else
       iret = nf_create(trim(output_flnm), NF_CLOBBER, ncid)
#endif

#ifdef HYDRO_D
       if (iret /= 0) then
         print*, "Problem nf_create"
         call hydro_stop("output_rt")
       endif
#endif

       iret = nf_def_dim(ncid, "time", NF_UNLIMITED, dimid_times)
       iret = nf_def_dim(ncid, "x", ixrtd, dimid_ix)  !-- make a decimated grid
       iret = nf_def_dim(ncid, "y", jxrtd, dimid_jx)

!--- define variables
       !- time definition, timeObs
       iret = nf_def_var(ncid,"time",NF_INT, 1, (/dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',34,sec_since_date)

       !- x-coordinate in cartesian system
       if (iret_x .eq. 0) then
       iret = nf_def_var(ncid,"x",NF_DOUBLE, 1, (/dimid_ix/), varid)
       iret = nf_put_att_text(ncid,varid,'long_name',26,'x coordinate of projection')
       iret = nf_put_att_text(ncid,varid,'standard_name',23,'projection_x_coordinate')
       iret = nf_put_att_text(ncid,varid,'units',5,'Meter')
       endif     
 
       !- y-coordinate in cartesian ssystem
       if (iret_y .eq. 0) then
       iret = nf_def_var(ncid,"y",NF_DOUBLE, 1, (/dimid_jx/), varid)
       iret = nf_put_att_text(ncid,varid,'long_name',26,'y coordinate of projection')
       iret = nf_put_att_text(ncid,varid,'standard_name',23,'projection_y_coordinate')
       iret = nf_put_att_text(ncid,varid,'units',5,'Meter')
       endif

       !- LATITUDE
       if(present(LATVAL)) then
       iret = nf_def_var(ncid,"LATITUDE",NF_FLOAT, 2, (/dimid_ix,dimid_jx/), varid)
       iret = nf_put_att_text(ncid,varid,'long_name',8,'LATITUDE')
       iret = nf_put_att_text(ncid,varid,'standard_name',8,'LATITUDE')
       iret = nf_put_att_text(ncid,varid,'units',5,'deg North')
       endif

       !- LONGITUDE
       if(present(LONVAL)) then
       iret = nf_def_var(ncid,"LONGITUDE",NF_FLOAT, 2, (/dimid_ix,dimid_jx/), varid)
       iret = nf_put_att_text(ncid,varid,'long_name',9,'LONGITUDE')
       iret = nf_put_att_text(ncid,varid,'standard_name',9,'LONGITUDE')
       iret = nf_put_att_text(ncid,varid,'units',5,'deg east')
       endif


       iret = nf_def_var(ncid,"GwHead",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',17,'groundwater head')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

       iret = nf_def_var(ncid,"GwConv",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'mm')
       iret = nf_put_att_text(ncid,varid,'long_name',12,'groundwater convergence')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       
       iret = nf_def_var(ncid,"GwExcess",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',17,'surface excess groundwater')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)

!-- place projection information
     if(hires_flag.eq.1) then !if/then hires_georef
      iret = nf_def_var(ncid,"lambert_conformal_conic",NF_INT,0, 0,varid)
      iret = nf_put_att_text(ncid,varid,'grid_mapping_name',23,'lambert_conformal_conic')
      iret = nf_put_att_real(ncid,varid,'longitude_of_central_meridian',NF_FLOAT,1,long_cm)
      iret = nf_put_att_real(ncid,varid,'latitude_of_projection_origin',NF_FLOAT,1,lat_po)
      iret = nf_put_att_real(ncid,varid,'false_easting',NF_FLOAT,1,fe)
      iret = nf_put_att_real(ncid,varid,'false_northing',NF_FLOAT,1,fn)
      iret = nf_put_att_real(ncid,varid,'standard_parallel',NF_FLOAT,2,sp)
     end if   !endif hires_georef

!      iret = nf_def_var(ncid,"Date",   NF_CHAR,  2, (/dimid_datelen,dimid_times/),     varid)

      date19(1:19) = "0000-00-00_00:00:00"
      date19(1:len_trim(startdate)) = startdate
      convention(1:32) = "CF-1.0"
      iret = nf_put_att_real(ncid, NF_GLOBAL, "missing_value", NF_FLOAT, 1, -9E15)
      iret = nf_put_att_text(ncid, NF_GLOBAL, "Conventions",6, convention)
      iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_start",19, date19)
      iret = nf_put_att_int(ncid,NF_GLOBAL,"output_decimation_factor",NF_INT, 1,decimation)

      iret = nf_enddef(ncid)

!!-- write latitude and longitude locations
        if (iret_x .eq. 0) then
!       xdumd = LATVAL
        iret = nf_inq_varid(ncid,"x", varid)
!       iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ixrtd,jxrtd/), xdumd)
	iret = nf_put_vara_double(ncid, varid, (/1/), (/ixrtd/), xcoord_d) !-- 1-d array
        endif

        if (iret_y .eq. 0) then
!       xdumd = LONVAL
        iret = nf_inq_varid(ncid,"y", varid)
!       iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ixrtd,jxrtd/), xdumd)
        iret = nf_put_vara_double(ncid, varid, (/1/), (/jxrtd/), ycoord_d) !-- 1-d array
        endif

        if(present(LATVAL)) then   !Qian, 2016.09
#ifdef MPP_LAND
        xdumd = gLATVAL
#else  
        xdumd = LATVAL
#endif
        iret = nf_inq_varid(ncid,"LATITUDE", varid)
        iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ixrtd,jxrtd/), xdumd)
        endif

        if(present(LONVAL)) then   !Qian, 2016.09
#ifdef MPP_LAND
        xdumd = gLONVAL
#else  
        xdumd = LONVAL
#endif
        iret = nf_inq_varid(ncid,"LONGITUDE", varid)
        iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ixrtd,jxrtd/), xdumd)
        endif

    endif

    output_count = output_count + 1

!!-- time
        iret = nf_inq_varid(ncid,"time", varid)
        iret = nf_put_vara_int(ncid, varid, (/output_count/), (/1/), seconds_since)


#ifdef MPP_LAND
        xdumd = gHead
#else  
        xdumd = head
#endif

     iret = nf_inq_varid(ncid,  "GwHead", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)

#ifdef MPP_LAND
        xdumd = gConvgw
#else  
        xdumd = convgw
#endif
     iret = nf_inq_varid(ncid,  "GwConv", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)

     
#ifdef MPP_LAND
        xdumd = gExcess
#else  
        xdumd = excess
#endif
     iret = nf_inq_varid(ncid,  "GwExcess", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    

!!time in seconds since startdate

       iret = nf_redef(ncid)
       date19(1:len_trim(date)) = date
       iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_end", 19, date19)
 
       iret = nf_enddef(ncid)
      iret = nf_sync(ncid)
      if (output_count == split_output_count) then
        output_count = 0
        iret = nf_close(ncid)
      endif

     deallocate(xdumd)
     deallocate(xcoord_d)
     deallocate(xcoord)
     deallocate(ycoord_d)
     deallocate(ycoord)
 
#ifdef HYDRO_D 
    write(6,*) "finish output for gw2d-spinup" 
#endif
   
#ifdef MPP_LAND
    endif
#endif

  end subroutine output_gw_spinup


  ! modified by Qian, 2016.09
  subroutine sub_output_gw(igrid,      split_output_count, ixrt,   jxrt,   &
                           nsoil,      startdate,          date,   HEAD,   &
                           SMCRT,      convgw,             excess, qsgwrt, &
                           qgw_chanrt, geo_finegrid_flnm,  dt,     sldpth, &
                           LATVAL,     LONVAL,             dist,           &
                           output_gw,  outname )

#ifdef MPP_LAND
    USE module_mpp_land
#endif
    ! output the routing variables over routing grid.
    implicit none
#include <netcdf.inc>

    integer,                            intent(in) :: igrid
    integer,                            intent(in) :: split_output_count
    integer,                            intent(in) :: ixrt,jxrt
    real,                               intent(in) :: dt
    real,                               intent(in) :: dist(ixrt,jxrt,9)
    integer,                            intent(in) :: nsoil
    integer,                            intent(in) :: output_gw
    character(len=*),                   intent(in) :: startdate
    character(len=*),                   intent(in) :: date
    character(len=*),                   intent(in) :: geo_finegrid_flnm
    character(len=*),                   intent(in) :: outname
    real,   dimension(nsoil),           intent(in) :: sldpth
    real,   OPTIONAL, dimension(ixrt,jxrt),       intent(in) :: head
    real,   OPTIONAL, dimension(ixrt,jxrt,nsoil), intent(in) :: SMCRT
    real,   OPTIONAL, dimension(ixrt,jxrt),       intent(in) :: convgw
    real,   OPTIONAL, dimension(ixrt,jxrt),       intent(in) :: excess
    real,   OPTIONAL, dimension(ixrt,jxrt),       intent(in) :: qsgwrt
    real,   OPTIONAL, dimension(ixrt,jxrt),       intent(in) :: qgw_chanrt
    real,   OPTIONAL, dimension(ixrt,jxrt),       intent(in) :: latval  !(lat,lon)
    real,   OPTIONAL, dimension(ixrt,jxrt),       intent(in) :: lonval  !(lat,lon)

    ! LOCAL VARIABLES
    real,   allocatable, DIMENSION(:,:)  :: xdumd  !-- decimated variable
    real*8, allocatable, DIMENSION(:)    :: xcoord_d, xcoord
    real*8, allocatable, DIMENSION(:)    :: ycoord_d, ycoord

    integer, save :: ncid,ncstatic
    integer, save :: output_count
    real,    dimension(nsoil) :: asldpth
    real,    dimension(2)     :: sp

    integer :: dimid_ix, dimid_jx, dimid_times, dimid_datelen, varid, n
    integer :: iret, dimid_soil, i,j,ii,jj
    real    :: long_cm,lat_po,fe,fn, chan_in
    integer :: seconds_since, decimation, ixrtd,jxrtd, hires_flag
    integer :: iret_proj, iret_x, iret_y    
    character(len=256) :: output_flnm
    character(len=19)  :: date19
    character(len=32)  :: convention
    character(len=34)  :: sec_since_date
!    character(len=30)  :: soilm
    
#ifdef MPP_LAND
    real, allocatable :: gHead(:,:)
    real, allocatable :: gConvgw(:,:)
    real, allocatable :: gqsgwrt(:,:)
    real, allocatable :: gExcess(:,:)
    real, allocatable :: gQgw_chanrt(:,:)
    real, allocatable :: gLatval(:,:)
    real, allocatable :: gLonval(:,:)
    real, allocatable :: gSMCRT(:,:,:)

    if (my_id .eq. IO_id) then
       if (present(head)) allocate(gHead(global_rt_nx,global_rt_ny))
       if (present(convgw)) allocate(gConvgw(global_rt_nx,global_rt_ny))
       if (present(qsgwrt)) allocate(gqsgwrt(global_rt_nx,global_rt_ny))
       if (present(excess)) allocate(gExcess(global_rt_nx,global_rt_ny))
       if (present(qgw_chanrt)) allocate(gQgw_chanrt(global_rt_nx,global_rt_ny))
       if (present(latval)) allocate(gLatval(global_rt_nx,global_rt_ny))
       if (present(lonval)) allocate(gLonval(global_rt_nx,global_rt_ny))
       if (present(SMCRT)) allocate(gSMCRT(global_rt_nx,global_rt_ny,nsoil))
   else
       if (present(head)) allocate(gHead(1,1))
       if (present(convgw)) allocate(gConvgw(1,1))
       if (present(qsgwrt)) allocate(gqsgwrt(1,1))
       if (present(excess)) allocate(gExcess(1,1))
       if (present(qgw_chanrt)) allocate(gQgw_chanrt(1,1))
       if (present(latval)) allocate(gLatval(1,1))
       if (present(lonval)) allocate(gLonval(1,1))
       if (present(SMCRT)) allocate(gSMCRT(1,1,1))
   endif
#endif
    
#ifdef MPP_LAND
    if (present(convgw)) then
       call MPP_LAND_COM_REAL(convgw, ixrt, jxrt, 99)
       call write_IO_rt_real(convgw,gConvgw)
    endif

    if (present(qsgwrt)) then
       call MPP_LAND_COM_REAL(qsgwrt, ixrt, jxrt, 99)
       call write_IO_rt_real(qsgwrt,gqsgwrt)
    endif

    if (present(qgw_chanrt)) then
       call MPP_LAND_COM_REAL(qgw_chanrt, ixrt, jxrt, 99)
       call write_IO_rt_real(qgw_chanrt,gQgw_chanrt)
    endif

    if (present(latval)) &
       call write_IO_rt_real(latval,gLatval)

    if (present(lonval)) &
       call write_IO_rt_real(lonval,gLonval)

    if (present(head)) then
       call MPP_LAND_COM_REAL(head, ixrt, jxrt, 99)
       call write_IO_rt_real(head,gHead)
    endif
 
    if (present(excess)) then
       call MPP_LAND_COM_REAL(excess, ixrt, jxrt, 99)
       call write_IO_rt_real(excess,gExcess)
    endif    

    if (present(smcrt)) then
    do i = 1, NSOIL
       call MPP_LAND_COM_REAL(smcrt(:,:,i), ixrt, jxrt, 99)
       call write_IO_rt_real(SMCRT(:,:,i),gSMCRT(:,:,i))
    end do
    endif

   if(my_id.eq.IO_id) then

#endif
    sec_since_date = 'seconds since '//date(1:4)//'-'//date(6:7)//'-'//date(9:10)//' '//date(12:13)//':'//date(15:16)//' UTC'
    seconds_since = int(dt)*output_count

    decimation = 1 !-- decimation factor
#ifdef MPP_LAND
    ixrtd = int(global_rt_nx/decimation)
    jxrtd = int(global_rt_ny/decimation)
#else
    ixrtd = int(ixrt/decimation)
    jxrtd = int(jxrt/decimation)
#endif
    allocate(xdumd(ixrtd,jxrtd))
    allocate(xcoord_d(ixrtd))
    allocate(ycoord_d(jxrtd))
    allocate(xcoord(ixrtd))
    allocate(ycoord(jxrtd))
    ii = 0
    jj = 0

    if (output_gw .eq. 0) return  ! return if hires flag eq 0, if =1 output full grid

    if (output_count == 0) then
     !-- Open the  finemesh static files to obtain projection information
#ifdef HYDRO_D
     write(*,'("geo_finegrid_flnm: ''", A, "''")') trim(geo_finegrid_flnm)
#endif
     iret = nf_open(trim(geo_finegrid_flnm), NF_NOWRITE, ncstatic)

     if (iret /= 0) then
#ifdef HYDRO_D
       write(*,'("Problem opening geo_finegrid file: ''", A, "''")') &
       trim(geo_finegrid_flnm)
       write(*,*) "HIRES_OUTPUT will not be georeferenced..."
#endif
       hires_flag = 0
     else
       hires_flag = 1
     endif

     ! obtain projection coordinates
     if(hires_flag.eq.1) then !if/then hires_georef
       ! Get Latitude (X)
       iret_x = NF_INQ_VARID(ncstatic, 'x',varid)
       if(iret_x .eq. 0) iret = NF_GET_VAR_DOUBLE(ncstatic, varid, xcoord)
       ! Get Longitude (Y)
       iret_y = NF_INQ_VARID(ncstatic,'y',varid)
       if(iret_y .eq. 0) iret = NF_GET_VAR_DOUBLE(ncstatic, varid, ycoord)
     else
       xcoord_d = 0.
       ycoord_d = 0.
     end if  !endif hires_georef 

     do j=jxrtd,1,-1*decimation
        jj = jj+1
        if (jj<= jxrtd) then
         ycoord_d(jj) = ycoord(j)
        endif
     enddo

!yw     do i = 1,ixrt,decimation
!yw        ii = ii + 1
!yw        if (ii <= ixrtd) then 
!yw         xcoord_d(ii) = xcoord(i)
         xcoord_d = xcoord
!yw        endif
!yw     enddo
       
     if(hires_flag .eq. 1) then !if/then hires_georef
       ! Get projection information from finegrid netcdf file
       iret_proj = NF_INQ_VARID(ncstatic,'lambert_conformal_conic',varid)
       if(iret_proj .eq. 0) then
         iret = NF_GET_ATT_REAL(ncstatic, varid, 'longitude_of_central_meridian', long_cm) !-- read it from the static file
         iret = NF_GET_ATT_REAL(ncstatic, varid, 'latitude_of_projection_origin', lat_po)  !-- read it from the static file
         iret = NF_GET_ATT_REAL(ncstatic, varid, 'false_easting', fe)                      !-- read it from the static file
         iret = NF_GET_ATT_REAL(ncstatic, varid, 'false_northing', fn)                     !-- read it from the static file
         iret = NF_GET_ATT_REAL(ncstatic, varid, 'standard_parallel', sp)                  !-- read it from the static file
       endif
     end if  !endif hires_georef 
     iret = nf_close(ncstatic)

     !-- create the fine grid routing file
     ! outname option is added by Qian, 2016.11
     if (outname == "") then  
        write(output_flnm, '(A12,".GW_DOMAIN",I1)') date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16), igrid
     else
        write(output_flnm, '(A, I1, A, A12)') trim(outname), igrid, "_",date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16)
     endif
#ifdef HYDRO_D
     print*, 'output_flnm = "'//trim(output_flnm)//'"'
#endif

#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
       iret = nf_create(trim(output_flnm), IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
#else
       iret = nf_create(trim(output_flnm), NF_CLOBBER, ncid)
#endif

#ifdef HYDRO_D
       if (iret /= 0) then
         print*, "Problem nf_create"
         call hydro_stop("output_rt")
       endif
#endif

       iret = nf_def_dim(ncid, "time",  NF_UNLIMITED, dimid_times)
       iret = nf_def_dim(ncid, "x",     ixrtd,        dimid_ix)    !-- make a decimated grid
       iret = nf_def_dim(ncid, "y",     jxrtd,        dimid_jx)
       iret = nf_def_dim(ncid, "depth", nsoil,        dimid_soil)  !-- 3-d soils

       !--- define variables
       !- time definition, timeObs
       iret = nf_def_var(ncid,"time",NF_INT, 1, (/dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',34,sec_since_date)

       !- x-coordinate in cartesian system
       if (iret_x .eq. 0) then    ! Qian, 2016.09
       iret = nf_def_var(ncid,"x",NF_DOUBLE, 1, (/dimid_ix/), varid)
       iret = nf_put_att_text(ncid,varid,'long_name',26,'x coordinate of projection')
       iret = nf_put_att_text(ncid,varid,'standard_name',23,'projection_x_coordinate')
       iret = nf_put_att_text(ncid,varid,'units',5,'Meter')
       endif

       !- y-coordinate in cartesian ssystem
       if(iret_y .eq. 0) then
       iret = nf_def_var(ncid,"y",NF_DOUBLE, 1, (/dimid_jx/), varid)
       iret = nf_put_att_text(ncid,varid,'long_name',26,'y coordinate of projection')
       iret = nf_put_att_text(ncid,varid,'standard_name',23,'projection_y_coordinate')
       iret = nf_put_att_text(ncid,varid,'units',5,'Meter')
       endif

       !- LATITUDE
       if(present(LATVAL)) then
       iret = nf_def_var(ncid,"LATITUDE",NF_FLOAT, 2, (/dimid_ix,dimid_jx/), varid)
       iret = nf_put_att_text(ncid,varid,'long_name',8,'LATITUDE')
       iret = nf_put_att_text(ncid,varid,'standard_name',8,'LATITUDE')
       iret = nf_put_att_text(ncid,varid,'units',5,'deg North')
       endif

       !- LONGITUDE
       if(present(LONVAL)) then
       iret = nf_def_var(ncid,"LONGITUDE",NF_FLOAT, 2, (/dimid_ix,dimid_jx/), varid)
       iret = nf_put_att_text(ncid,varid,'long_name',9,'LONGITUDE')
       iret = nf_put_att_text(ncid,varid,'standard_name',9,'LONGITUDE')
       iret = nf_put_att_text(ncid,varid,'units',5,'deg east')
       endif

       !-- z-level is soil
       iret = nf_def_var(ncid,"depth", NF_FLOAT, 1, (/dimid_soil/),varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'cm')
       iret = nf_put_att_text(ncid,varid,'long_name',19,'depth of soil layer')

       if (present(SMCRT)) then
       iret = nf_def_var(ncid,"SOIL_M",  NF_FLOAT, 4, (/dimid_ix,dimid_jx,dimid_soil,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',6,'kg m-2')
!       iret = nf_put_att_text(ncid,varid,'description',16,'moisture content')
       iret = nf_put_att_text(ncid,varid,'description',16,'soil water volumetric content')
!       iret = nf_put_att_text(ncid,varid,'long_name',26,soilm)
       iret = nf_put_att_text(ncid,varid,'long_name',26,'soil moisture') 
!       iret = nf_put_att_text(ncid,varid,'coordinates',5,'x y z')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(HEAD)) then
       iret = nf_def_var(ncid,"HEAD",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',17,'groundwater head')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(convgw)) then
       iret = nf_def_var(ncid,"CONVGW",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'mm')
       iret = nf_put_att_text(ncid,varid,'long_name',12,'channel flux')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif
 
       if(present(excess)) then
       iret = nf_def_var(ncid,"GwExcess",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'mm')
       iret = nf_put_att_text(ncid,varid,'long_name',17,'surface excess groundwater')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(qsgwrt)) then
       iret = nf_def_var(ncid,"QSGWRT",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'mm')
       iret = nf_put_att_text(ncid,varid,'long_name',12,'surface head')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if (present(qgw_chanrt)) then
       iret = nf_def_var(ncid,"QGW_CHANRT",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'m3 s-1')
       iret = nf_put_att_text(ncid,varid,'long_name',12,'surface head')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       !-- place projection information
       if(hires_flag .eq. 1) then !if/then hires_georef
        if(iret_proj .eq. 0) then
        iret = nf_def_var(ncid,"lambert_conformal_conic",NF_INT,0, 0,varid)
        iret = nf_put_att_text(ncid,varid,'grid_mapping_name',23,'lambert_conformal_conic')
        iret = nf_put_att_real(ncid,varid,'longitude_of_central_meridian',NF_FLOAT,1,long_cm)
        iret = nf_put_att_real(ncid,varid,'latitude_of_projection_origin',NF_FLOAT,1,lat_po)
        iret = nf_put_att_real(ncid,varid,'false_easting',NF_FLOAT,1,fe)
        iret = nf_put_att_real(ncid,varid,'false_northing',NF_FLOAT,1,fn)
        iret = nf_put_att_real(ncid,varid,'standard_parallel',NF_FLOAT,2,sp)
        endif
       end if   !endif hires_georef

!      iret = nf_def_var(ncid,"Date",   NF_CHAR,  2, (/dimid_datelen,dimid_times/),     varid)

       date19(1:19) = "0000-00-00_00:00:00"
       date19(1:len_trim(startdate)) = startdate
       convention(1:32) = "CF-1.0"
       iret = nf_put_att_real(ncid, NF_GLOBAL, "missing_value", NF_FLOAT, 1, -9E15)
       iret = nf_put_att_text(ncid, NF_GLOBAL, "Conventions",6, convention)
       iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_start",19, date19)
       iret = nf_put_att_int(ncid,NF_GLOBAL,   "output_decimation_factor",NF_INT, 1,decimation)

       iret = nf_enddef(ncid)

        !-- write projection locations
        if (iret_x .eq. 0) then    !Qian, 2016.09
!       xdumd = LATVAL
        iret = nf_inq_varid(ncid,"x", varid)
!       iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ixrtd,jxrtd/), xdumd)
	iret = nf_put_vara_double(ncid, varid, (/1/), (/ixrtd/), xcoord_d) !-- 1-d array
        endif

        if (iret_y .eq. 0) then    !Qian, 2016.09
!       xdumd = LONVAL
        iret = nf_inq_varid(ncid,"y", varid)
!       iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ixrtd,jxrtd/), xdumd)
        iret = nf_put_vara_double(ncid, varid, (/1/), (/jxrtd/), ycoord_d) !-- 1-d array
        endif

        !-- write latitude and longitude locations
        if(present(LATVAL)) then   !Qian, 2016.09
#ifdef MPP_LAND
        xdumd = gLATVAL
#else  
        xdumd = LATVAL
#endif
        iret = nf_inq_varid(ncid,"LATITUDE", varid)
        iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ixrtd,jxrtd/), xdumd)
        endif

        if(present(LONVAL)) then   !Qian, 2016.09
#ifdef MPP_LAND
        xdumd = gLONVAL
#else  
        xdumd = LONVAL
#endif
        iret = nf_inq_varid(ncid,"LONGITUDE", varid)
        iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ixrtd,jxrtd/), xdumd)
        endif

        do n = 1,nsoil
           if(n == 1) then
             asldpth(n) = -sldpth(n)
           else
             asldpth(n) = asldpth(n-1) - sldpth(n)
           endif
        end do

       iret = nf_inq_varid(ncid,"depth", varid)
       iret = nf_put_vara_real(ncid, varid, (/1/), (/nsoil/), asldpth)

    endif

    output_count = output_count + 1

    !-- time
    iret = nf_inq_varid(ncid,"time", varid)
    iret = nf_put_vara_int(ncid, varid, (/output_count/), (/1/), seconds_since)

    !-- 3-d soils
    if (present(SMCRT)) then
      do n = 1, nsoil
#ifdef MPP_LAND
       xdumd = gSMCRT(:,:,n)
#else  
       xdumd = SMCRT(:,:,n)
#endif
! !DJG inv      jj = int(jxrt/decimation)
!       jj = 1
!       ii = 0
! !DJG inv      do j = jxrt,1,-decimation
!        do j = 1,jxrt,decimation
!        do i = 1,ixrt,decimation
!         ii = ii + 1  
!         if(ii <= ixrtd .and. jj <= jxrtd .and. jj >0) then
!          xdumd(ii,jj) = smcrt(i,j,n)
!         endif
!       enddo 
!        ii = 0
! !DJG inv       jj = jj -1
!        jj = jj + 1
!      enddo
!       where (vegtyp(:,:) == 16) xdum = -1.E33
          iret = nf_inq_varid(ncid,  "SOIL_M", varid)
          iret = nf_put_vara_real(ncid, varid, (/1,1,n,output_count/), (/ixrtd,jxrtd,1,1/), xdumd)
      enddo !-n soils
    endif

    if (present(HEAD)) then
#ifdef MPP_LAND
     xdumd = gHead
#else  
     xdumd = head
#endif
     iret = nf_inq_varid(ncid,  "HEAD", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(convgw)) then
#ifdef MPP_LAND
     xdumd = gConvgw
#else  
     xdumd = convgw
#endif
     iret = nf_inq_varid(ncid,  "CONVGW", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(excess)) then     
#ifdef MPP_LAND
     xdumd = gExcess
#else  
     xdumd = excess
#endif
     iret = nf_inq_varid(ncid,  "GwExcess", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(qsgwrt)) then
#ifdef MPP_LAND
     xdumd = gqsgwrt
#else  
     xdumd = qsgwrt
#endif
     iret = nf_inq_varid(ncid,  "QSGWRT", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(qgw_chanrt)) then 
#ifdef MPP_LAND
     xdumd = gQgw_chanrt
#else  
     xdumd = qgw_chanrt
#endif
     iret = nf_inq_varid(ncid,  "QGW_CHANRT", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif
     
    ! time in seconds since startdate
    iret = nf_redef(ncid)
    date19(1:len_trim(date)) = date
    iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_end", 19, date19)
 
    iret = nf_enddef(ncid)
    iret = nf_sync(ncid)
    if (output_count == split_output_count) then
        output_count = 0
        iret = nf_close(ncid)
    endif

    deallocate(xdumd)
    deallocate(xcoord_d)
    deallocate(xcoord)
    deallocate(ycoord_d)
    deallocate(ycoord)
    
#ifdef HYDRO_D 
    write(6,*) "finish output for gw2d" 
#endif

#ifdef MPP_LAND
    endif    ! if (my_id .eq. IO_id) then
#endif

  end subroutine sub_output_gw

  ! output overland-channel routing
  ! added by Qian, 2016.09 
  subroutine sub_output_rt2d(igrid,      split_output_count, ixrt,      jxrt, &
                             nsoil,      startdate,          date,        dt, &
                             output_rt2d,outname,                             &
                             gwhead,     runoff,            r_sat,     sinfil,&
                             tinfil,     fwater,             fsat,    cuiflux,&
                             cgiflux,    cglflux,           hlake, streamflux,&
                         istreamflux,    w_us,               w_gw,       w_sf,&
                             ac_pmewsf,  ac_infil_grnd, ac_drain_grnd,        &
                             ac_roff_grnd,ac_rsat_grnd, ac_precip,    ac_evap,&
                             ac_cgi,     ac_cgl,           ac_cui,   ac_ocngw,&
                             ac_ocnsf,   ac_fixgw,       ac_fixsf,      ltype,&
                             elev,       bank,                width,      bot,&
                             fbed,topo,  irg,               sabs, ac_precevap,&
                             ac_irgrw,   ac_irggw,            veloc,          &
                             qsKR,       NumFracts,         SusSed,    DepSed,&
                             ErodSed,    SedFlux,        t_SedFlux,      qsEH,&
                             SedConc,    T_SedConc)
#ifdef MPP_LAND
    USE module_mpp_land
#endif
                             
    implicit none
#include <netcdf.inc>

    integer,  intent(in)                       :: igrid
    integer,  intent(in)                       :: split_output_count
    integer,  intent(in)                       :: ixrt, jxrt
    integer,  intent(in)                       :: nsoil
    integer,  intent(in)                       :: NumFracts
    real,     intent(in)                       :: dt
    character(len=*), intent(in)               :: startdate
    character(len=*), intent(in)               :: date
    character(len=*), intent(in)               :: outname
    integer,  intent(in)                       :: output_rt2d
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: gwhead
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: runoff
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: r_sat
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: sinfil
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: tinfil
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: fwater
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: fsat
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: cuiflux
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: cgiflux
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: cglflux
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: hlake
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: streamflux
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: veloc
    !
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: qsKR
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: T_SedConc
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: t_SedFlux
    real,     intent(in), OPTIONAL, dimension(NumFracts,ixrt,jxrt) :: &
                                             qsEH,SusSed,DepSed,ErodSed,SedFlux,SedConc

    integer,  intent(in), OPTIONAL, dimension(ixrt,jxrt) :: istreamflux 
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: w_us
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: w_gw
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: w_sf
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_pmewsf
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_infil_grnd
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_drain_grnd
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_roff_grnd
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_rsat_grnd
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_precip
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_evap
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_cgi
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_cgl
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_cui
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_precevap
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_irgrw
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_irggw
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_ocngw
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_ocnsf
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_fixgw
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ac_fixsf
    integer,  intent(in), OPTIONAL, dimension(ixrt,jxrt) :: ltype
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: elev
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: bank
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: width
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: bot
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: fbed
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: topo
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: irg
    real,     intent(in), OPTIONAL, dimension(ixrt,jxrt) :: sabs

    ! LOCAL VARIABLES
    real,   allocatable, DIMENSION(:,:)  :: xdumd  !-- decimated variable
    integer,allocatable, DIMENSION(:,:)  :: idumd 

    integer, save :: ncid,ncstatic
    integer, save :: output_count
    real,    dimension(nsoil) :: asldpth
    real,    dimension(2)     :: sp

    integer :: dimid_ix, dimid_jx, dimid_times, dimid_datelen, varid, n
    integer :: iret, dimid_soil, i,j,ii,jj, dimid_sFracts
    real    :: long_cm,lat_po,fe,fn, chan_in
    integer :: seconds_since, decimation, ixrtd,jxrtd, hires_flag
    integer :: iret_proj, iret_x, iret_y
    character(len=256) :: output_flnm
    character(len=19)  :: date19
    character(len=32)  :: convention
    character(len=34)  :: sec_since_date

#ifdef MPP_LAND
    real, allocatable :: gGwhead(:,:)
    real, allocatable :: gRunoff(:,:)
    real, allocatable :: gR_sat(:,:)
    real, allocatable :: gSinfil(:,:)
    real, allocatable :: gTinfil(:,:)
    real, allocatable :: gFwater(:,:)
    real, allocatable :: gFsat(:,:) 
    real, allocatable :: gCuiflux(:,:)
    real, allocatable :: gCgiflux(:,:)
    real, allocatable :: gCglflux(:,:)
    real, allocatable :: gHlake(:,:)
    real, allocatable :: gStreamflux(:,:)
    real, allocatable :: gVeloc(:,:)
    integer, allocatable :: gIstreamflux(:,:)
    real, allocatable :: gW_us(:,:)
    real, allocatable :: gW_gw(:,:)
    real, allocatable :: gW_sf(:,:)
    real, allocatable :: gAc_pmewsf(:,:)
    real, allocatable :: gAc_infil_grnd(:,:)
    real, allocatable :: gAc_drain_grnd(:,:)
    real, allocatable :: gAc_roff_grnd(:,:)
    real, allocatable :: gAc_rsat_grnd(:,:)
    real, allocatable :: gAc_precip(:,:)
    real, allocatable :: gAc_evap(:,:)
    real, allocatable :: gAc_cgi(:,:)
    real, allocatable :: gAc_cgl(:,:)
    real, allocatable :: gAc_cui(:,:)
    real, allocatable :: gAc_precevap(:,:)
    real, allocatable :: gAc_ocngw(:,:)
    real, allocatable :: gAc_ocnsf(:,:)
    real, allocatable :: gAc_fixgw(:,:)
    real, allocatable :: gAc_fixsf(:,:)
    integer, allocatable :: gLtype(:,:)
    real, allocatable :: gElev(:,:)
    real, allocatable :: gBank(:,:)
    real, allocatable :: gWidth(:,:)
    real, allocatable :: gBot(:,:)
    real, allocatable :: gFbed(:,:)
    real, allocatable :: gtopo(:,:)
    real, allocatable :: girg(:,:)
#endif

#ifdef MPP_LAND
    if (my_id .eq. IO_id) then
       if (present(gwhead)) allocate(gGwhead(global_rt_nx,global_rt_ny))
       if (present(runoff)) allocate(gRunoff(global_rt_nx,global_rt_ny))
       if (present(r_sat))  allocate(gR_sat(global_rt_nx,global_rt_ny))
       if (present(sinfil)) allocate(gSinfil(global_rt_nx,global_rt_ny)) 
       if (present(tinfil)) allocate(gTinfil(global_rt_nx,global_rt_ny))
       if (present(fwater)) allocate(gFwater(global_rt_nx,global_rt_ny))
       if (present(fsat))   allocate(gFsat(global_rt_nx,global_rt_ny))
       if (present(cuiflux)) allocate(gCuiflux(global_rt_nx,global_rt_ny))
       if (present(cgiflux)) allocate(gCgiflux(global_rt_nx,global_rt_ny))
       if (present(cglflux)) allocate(gCglflux(global_rt_nx,global_rt_ny))       
       if (present(hlake))   allocate(gHlake(global_rt_nx,global_rt_ny))
       if (present(streamflux))  allocate(gStreamflux(global_rt_nx,global_rt_ny))
       if (present(veloc))  allocate(gVeloc(global_rt_nx,global_rt_ny))
       if (present(istreamflux)) allocate(gIstreamflux(global_rt_nx,global_rt_ny))
       if (present(w_us)) allocate(gW_us(global_rt_nx,global_rt_ny))
       if (present(w_gw)) allocate(gW_gw(global_rt_nx,global_rt_ny))
       if (present(w_sf)) allocate(gW_sf(global_rt_nx,global_rt_ny))
       if (present(ac_pmewsf))     allocate(gAc_pmewsf(global_rt_nx,global_rt_ny))
       if (present(ac_infil_grnd)) allocate(gAc_infil_grnd(global_rt_nx,global_rt_ny))
       if (present(ac_drain_grnd)) allocate(gAc_drain_grnd(global_rt_nx,global_rt_ny))
       if (present(ac_roff_grnd))  allocate(gAc_roff_grnd(global_rt_nx,global_rt_ny))
       if (present(ac_rsat_grnd))  allocate(gAc_rsat_grnd(global_rt_nx,global_rt_ny))
       if (present(ac_precip)) allocate(gAc_precip(global_rt_nx,global_rt_ny))
       if (present(ac_evap))   allocate(gAc_evap(global_rt_nx,global_rt_ny))
       if (present(ac_cgi))    allocate(gAc_cgi(global_rt_nx,global_rt_ny))
       if (present(ac_cgl))    allocate(gAc_cgl(global_rt_nx,global_rt_ny)) 
       if (present(ac_cui))    allocate(gAc_cui(global_rt_nx,global_rt_ny))
       if (present(ac_precevap))    allocate(gAc_precevap(global_rt_nx,global_rt_ny))
       if (present(ac_ocngw))  allocate(gAc_ocngw(global_rt_nx,global_rt_ny))
       if (present(ac_ocnsf))  allocate(gAc_ocnsf(global_rt_nx,global_rt_ny))
       if (present(ac_fixgw))  allocate(gAc_fixgw(global_rt_nx,global_rt_ny))
       if (present(ac_fixsf))  allocate(gAc_fixsf(global_rt_nx,global_rt_ny))
       if (present(ltype))     allocate(gLtype(global_rt_nx,global_rt_ny))
       if (present(elev))      allocate(gElev(global_rt_nx,global_rt_ny))
       if (present(bank))      allocate(gBank(global_rt_nx,global_rt_ny))
       if (present(width))     allocate(gWidth(global_rt_nx,global_rt_ny))
       if (present(bot))       allocate(gBot(global_rt_nx,global_rt_ny))
       if (present(fbed))      allocate(gFbed(global_rt_nx,global_rt_ny))
       if (present(topo))      allocate(gtopo(global_rt_nx,global_rt_ny))
       if (present(irg))      allocate(girg(global_rt_nx,global_rt_ny))
    else
       if (present(gwhead)) allocate(gGwhead(1,1))
       if (present(runoff)) allocate(gRunoff(1,1))
       if (present(r_sat))  allocate(gR_sat(1,1))
       if (present(sinfil)) allocate(gSinfil(1,1))
       if (present(tinfil)) allocate(gTinfil(1,1))
       if (present(fwater)) allocate(gFwater(1,1))
       if (present(fsat))   allocate(gFsat(1,1))
       if (present(cuiflux)) allocate(gCuiflux(1,1))
       if (present(cgiflux)) allocate(gCgiflux(1,1))
       if (present(cglflux)) allocate(gCglflux(1,1))
       if (present(hlake))   allocate(gHlake(1,1))
       if (present(streamflux))  allocate(gStreamflux(1,1))
       if (present(veloc))  allocate(gVeloc(1,1))
       if (present(istreamflux)) allocate(gIstreamflux(1,1))
       if (present(w_us)) allocate(gW_us(1,1))
       if (present(w_gw)) allocate(gW_gw(1,1))
       if (present(w_sf)) allocate(gW_sf(1,1))
       if (present(ac_pmewsf))     allocate(gAc_pmewsf(1,1))
       if (present(ac_infil_grnd)) allocate(gAc_infil_grnd(1,1))
       if (present(ac_drain_grnd)) allocate(gAc_drain_grnd(1,1))
       if (present(ac_roff_grnd))  allocate(gAc_roff_grnd(1,1))
       if (present(ac_rsat_grnd))  allocate(gAc_rsat_grnd(1,1))
       if (present(ac_precip)) allocate(gAc_precip(1,1))
       if (present(ac_evap))   allocate(gAc_evap(1,1))
       if (present(ac_cgi))    allocate(gAc_cgi(1,1))
       if (present(ac_cgl))    allocate(gAc_cgl(1,1))
       if (present(ac_cui))    allocate(gAc_cui(1,1))
       if (present(ac_precevap))    allocate(gAc_precevap(1,1))
       if (present(ac_ocngw))  allocate(gAc_ocngw(1,1))
       if (present(ac_ocnsf))  allocate(gAc_ocnsf(1,1))
       if (present(ac_fixgw))  allocate(gAc_fixgw(1,1))
       if (present(ac_fixsf))  allocate(gAc_fixsf(1,1))
       if (present(ltype))     allocate(gLtype(1,1))
       if (present(elev))      allocate(gElev(1,1))
       if (present(bank))      allocate(gBank(1,1))
       if (present(width))     allocate(gWidth(1,1))
       if (present(bot))       allocate(gBot(1,1))
       if (present(fbed))      allocate(gFbed(1,1))
       if (present(topo))      allocate(topo(1,1))
       if (present(irg))      allocate(irg(1,1))
    endif

    if (present(gwhead)) then
       call MPP_LAND_COM_REAL(gwhead, ixrt, jxrt, 99)
       call write_IO_rt_real(gwhead,gGwhead)
    endif
    if (present(runoff)) then
       call MPP_LAND_COM_REAL(runoff, ixrt, jxrt, 99)
       call write_IO_rt_real(runoff,gRunoff)
    endif
    if (present(r_sat)) then
       call MPP_LAND_COM_REAL(r_sat, ixrt, jxrt, 99)
       call write_IO_rt_real(r_sat, gR_sat)
    endif
    if (present(sinfil)) then
       call MPP_LAND_COM_REAL(sinfil, ixrt, jxrt, 99)
       call write_IO_rt_real(sinfil,gSinfil)
    endif
    if (present(tinfil)) then
       call MPP_LAND_COM_REAL(tinfil, ixrt, jxrt, 99)
       call write_IO_rt_real(tinfil,gTinfil)
    endif
    if (present(fwater)) then
       call MPP_LAND_COM_REAL(fwater, ixrt, jxrt, 99)
       call write_IO_rt_real(fwater,gFwater)
    endif
    if (present(fsat)) then
       call MPP_LAND_COM_REAL(fsat, ixrt, jxrt, 99)
       call write_IO_rt_real(fsat,gFsat)
    endif
    if (present(cuiflux)) then
       call MPP_LAND_COM_REAL(cuiflux, ixrt, jxrt, 99)
       call write_IO_rt_real(cuiflux,gCuiflux)
    endif
    if (present(cgiflux)) then
       call MPP_LAND_COM_REAL(cgiflux, ixrt, jxrt, 99)
       call write_IO_rt_real(cgiflux,gCgiflux)
    endif
    if (present(cglflux)) then
       call MPP_LAND_COM_REAL(cglflux, ixrt, jxrt, 99)
       call write_IO_rt_real(cglflux,gCglflux)
    endif
    if (present(hlake)) then
       call MPP_LAND_COM_REAL(hlake, ixrt, jxrt, 99)
       call write_IO_rt_real(hlake,gHlake)
    endif
    if (present(streamflux)) then
       call MPP_LAND_COM_REAL(streamflux, ixrt, jxrt, 99)
       call write_IO_rt_real(streamflux,gStreamflux)
    endif
    if (present(veloc)) then
        call MPP_LAND_COM_REAL(veloc, ixrt, jxrt, 99)
        call write_IO_rt_real(veloc, gVeloc)
    endif
    if (present(istreamflux)) then
       call MPP_LAND_COM_INTEGER(istreamflux, ixrt, jxrt, 99)
       call write_IO_rt_int(istreamflux,gIstreamflux)
    endif
    if (present(w_us)) then
       call MPP_LAND_COM_REAL(w_us, ixrt, jxrt, 99)
       call write_IO_rt_real(w_us,gW_us)
    endif
    if (present(w_gw)) then
       call MPP_LAND_COM_REAL(w_gw, ixrt, jxrt, 99)
       call write_IO_rt_real(w_gw,gW_gw)
    endif
    if (present(w_sf)) then
       call MPP_LAND_COM_REAL(w_sf, ixrt, jxrt, 99)
       call write_IO_rt_real(w_sf,gW_sf)
    endif
    if (present(ac_pmewsf)) then
       call MPP_LAND_COM_REAL(ac_pmewsf, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_pmewsf,gAc_pmewsf)
    endif
    if (present(ac_infil_grnd)) then
       call MPP_LAND_COM_REAL(ac_infil_grnd, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_infil_grnd,gAc_infil_grnd)
    endif
    if (present(ac_drain_grnd)) then
       call MPP_LAND_COM_REAL(ac_drain_grnd, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_drain_grnd,gAc_drain_grnd)
    endif
    if (present(ac_roff_grnd)) then
       call MPP_LAND_COM_REAL(ac_roff_grnd, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_roff_grnd,gAc_roff_grnd)
    endif
    if (present(ac_rsat_grnd)) then
       call MPP_LAND_COM_REAL(ac_rsat_grnd, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_rsat_grnd,gAc_rsat_grnd)
    endif
    if (present(ac_precip)) then
       call MPP_LAND_COM_REAL(ac_precip, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_precip,gAc_precip)
    endif
    if (present(ac_evap)) then
       call MPP_LAND_COM_REAL(ac_evap, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_evap,gAc_evap)
    endif
    if (present(ac_cgi)) then
       call MPP_LAND_COM_REAL(ac_cgi, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_cgi,gAc_cgi)
    endif
    if (present(ac_cgl)) then
       call MPP_LAND_COM_REAL(ac_cgl, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_cgl,gAc_cgl)
    endif
    if (present(ac_cui)) then
       call MPP_LAND_COM_REAL(ac_cui, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_cui,gAc_cui)
    endif
    if (present(ac_precevap)) then
       call MPP_LAND_COM_REAL(ac_precevap, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_precevap,gAc_precevap)
    endif
    if (present(ac_ocngw)) then
       call MPP_LAND_COM_REAL(ac_ocngw, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_ocngw,gAc_ocngw)
    endif
    if (present(ac_ocnsf)) then
       call MPP_LAND_COM_REAL(ac_ocnsf, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_ocnsf,gAc_ocnsf)
    endif
    if (present(ac_fixgw)) then
       call MPP_LAND_COM_REAL(ac_fixgw, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_fixgw,gAc_fixgw)
    endif
    if (present(ac_fixsf)) then
       call MPP_LAND_COM_REAL(ac_fixsf, ixrt, jxrt, 99)
       call write_IO_rt_real(ac_fixsf,gAc_fixsf)
    endif
    if (present(ltype)) call write_IO_rt_int(ltype,gLtype)
    if (present(elev))  call write_IO_rt_real(elev,gElev)
    if (present(bank))  call write_IO_rt_real(bank,gBank)
    if (present(width)) call write_IO_rt_real(width,gWidth)
    if (present(bot))   call write_IO_rt_real(bot,gBot)
    if (present(fbed))  call write_IO_rt_real(fbed,gFbed)

  if(my_id .eq. IO_id) then

#endif

    sec_since_date = 'seconds since '//date(1:4)//'-'//date(6:7)//'-'//date(9:10)//' '//date(12:13)//':'//date(15:16)//' UTC'
    seconds_since = int(dt)*output_count

    decimation = 1 !-- decimation factor
#ifdef MPP_LAND
    ixrtd = int(global_rt_nx/decimation)
    jxrtd = int(global_rt_ny/decimation)
#else
    ixrtd = int(ixrt/decimation)
    jxrtd = int(jxrt/decimation)
#endif
    allocate(xdumd(ixrtd,jxrtd))
    allocate(idumd(ixrtd,jxrtd))
    xdumd = -9999.
    idumd = -9999

    if (output_rt2d .eq. 0) return  ! return if hires flag eq 0, if =1 output full grid

    if (output_count == 0) then
     !-- create the fine grid routing file
     if (outname == "") then
        write(output_flnm, '(A12,"_RT2D_DOMAIN",I1)') date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16), igrid
     else
        write(output_flnm,'(A, I1, A, A12)') trim(outname), igrid,"_", date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16)
     endif
#ifdef HYDRO_D
     print*, 'output_flnm = "'//trim(output_flnm)//'"'
#endif

#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
       iret = nf_create(trim(output_flnm), IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
#else
       iret = nf_create(trim(output_flnm), NF_CLOBBER, ncid)
#endif
#ifdef HYDRO_D
       if (iret /= 0) then
         print*, "Problem nf_create"
         call hydro_stop("output_rt")
       endif
#endif

       iret = nf_def_dim(ncid, "time",  NF_UNLIMITED, dimid_times)
       iret = nf_def_dim(ncid, "x",     ixrtd,        dimid_ix)    !-- make a decimated grid
       iret = nf_def_dim(ncid, "y",     jxrtd,        dimid_jx)
       iret = nf_def_dim(ncid, "depth", nsoil,        dimid_soil)  !-- 3-d soils
       !C.Jiang, add to output sediment variable
       iret = nf_def_dim(ncid, "numfracts", NumFracts, dimid_sFracts)  !-- 3 soil fraction

       !--- define variables
       !- time definition, timeObs
       iret = nf_def_var(ncid,"time",NF_INT, 1, (/dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',34,sec_since_date)

       if(present(ltype)) then
       iret = nf_def_var(ncid,"grid_type",NF_INT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',8,'unitless')
       iret = nf_put_att_text(ncid,varid,'long_name',56,'grid type(1-land basin, 2-ocean basin, 0-ocean, -1-land)')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_INT,1,-9E15)
       endif

       if(present(elev)) then
       iret = nf_def_var(ncid,"elev",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',27,'elevation of ground surface')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(bank)) then
       iret = nf_def_var(ncid,"bank",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',20,'depth of stream bank')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(sabs)) then
       iret = nf_def_var(ncid,"sabs",NF_FLOAT, 3,(/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',3,'m/m')
       iret = nf_put_att_text(ncid,varid,'long_name',14,'sabs slope abs')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(width)) then
       iret = nf_def_var(ncid,"width",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',20,'width of stream bank')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(bot)) then
       iret = nf_def_var(ncid,"bot",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',27,'elevation of aquifer bottom')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(fbed)) then
       iret = nf_def_var(ncid,"f_bed",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',5,'m2/m2')
       iret = nf_put_att_text(ncid,varid,'long_name',28,'fractional area of streambed')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(gwhead) .and. present(elev)) then
       iret = nf_def_var(ncid,"d_gwhead",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',46,'depth of groundwater head below ground surface')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(runoff)) then
       iret = nf_def_var(ncid,"runoff",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',3,'m/s')
       iret = nf_put_att_text(ncid,varid,'long_name',19,'surface runoff rate')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif
 
       if(present(irg)) then
       iret = nf_def_var(ncid,"irg",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',3,'m/s')
       iret = nf_put_att_text(ncid,varid,'long_name',19,'irrigation water use rate')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(r_sat)) then
       iret = nf_def_var(ncid,"sat_runoff",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',3,'m/s')
       iret = nf_put_att_text(ncid,varid,'long_name',24,'saturation excess runoff rate')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(sinfil)) then
       iret = nf_def_var(ncid,"recharge",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',3,'m/s')
       iret = nf_put_att_text(ncid,varid,'long_name',28,'soil drainage (+ to aquifer)')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(tinfil)) then
       iret = nf_def_var(ncid,"infil",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',3,'m/s')
       iret = nf_put_att_text(ncid,varid,'long_name',24,'infiltration')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(fwater)) then
       iret = nf_def_var(ncid,"f_water",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',5,'m2/m2')
       iret = nf_put_att_text(ncid,varid,'long_name',41,'fractional area of water surface per grid')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(fsat)) then
       iret = nf_def_var(ncid,"f_sat",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',5,'m2/m2')
       iret = nf_put_att_text(ncid,varid,'long_name',42,'fractional area of saturated soil per grid')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(cuiflux)) then
       iret = nf_def_var(ncid,"cuiflux",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',3,'m/s')
       iret = nf_put_att_text(ncid,varid,'long_name',33,'flux from surface water to vadose')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(cgiflux)) then
       iret = nf_def_var(ncid,"cgiflux",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',3,'m/s')
       iret = nf_put_att_text(ncid,varid,'long_name',45,'flux between aquifer and river (+ to aquifer)')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(cglflux)) then
       iret = nf_def_var(ncid,"cglflux",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',3,'m/s')
       iret = nf_put_att_text(ncid,varid,'long_name',44,'flux between aquifer and lake (+ to aquifer)')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(hlake)) then
       iret = nf_def_var(ncid,"h_sfwater",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',19,'surface water level')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(hlake) .and. present(elev) .and. present(bank)) then
       iret = nf_def_var(ncid,"d_sfwater",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',22,'depth of surface water')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(streamflux)) then
       iret = nf_def_var(ncid,"streamflux",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',4,'m3/s')
       iret = nf_put_att_text(ncid,varid,'long_name',20,'stream flow per grid')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif
       if(present(veloc)) then 
       iret = nf_def_var(ncid,"velocity",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',3,'m/s')
       iret = nf_put_att_text(ncid,varid,'long_name',24,'stream velocity per grid')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if(present(qsKR)) then
       iret = nf_def_var(ncid,"qsKR",NF_FLOAT, 3,(/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'m3')
       iret = nf_put_att_text(ncid,varid,'long_name',62,'sediment transport capacity of overland flow using KR equation')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(qsEH)) then
       iret = nf_def_var(ncid,"qsEH",NF_FLOAT,4,(/dimid_sFracts,dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'m3')
       iret = nf_put_att_text(ncid,varid,'long_name',61,'sediment transport capacity of channel flow using EH equation')
       iret = nf_put_att_text(ncid,varid,'coordinates',4,'numfracts x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(T_SedConc)) then
       iret = nf_def_var(ncid,"T_SedConc",NF_FLOAT,3,(/dimid_ix,dimid_jx,dimid_times/),varid)
       iret = nf_put_att_text(ncid,varid,'units',4,'mg/l')
       iret = nf_put_att_text(ncid,varid,'long_name',28,'Total sediment concentration')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(SedConc)) then
       iret = nf_def_var(ncid,"SedConc",NF_FLOAT,4,(/dimid_sFracts,dimid_ix,dimid_jx,dimid_times/),varid)
       iret = nf_put_att_text(ncid,varid,'units',4,'mg/l')
       iret = nf_put_att_text(ncid,varid,'long_name',40,'Sediment concentraion of size fraction i')
       iret = nf_put_att_text(ncid,varid,'coordinates',13,'numfracts x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(SusSed)) then
       iret = nf_def_var(ncid,"SusSed",NF_FLOAT,4,(/dimid_sFracts,dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'m3')
       iret = nf_put_att_text(ncid,varid,'long_name',18,'Suspended seidment')
       iret = nf_put_att_text(ncid,varid,'coordinates',13,'numfracts x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(DepSed)) then
       iret = nf_def_var(ncid,"DepSed",NF_FLOAT,4,(/dimid_sFracts,dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'m3')
       iret = nf_put_att_text(ncid,varid,'long_name',18,'Deposited seidment')
       iret = nf_put_att_text(ncid,varid,'coordinates',13,'numfracts x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ErodSed)) then
       iret = nf_def_var(ncid,"ErodSed",NF_FLOAT,4,(/dimid_sFracts,dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',2,'m3')
       iret = nf_put_att_text(ncid,varid,'long_name',32,'Eroded seidment of overland flow')
       iret = nf_put_att_text(ncid,varid,'coordinates',13,'numfracts x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(SedFlux)) then
       iret = nf_def_var(ncid,"SedFlux",NF_FLOAT,4,(/dimid_sFracts,dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',4,'m3/s')
       iret = nf_put_att_text(ncid,varid,'long_name',13,'Sediment Flux')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(t_SedFlux)) then
       iret = nf_def_var(ncid,"t_SedFlux",NF_FLOAT,3,(/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',4,'m3/s')
       iret = nf_put_att_text(ncid,varid,'long_name',19,'Total Sediment Flux')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if(present(istreamflux)) then
       iret = nf_def_var(ncid,"istreamflux",NF_INT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',8,'unitless')
       iret = nf_put_att_text(ncid,varid,'long_name',27,'stream flow direction (1-8)')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_INT,1,-9999)
       endif

       if(present(w_us)) then
       iret = nf_def_var(ncid,"water_us",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',32,'water storage in unsaturated soil')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(w_gw)) then
       iret = nf_def_var(ncid,"water_gw",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',19,'groundwater storage')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(w_sf)) then
       iret = nf_def_var(ncid,"water_sf",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',31,'water storage in rivers & lakes')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_pmewsf)) then
       iret = nf_def_var(ncid,"ac_pme_wsf",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',49,'accumulated net rainfall (P-E) over water surface')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_infil_grnd)) then
       iret = nf_def_var(ncid,"ac_infil_grnd",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',42,'accumulated infiltration at ground surface')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_drain_grnd)) then
       iret = nf_def_var(ncid,"ac_drain_grnd",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',39,'accumulated soil drainage for soil column')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_roff_grnd)) then
       iret = nf_def_var(ncid,"ac_runoff_grnd",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',36,'accumulated runoff at ground surface')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_rsat_grnd)) then
       iret = nf_def_var(ncid,"ac_sat_runoff_grnd",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',54,'accumulated saturation-excess runoff at ground surface')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_precip)) then
       iret = nf_def_var(ncid,"ac_precip",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',36,'accumulated grid-mean precipitation')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_evap)) then
       iret = nf_def_var(ncid,"ac_evap",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',33,'accumulated grid-mean evaporation')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_cgi)) then
       iret = nf_def_var(ncid,"ac_cgi",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',53,'accumulated water between aquifer and river(+ to aq.)')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif 

       if(present(ac_cgl)) then
       iret = nf_def_var(ncid,"ac_cgl",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',52,'accumulated water between aquifer and lake(+ to aq.)')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_cui)) then
       iret = nf_def_var(ncid,"ac_cui",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',54,'accumulated water between vadose and sfwater(+ to va.)')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif
       if(present(ac_precevap)) then
       iret = nf_def_var(ncid,"ac_precevap",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',73,'accumulated precipitation minus evaporation on water surfaces(+ downward)')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_irgrw)) then 
       iret = nf_def_var(ncid,"ac_irgrw",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',39,'accumulated irrigation water from river')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif


       if(present(ac_irggw)) then 
       iret = nf_def_var(ncid,"ac_irggw",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',45,'accumulated irrigation water from groundwater')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif



       if(present(ac_ocngw)) then
       iret = nf_def_var(ncid,"ac_gw_ocn",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',32,'accumulated groundwater to ocean')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_ocnsf)) then
       iret = nf_def_var(ncid,"ac_sfw_ocn",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',34,'accumulated surface water to ocean')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_fixgw)) then
       iret = nf_def_var(ncid,"ac_fix_gw",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',41,'accumulated water to fix groundwater head')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       if(present(ac_fixsf)) then
       iret = nf_def_var(ncid,"ac_fix_sfw",NF_FLOAT, 3, (/dimid_ix,dimid_jx,dimid_times/), varid)
       iret = nf_put_att_text(ncid,varid,'units',1,'m')
       iret = nf_put_att_text(ncid,varid,'long_name',44,'accumulated water to fix surface water level')
       iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
       iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
       iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
       endif

       date19(1:19) = "0000-00-00_00:00:00"
       date19(1:len_trim(startdate)) = startdate
       convention(1:32) = "CF-1.0"
       iret = nf_put_att_real(ncid, NF_GLOBAL, "missing_value", NF_FLOAT, 1, -9E15)
       iret = nf_put_att_text(ncid, NF_GLOBAL, "Conventions",6, convention)
       iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_start",19, date19)
       iret = nf_put_att_int(ncid,  NF_GLOBAL, "output_decimation_factor",NF_INT, 1,decimation)

       iret = nf_enddef(ncid)

    endif

    output_count = output_count + 1

    !-- time
    iret = nf_inq_varid(ncid,"time", varid)
    iret = nf_put_vara_int(ncid, varid, (/output_count/), (/1/), seconds_since)

       if (present(ltype)) then
#ifdef MPP_LAND
          idumd = gLtype
#else  
          idumd = ltype
#endif
          iret = nf_inq_varid(ncid, "grid_type", varid)
          iret = nf_put_vara_int(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), idumd)
       endif

       if (present(elev)) then
#ifdef MPP_LAND
          xdumd = gElev
#else  
          xdumd = elev
#endif
          iret = nf_inq_varid(ncid, "elev", varid)
          iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
       endif

       if (present(bank)) then
#ifdef MPP_LAND
          xdumd = gBank
#else 
          xdumd = bank
#endif
          iret = nf_inq_varid(ncid, "bank", varid)
          iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
       endif

       if (present(sabs)) then
!#ifdef MPP_LAND
!          xdumd = sabs
!#else
          xdumd = sabs
!#endif
          iret = nf_inq_varid(ncid, "sabs", varid)
          iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/),(/ixrtd,jxrtd,1/), xdumd)
       endif


       if (present(width)) then
#ifdef MPP_LAND
          xdumd = gWidth
#else 
          xdumd = width
#endif
          iret = nf_inq_varid(ncid, "width", varid)
          iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
       endif


       if (present(bot)) then
#ifdef MPP_LAND
          xdumd = gBot
#else  
          xdumd = bot
#endif
          iret = nf_inq_varid(ncid, "bot", varid)
          iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
       endif

       if (present(fbed)) then
#ifdef MPP_LAND
          xdumd = gFbed
#else  
          xdumd = fbed
#endif
          iret = nf_inq_varid(ncid, "f_bed", varid)
          iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
       endif

    if (present(gwhead) .and. present(elev)) then
#ifdef MPP_LAND
     xdumd = gElev - gGwhead 
#else  
     xdumd = elev - gwhead
!C.Jiang,2020,
!      xdumd = topo - gwhead
#endif
     iret = nf_inq_varid(ncid, "d_gwhead", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(irg)) then
#ifdef MPP_LAND
     xdumd = girg
#else
     xdumd = irg
#endif
     iret = nf_inq_varid(ncid, "irg", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(runoff)) then
#ifdef MPP_LAND
     xdumd = gRunoff
#else
     xdumd = runoff
#endif
     iret = nf_inq_varid(ncid, "runoff", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(r_sat)) then
#ifdef MPP_LAND
     xdumd = gR_sat
#else  
     xdumd = r_sat
#endif
     iret = nf_inq_varid(ncid, "sat_runoff", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(sinfil)) then
#ifdef MPP_LAND
     xdumd = gSinfil
#else  
     xdumd = sinfil
#endif
     iret = nf_inq_varid(ncid, "recharge", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(tinfil)) then
#ifdef MPP_LAND
     xdumd = gTinfil
#else  
     xdumd = tinfil
#endif
     iret = nf_inq_varid(ncid, "infil", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(fwater)) then
#ifdef MPP_LAND
     xdumd = gFwater
#else  
     xdumd = fwater
#endif
     iret = nf_inq_varid(ncid, "f_water", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(fsat)) then
#ifdef MPP_LAND
     xdumd = gFsat
#else  
     xdumd = fsat
#endif
     iret = nf_inq_varid(ncid, "f_sat", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(cuiflux)) then
#ifdef MPP_LAND
     xdumd = gCuiflux
#else  
     xdumd = cuiflux
#endif
     iret = nf_inq_varid(ncid, "cuiflux", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(cgiflux)) then
#ifdef MPP_LAND
     xdumd = gCgiflux
#else  
     xdumd = cgiflux
#endif
     iret = nf_inq_varid(ncid, "cgiflux", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(cglflux)) then
#ifdef MPP_LAND
     xdumd = gCglflux
#else  
     xdumd = cglflux
#endif
     iret = nf_inq_varid(ncid, "cglflux", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(hlake)) then
#ifdef MPP_LAND
     xdumd = gHlake
#else  
     xdumd = hlake
#endif
     iret = nf_inq_varid(ncid, "h_sfwater", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if(present(hlake) .and. present(elev) .and. present(bank)) then
#ifdef MPP_LAND
     xdumd = gHlake - (gElev - gBank)
#else  
     xdumd = hlake - (elev - bank)
#endif
     iret = nf_inq_varid(ncid, "d_sfwater", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(streamflux)) then
#ifdef MPP_LAND
     xdumd = gStreamflux
#else  
     xdumd = streamflux
#endif
     iret = nf_inq_varid(ncid, "streamflux", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(veloc)) then
#ifdef MPP_LAND
     xdumd = gVeloc
#else  
     xdumd = veloc
#endif
     iret = nf_inq_varid(ncid, "velocity", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (present(qsKR)) then
#ifdef MPP_LAND
     xdumd = gqsKR
#else
     xdumd = qsKR
#endif
     iret = nf_inq_varid(ncid, "qsKR", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/),(/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(T_SedConc)) then
#ifdef MPP_LAND
     xdumd = gT_SedConc
#else
     xdumd = T_SedConc
#endif
     iret = nf_inq_varid(ncid, "T_SedConc", varid)
     iret = nf_put_vara_real(ncid, varid,(/1,1,output_count/),(/ixrtd,jxrtd,1/), xdumd)
    endif

    !-- 3-d soil numFracts
    if (present(SedConc)) then
      do n = 1, NumFracts
#ifdef MPP_LAND
     xdumd = gSedConc(n,:,:)
#else
     xdumd = SedConc(n,:,:)
#endif
     iret = nf_inq_varid(ncid, "SedConc", varid)
     iret = nf_put_vara_real(ncid,varid,(/n,1,1,output_count/),(/1,ixrtd,jxrtd,1/), xdumd)
      enddo !-n soil numFracts
    endif

    !-- 3-d soil numFracts
    if (present(qsEH)) then
      do n = 1, NumFracts
#ifdef MPP_LAND
     xdumd = gqsEH(n,:,:)
#else
     xdumd = qsEH(n,:,:)
#endif
     iret = nf_inq_varid(ncid, "qsEH", varid)
     iret = nf_put_vara_real(ncid,varid,(/n,1,1,output_count/),(/1,ixrtd,jxrtd,1/), xdumd)
      enddo !-n soil numFracts
    endif

    !-- 3-d soil numFracts
    if (present(SusSed)) then
      do n = 1, NumFracts
#ifdef MPP_LAND
     xdumd = gSusSed(n,:,:)
#else
     xdumd = SusSed(n,:,:)
#endif
     iret = nf_inq_varid(ncid, "SusSed", varid)
     iret = nf_put_vara_real(ncid, varid,(/n,1,1,output_count/),(/1,ixrtd,jxrtd,1/), xdumd)
      enddo !-n soil numFracts
    endif

    !--3-d soil numFracts
    if (present(DepSed)) then
      do n = 1, NumFracts
#ifdef MPP_LAND
     xdumd = gDepSed(n,:,:)
#else
     xdumd = DepSed(n,:,:)
#endif
     iret = nf_inq_varid(ncid, "DepSed", varid)
     iret = nf_put_vara_real(ncid, varid,(/n,1,1,output_count/),(/1,ixrtd,jxrtd,1/),xdumd)
      enddo !-n soil numFracts
    endif

    !--3-d soil numFracts
    if (present(ErodSed)) then
      do n = 1, NumFracts
#ifdef MPP_LAND
     xdumd = gErodSed(n,:,:)
#else
     xdumd = ErodSed(n,:,:)
#endif
     iret = nf_inq_varid(ncid, "ErodSed", varid)
     iret = nf_put_vara_real(ncid,varid,(/n,1,1,output_count/),(/1,ixrtd,jxrtd,1/),xdumd)
      enddo !-n soil numFracts
    endif

    !--3-d soil numFracts
    if (present(SedFlux)) then
      do n = 1, NumFracts
#ifdef MPP_LAND
     xdumd = gSedFlux(n,:,:)
#else
     xdumd = SedFlux(n,:,:)
#endif
     iret = nf_inq_varid(ncid, "SedFlux", varid)
     iret = nf_put_vara_real(ncid,varid,(/n,1,1,output_count/),(/1,ixrtd,jxrtd,1/),xdumd)
      enddo !-n soil numFracts
    endif

    if (present(t_SedFlux)) then
#ifdef MPP_LAND
     xdumd = gt_SedFlux
#else
     xdumd = t_SedFlux
#endif
     iret = nf_inq_varid(ncid, "t_SedFlux", varid)
     iret = nf_put_vara_real(ncid,varid,(/1,1,output_count/),(/ixrtd,jxrtd,1/),xdumd)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(istreamflux)) then
#ifdef MPP_LAND
     idumd = gIstreamflux
#else  
     idumd = istreamflux
#endif
     iret = nf_inq_varid(ncid, "istreamflux", varid)
     iret = nf_put_vara_int(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), idumd)
    endif

    if (present(w_us)) then
#ifdef MPP_LAND
     xdumd = gW_us
#else  
     xdumd = w_us
#endif
     iret = nf_inq_varid(ncid, "water_us", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(w_gw)) then
#ifdef MPP_LAND
     xdumd = gW_gw
#else  
     xdumd = w_gw
#endif
     iret = nf_inq_varid(ncid, "water_gw", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(w_sf)) then
#ifdef MPP_LAND
     xdumd = gW_sf
#else  
     xdumd = w_sf
#endif
     iret = nf_inq_varid(ncid, "water_sf", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_pmewsf)) then
#ifdef MPP_LAND
     xdumd = gAc_pmewsf
#else  
     xdumd = ac_pmewsf
#endif
     iret = nf_inq_varid(ncid, "ac_pme_wsf", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_infil_grnd)) then
#ifdef MPP_LAND
     xdumd = gAc_infil_grnd
#else  
     xdumd = ac_infil_grnd
#endif
     iret = nf_inq_varid(ncid, "ac_infil_grnd", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_drain_grnd)) then
#ifdef MPP_LAND
     xdumd = gAc_drain_grnd
#else 
     xdumd = ac_drain_grnd
#endif
     iret = nf_inq_varid(ncid, "ac_drain_grnd", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_roff_grnd)) then
#ifdef MPP_LAND
     xdumd = gAc_roff_grnd
#else  
     xdumd = ac_roff_grnd
#endif
     iret = nf_inq_varid(ncid, "ac_runoff_grnd", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_rsat_grnd)) then
#ifdef MPP_LAND
     xdumd = gAc_rsat_grnd
#else  
     xdumd = ac_rsat_grnd
#endif
     iret = nf_inq_varid(ncid, "ac_sat_runoff_grnd", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_precip)) then
#ifdef MPP_LAND
     xdumd = gAc_precip
#else  
     xdumd = ac_precip
#endif
     iret = nf_inq_varid(ncid, "ac_precip", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_evap)) then
#ifdef MPP_LAND
     xdumd = gAc_evap
#else  
     xdumd = ac_evap
#endif
     iret = nf_inq_varid(ncid, "ac_evap", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_irgrw)) then
#ifdef MPP_LAND
!     xdumd = gAc_evap
#else  
     xdumd = ac_irgrw
#endif
     iret = nf_inq_varid(ncid, "ac_irgrw", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_irggw)) then
#ifdef MPP_LAND
!     xdumd = gAc_evap
#else  
     xdumd = ac_irggw
#endif
     iret = nf_inq_varid(ncid, "ac_irggw", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_cgi)) then
#ifdef MPP_LAND
     xdumd = gAc_cgi
#else  
     xdumd = ac_cgi
#endif
     iret = nf_inq_varid(ncid, "ac_cgi", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_cgl)) then
#ifdef MPP_LAND
     xdumd = gAc_cgl
#else  
     xdumd = ac_cgl
#endif
     iret = nf_inq_varid(ncid, "ac_cgl", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_cui)) then
#ifdef MPP_LAND
     xdumd = gAc_cui
#else  
     xdumd = ac_cui
#endif
     iret = nf_inq_varid(ncid, "ac_cui", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_precevap)) then
#ifdef MPP_LAND
     xdumd = gAc_precevap
#else
     xdumd = ac_precevap
#endif
     iret = nf_inq_varid(ncid, "ac_precevap", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif


    if (present(ac_ocngw)) then
#ifdef MPP_LAND
     xdumd = gAc_ocngw
#else  
     xdumd = ac_ocngw
#endif
     iret = nf_inq_varid(ncid, "ac_gw_ocn", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_ocnsf)) then
#ifdef MPP_LAND
     xdumd = gAc_ocnsf
#else  
     xdumd = ac_ocnsf
#endif
     iret = nf_inq_varid(ncid, "ac_sfw_ocn", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_fixgw)) then
#ifdef MPP_LAND
     xdumd = gAc_fixgw
#else  
     xdumd = ac_fixgw
#endif
     iret = nf_inq_varid(ncid, "ac_fix_gw", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    if (present(ac_fixsf)) then
#ifdef MPP_LAND
     xdumd = gAc_fixsf
#else  
     xdumd = ac_fixsf
#endif
     iret = nf_inq_varid(ncid, "ac_fix_sfw", varid)
     iret = nf_put_vara_real(ncid, varid, (/1,1,output_count/), (/ixrtd,jxrtd,1/), xdumd)
    endif

    ! time in seconds since startdate
    iret = nf_redef(ncid)
    date19(1:len_trim(date)) = date
    iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_end", 19, date19)

    iret = nf_enddef(ncid)
    iret = nf_sync(ncid)
    if (output_count == split_output_count) then
        output_count = 0
        iret = nf_close(ncid)
    endif

     deallocate(xdumd)
     deallocate(idumd)

#ifdef HYDRO_D 
    write(6,*) "finish output for rt2d"
#endif

#ifdef MPP_LAND
    endif      ! if (my_id .eq. IO_id) then
#endif

  end subroutine sub_output_rt2d

!-- output the channel route in an IDV 'station' compatible format
   subroutine output_chrt(igrid, split_output_count, NLINKS, ORDER, &
        startdate,date,chlon, chlat, hlink,zelev,qlink,dtrt,K,   &
        STRMFRXSTPTS,order_to_write)
     
     implicit none
#include <netcdf.inc>
!!output the routing variables over just channel
     integer,                                  intent(in) :: igrid,K
     integer,                                  intent(in) :: split_output_count
     integer,                                  intent(in) :: NLINKS
     real, dimension(NLINKS),                  intent(in) :: chlon,chlat
     real, dimension(NLINKS),                  intent(in) :: hlink,zelev
     integer, dimension(NLINKS),               intent(in) :: ORDER
     integer, dimension(NLINKS),               intent(inout) :: STRMFRXSTPTS

     real,                                     intent(in) :: dtrt
     real, dimension(NLINKS,2),                intent(in) :: qlink

     character(len=*),                         intent(in) :: startdate
     character(len=*),                         intent(in) :: date

     real, allocatable, DIMENSION(:)            :: chanlat,chanlon
     real, allocatable, DIMENSION(:)            :: chanlatO,chanlonO

     real, allocatable, DIMENSION(:)            :: elevation
     real, allocatable, DIMENSION(:)            :: elevationO

     integer, allocatable, DIMENSION(:)         :: station_id
     integer, allocatable, DIMENSION(:)         :: station_idO

     integer, allocatable, DIMENSION(:)         :: rec_num_of_station
     integer, allocatable, DIMENSION(:)         :: rec_num_of_stationO

     integer, allocatable, DIMENSION(:)         :: lOrder !- local stream order
     integer, allocatable, DIMENSION(:)         :: lOrderO !- local stream order

     integer, save  :: output_count
     integer, save  :: ncid,ncid2

     integer :: stationdim, dimdata, varid, charid, n
     integer :: obsdim, dimdataO, charidO

     integer :: iret,i, start_pos, prev_pos, order_to_write!-- order_to_write is the lowest stream order to output
     integer :: start_posO, prev_posO

     integer :: previous_pos  !-- used for the station model
     character(len=256) :: output_flnm,output_flnm2
     character(len=19)  :: date19,date19start
     character(len=34)  :: sec_since_date
     integer :: seconds_since,nstations,cnt,ObsStation,nobs
     character(len=32)  :: convention
     character(len=11),allocatable, DIMENSION(:)  :: stname
     character(len=11),allocatable, DIMENSION(:)  :: stnameO

    !--- all this for writing the station id string
     INTEGER   TDIMS, TXLEN
     PARAMETER (TDIMS=2)    ! number of TX dimensions
     PARAMETER (TXLEN = 11) ! length of example string
     INTEGER  TIMEID        ! record dimension id
     INTEGER  TXID          ! variable ID
     INTEGER  TXDIMS(TDIMS) ! variable shape
     INTEGER  TSTART(TDIMS), TCOUNT(TDIMS)

     !--  observation point  ids
     INTEGER   OTDIMS, OTXLEN
     PARAMETER (OTDIMS=2)    ! number of TX dimensions
     PARAMETER (OTXLEN = 11) ! length of example string
     INTEGER  OTIMEID        ! record dimension id
     INTEGER  OTXID          ! variable ID
     INTEGER  OTXDIMS(OTDIMS) ! variable shape
     INTEGER  OTSTART(OTDIMS), OTCOUNT(OTDIMS)

#ifdef HYDRO_D
     write(6,*) "yyww dtrt =", dtrt , "k =", k
#endif
   
     seconds_since = int(dtrt)*K

!    order_to_write = 2  !-- 1 all; 6 feweest

      nstations = 0  ! total number of channel points to display
      nobs      = 0  ! number of observation points

!-- output only the higher oder streamflows  and only observation points
       do i=1,NLINKS
        if(ORDER(i) .ge. order_to_write) then 
         nstations = nstations + 1
        endif
        if(STRMFRXSTPTS(i) .ne. -9999) then 
         nobs = nobs + 1
        endif
       enddo 

       if (nobs .eq. 0) then ! let's at least make one obs point
         nobs = 1
         STRMFRXSTPTS(1) = 1
       endif

       allocate(chanlat(nstations))
       allocate(chanlon(nstations))
       allocate(elevation(nstations))
       allocate(station_id(nstations))
       allocate(lOrder(nstations))
       allocate(rec_num_of_station(nstations))
       allocate(stname(nstations))

       allocate(chanlatO(nobs))
       allocate(chanlonO(nobs))
       allocate(elevationO(nobs))
       allocate(station_idO(nobs))
       allocate(lOrderO(nobs))
       allocate(rec_num_of_stationO(nobs))
       allocate(stnameO(nobs))

       if(output_count == 0) then 
!-- have moved sec_since_date from above here..
        sec_since_date = 'seconds since '//startdate(1:4)//'-'//startdate(6:7)//'-'//startdate(9:10) &
                  //' '//startdate(12:13)//':'//startdate(15:16)//' UTC'

        date19start(1:len_trim(startdate)) = startdate(1:4)//'-'//startdate(6:7)//'-'//startdate(9:10)//'_' &
                  //startdate(12:13)//':'//startdate(15:16)//':00'

        nstations = 0
        nobs = 0

        write(output_flnm, '(A12,".CHRTOUT_DOMAIN",I1)') date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16), igrid
        write(output_flnm2,'(A12,".CHANOBS_DOMAIN",I1)') date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16), igrid

#ifdef HYDRO_D
        print*, 'output_flnm = "'//trim(output_flnm)//'"'
#endif

#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
       iret = nf_create(trim(output_flnm), IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
#else
       iret = nf_create(trim(output_flnm), NF_CLOBBER, ncid)
#endif
#ifdef HYDRO_D
        if (iret /= 0) then
           print*,  "Problem nf_create points"
           call hydro_stop("output_chrt")
        endif
#endif

#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
       iret = nf_create(trim(output_flnm2), IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid2)
#else
       iret = nf_create(trim(output_flnm2), NF_CLOBBER, ncid2)
#endif
#ifdef HYDRO_D
        if (iret /= 0) then
            print*, "Problem nf_create observation"
            call hydro_stop("output_chrt")
        endif
#endif

       do i=1,NLINKS
        if(ORDER(i) .ge. order_to_write) then 
         nstations = nstations + 1
         chanlat(nstations) = chlat(i)
         chanlon(nstations) = chlon(i)
         elevation(nstations) = zelev(i)
         lOrder(nstations) = ORDER(i)
         station_id(nstations) = i
         if(STRMFRXSTPTS(nstations) .eq. -9999) then 
           ObsStation = 0
         else 
           ObsStation = 1
         endif
         write(stname(nstations),'(I6,"_",I1,"_S",I1)') nstations,lOrder(nstations),ObsStation
        endif
       enddo 


       do i=1,NLINKS
        if(STRMFRXSTPTS(i) .ne. -9999) then 
         nobs = nobs + 1
         chanlatO(nobs) = chlat(i)
         chanlonO(nobs) = chlon(i)
         elevationO(nobs) = zelev(i)
         lOrderO(nobs) = ORDER(i)
         station_idO(nobs) = i
         write(stnameO(nobs),'(I6,"_",I1)') nobs,lOrderO(nobs)
#ifdef HYDRO_D
!        print *,"stationobservation name",  stnameO(nobs)
#endif
        endif
       enddo 

       iret = nf_def_dim(ncid, "recNum", NF_UNLIMITED, dimdata)  !--for linked list approach


       iret = nf_def_dim(ncid, "station", nstations, stationdim)



       iret = nf_def_dim(ncid2, "recNum", NF_UNLIMITED, dimdataO)  !--for linked list approach
       iret = nf_def_dim(ncid2, "station", nobs, obsdim)


      !- station location definition all,  lat
        iret = nf_def_var(ncid,"latitude",NF_FLOAT, 1, (/stationdim/), varid)
#ifdef HYDRO_D
       write(6,*) "iret 2.1,  ", iret, stationdim
#endif
        iret = nf_put_att_text(ncid,varid,'long_name',16,'Station latitude')
#ifdef HYDRO_D
       write(6,*) "iret 2.2", iret
#endif
        iret = nf_put_att_text(ncid,varid,'units',13,'degrees_north')
#ifdef HYDRO_D
       write(6,*) "iret 2.3", iret
#endif


      !- station location definition obs,  lat
        iret = nf_def_var(ncid2,"latitude",NF_FLOAT, 1, (/obsdim/), varid)
        iret = nf_put_att_text(ncid2,varid,'long_name',20,'Observation latitude')
        iret = nf_put_att_text(ncid2,varid,'units',13,'degrees_north')


      !- station location definition,  long
        iret = nf_def_var(ncid,"longitude",NF_FLOAT, 1, (/stationdim/), varid)
        iret = nf_put_att_text(ncid,varid,'long_name',17,'Station longitude')
        iret = nf_put_att_text(ncid,varid,'units',12,'degrees_east')


      !- station location definition, obs long
        iret = nf_def_var(ncid2,"longitude",NF_FLOAT, 1, (/obsdim/), varid)
        iret = nf_put_att_text(ncid2,varid,'long_name',21,'Observation longitude')
        iret = nf_put_att_text(ncid2,varid,'units',12,'degrees_east')


!     !-- elevation is ZELEV
        iret = nf_def_var(ncid,"altitude",NF_FLOAT, 1, (/stationdim/), varid)
        iret = nf_put_att_text(ncid,varid,'long_name',16,'Station altitude')
        iret = nf_put_att_text(ncid,varid,'units',6,'meters')


!     !-- elevation is obs ZELEV
        iret = nf_def_var(ncid2,"altitude",NF_FLOAT, 1, (/obsdim/), varid)
        iret = nf_put_att_text(ncid2,varid,'long_name',20,'Observation altitude')
        iret = nf_put_att_text(ncid2,varid,'units',6,'meters')


!     !--  gage observation
!       iret = nf_def_var(ncid,"gages",NF_FLOAT, 1, (/stationdim/), varid)
!       iret = nf_put_att_text(ncid,varid,'long_name',20,'Stream Gage Location')
!       iret = nf_put_att_text(ncid,varid,'units',4,'none')

!-- parent index
        iret = nf_def_var(ncid,"parent_index",NF_INT,1,(/dimdata/), varid)
        iret = nf_put_att_text(ncid,varid,'long_name',36,'index of the station for this record')

        iret = nf_def_var(ncid2,"parent_index",NF_INT,1,(/dimdataO/), varid)
        iret = nf_put_att_text(ncid2,varid,'long_name',36,'index of the station for this record')

     !-- prevChild
        iret = nf_def_var(ncid,"prevChild",NF_INT,1,(/dimdata/), varid)
        iret = nf_put_att_text(ncid,varid,'long_name',57,'record number of the previous record for the same station')
!ywtmp        iret = nf_put_att_int(ncid,varid,'_FillValue',NF_INT,2,-1)
        iret = nf_put_att_int(ncid,varid,'_FillValue',2,-1)

        iret = nf_def_var(ncid2,"prevChild",NF_INT,1,(/dimdataO/), varid)
        iret = nf_put_att_text(ncid2,varid,'long_name',57,'record number of the previous record for the same station')
!ywtmp        iret = nf_put_att_int(ncid2,varid,'_FillValue',NF_INT,2,-1)
        iret = nf_put_att_int(ncid2,varid,'_FillValue',2,-1)

     !-- lastChild
        iret = nf_def_var(ncid,"lastChild",NF_INT,1,(/stationdim/), varid)
        iret = nf_put_att_text(ncid,varid,'long_name',30,'latest report for this station')
!ywtmp        iret = nf_put_att_int(ncid,varid,'_FillValue',NF_INT,2,-1)
        iret = nf_put_att_int(ncid,varid,'_FillValue',2,-1)

        iret = nf_def_var(ncid2,"lastChild",NF_INT,1,(/obsdim/), varid)
        iret = nf_put_att_text(ncid2,varid,'long_name',30,'latest report for this station')
!ywtmp        iret = nf_put_att_int(ncid2,varid,'_FillValue',NF_INT,2,-1)
        iret = nf_put_att_int(ncid2,varid,'_FillValue',2,-1)

!     !- flow definition, var
        iret = nf_def_var(ncid, "streamflow", NF_FLOAT, 1, (/dimdata/), varid)
        iret = nf_put_att_text(ncid,varid,'units',13,'meter^3 / sec')
        iret = nf_put_att_text(ncid,varid,'long_name',10,'River Flow')

        iret = nf_def_var(ncid2, "streamflow", NF_FLOAT, 1, (/dimdataO/), varid)
        iret = nf_put_att_text(ncid2,varid,'units',13,'meter^3 / sec')
        iret = nf_put_att_text(ncid2,varid,'long_name',10,'River Flow')

!     !- flow definition, var
!       iret = nf_def_var(ncid, "pos_streamflow", NF_FLOAT, 1, (/dimdata/), varid)
!       iret = nf_put_att_text(ncid,varid,'units',13,'meter^3 / sec')
!       iret = nf_put_att_text(ncid,varid,'long_name',14,'abs streamflow')

!     !- head definition, var
        iret = nf_def_var(ncid, "head", NF_FLOAT, 1, (/dimdata/), varid)
        iret = nf_put_att_text(ncid,varid,'units',5,'meter')
        iret = nf_put_att_text(ncid,varid,'long_name',11,'River Stage')

        iret = nf_def_var(ncid2, "head", NF_FLOAT, 1, (/dimdataO/), varid)
        iret = nf_put_att_text(ncid2,varid,'units',5,'meter')
        iret = nf_put_att_text(ncid2,varid,'long_name',11,'River Stage')


!     !- order definition, var
        iret = nf_def_var(ncid, "order", NF_INT, 1, (/dimdata/), varid)
        iret = nf_put_att_text(ncid,varid,'long_name',21,'Strahler Stream Order')
        iret = nf_put_att_int(ncid,varid,'_FillValue',2,-1)

        iret = nf_def_var(ncid2, "order", NF_INT, 1, (/dimdataO/), varid)
        iret = nf_put_att_text(ncid2,varid,'long_name',21,'Strahler Stream Order')
        iret = nf_put_att_int(ncid,varid,'_FillValue',2,-1)


     !-- station  id
     ! define character-position dimension for strings of max length 11
         iret = NF_DEF_DIM(ncid, "id_len", 11, charid)
         TXDIMS(1) = charid   ! define char-string variable and position dimension first
         TXDIMS(2) = stationdim
         iret = nf_def_var(ncid,"station_id",NF_CHAR, TDIMS, TXDIMS, varid)
         iret = nf_put_att_text(ncid,varid,'long_name',10,'Station id')

         iret = NF_DEF_DIM(ncid2, "id_len", 11, charidO)
         OTXDIMS(1) = charidO   ! define char-string variable and position dimension first
         OTXDIMS(2) = obsdim
         iret = nf_def_var(ncid2,"station_id",NF_CHAR, OTDIMS, OTXDIMS, varid)
         iret = nf_put_att_text(ncid2,varid,'long_name',14,'Observation id')


!     !- time definition, timeObs
         iret = nf_def_var(ncid,"time_observation",NF_INT, 1, (/dimdata/), varid)
         iret = nf_put_att_text(ncid,varid,'units',34,sec_since_date)
         iret = nf_put_att_text(ncid,varid,'long_name',19,'time of observation')

         iret = nf_def_var(ncid2,"time_observation",NF_INT, 1, (/dimdataO/), varid)
         iret = nf_put_att_text(ncid2,varid,'units',34,sec_since_date)
         iret = nf_put_att_text(ncid2,varid,'long_name',19,'time of observation')

         iret = nf_put_att_text(ncid, NF_GLOBAL, "Conventions",32, convention)
         iret = nf_put_att_text(ncid2, NF_GLOBAL, "Conventions",32, convention)

         convention(1:32) = "Unidata Observation Dataset v1.0"
         iret = nf_put_att_text(ncid, NF_GLOBAL, "Conventions",32, convention)
         iret = nf_put_att_text(ncid, NF_GLOBAL, "cdm_datatype",7, "Station")
         iret = nf_put_att_text(ncid, NF_GLOBAL, "geospatial_lat_max",4, "90.0")
         iret = nf_put_att_text(ncid, NF_GLOBAL, "geospatial_lat_min",5, "-90.0")
         iret = nf_put_att_text(ncid, NF_GLOBAL, "geospatial_lon_max",5, "180.0")
         iret = nf_put_att_text(ncid, NF_GLOBAL, "geospatial_lon_min",6, "-180.0")
         iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_start",19, date19start)
         iret = nf_put_att_text(ncid, NF_GLOBAL, "stationDimension",7, "station")
         iret = nf_put_att_real(ncid, NF_GLOBAL, "missing_value", NF_FLOAT, 1, -9E15)
         iret = nf_put_att_int(ncid, NF_GLOBAL, "stream order output",NF_INT,1,order_to_write)

         iret = nf_put_att_text(ncid2, NF_GLOBAL, "Conventions",32, convention)
         iret = nf_put_att_text(ncid2, NF_GLOBAL, "cdm_datatype",7, "Station")
         iret = nf_put_att_text(ncid2, NF_GLOBAL, "geospatial_lat_max",4, "90.0")
         iret = nf_put_att_text(ncid2, NF_GLOBAL, "geospatial_lat_min",5, "-90.0")
         iret = nf_put_att_text(ncid2, NF_GLOBAL, "geospatial_lon_max",5, "180.0")
         iret = nf_put_att_text(ncid2, NF_GLOBAL, "geospatial_lon_min",6, "-180.0")
         iret = nf_put_att_text(ncid2, NF_GLOBAL, "time_coverage_start",19, date19start)
         iret = nf_put_att_text(ncid2, NF_GLOBAL, "stationDimension",7, "station")
         iret = nf_put_att_real(ncid2, NF_GLOBAL, "missing_value", NF_FLOAT, 1, -9E15)
         iret = nf_put_att_int(ncid2, NF_GLOBAL, "stream order output",NF_INT,1,order_to_write)

         iret = nf_enddef(ncid)
         iret = nf_enddef(ncid2)
        
        !-- write latitudes
         iret = nf_inq_varid(ncid,"latitude", varid)
         iret = nf_put_vara_real(ncid, varid, (/1/), (/nstations/), chanlat)

         iret = nf_inq_varid(ncid2,"latitude", varid)
         iret = nf_put_vara_real(ncid2, varid, (/1/), (/nobs/), chanlatO)

        !-- write longitudes
         iret = nf_inq_varid(ncid,"longitude", varid)
         iret = nf_put_vara_real(ncid, varid, (/1/), (/nstations/), chanlon)

         iret = nf_inq_varid(ncid2,"longitude", varid)
         iret = nf_put_vara_real(ncid2, varid, (/1/), (/nobs/), chanlonO)

        !-- write elevations
         iret = nf_inq_varid(ncid,"altitude", varid)
         iret = nf_put_vara_real(ncid, varid, (/1/), (/nstations/), elevation)

         iret = nf_inq_varid(ncid2,"altitude", varid)
         iret = nf_put_vara_real(ncid2, varid, (/1/), (/nobs/), elevationO)

      !-- write gage location
!      iret = nf_inq_varid(ncid,"gages", varid)
!      iret = nf_put_vara_int(ncid, varid, (/1/), (/nstations/), STRMFRXSTPTS)

        !-- write number_of_stations, OPTIONAL
      !!  iret = nf_inq_varid(ncid,"number_stations", varid)
      !!  iret = nf_put_var_int(ncid, varid, nstations)

        !-- write station id's 
         do i=1,nstations
          TSTART(1) = 1
          TSTART(2) = i
          TCOUNT(1) = TXLEN
          TCOUNT(2) = 1
          iret = nf_inq_varid(ncid,"station_id", varid)
          iret = nf_put_vara_text(ncid, varid, TSTART, TCOUNT, stname(i))
         enddo

        !-- write observation id's 
         do i=1, nobs
          OTSTART(1) = 1
          OTSTART(2) = i
          OTCOUNT(1) = OTXLEN
          OTCOUNT(2) = 1
          iret = nf_inq_varid(ncid2,"station_id", varid)
          iret = nf_put_vara_text(ncid2, varid, OTSTART, OTCOUNT, stnameO(i))
         enddo

     endif

     output_count = output_count + 1

     open (unit=999,file='frxst_pts_out.txt',status='unknown',position='append')

     cnt=0
     do i=1,NLINKS

       if(ORDER(i) .ge. order_to_write) then 
         start_pos = (cnt+1)+(nstations*(output_count-1))

         !!--time in seconds since startdate
          iret = nf_inq_varid(ncid,"time_observation", varid)
          iret = nf_put_vara_int(ncid, varid, (/start_pos/), (/1/), seconds_since)

         iret = nf_inq_varid(ncid,"streamflow", varid)
         iret = nf_put_vara_real(ncid, varid, (/start_pos/), (/1/), qlink(i,1))

!        iret = nf_inq_varid(ncid,"pos_streamflow", varid)
!        iret = nf_put_vara_real(ncid, varid, (/start_pos/), (/1/), abs(qlink(i,1)))

         iret = nf_inq_varid(ncid,"head", varid)
         iret = nf_put_vara_real(ncid, varid, (/start_pos/), (/1/), hlink(i))

         iret = nf_inq_varid(ncid,"order", varid)
         iret = nf_put_vara_int(ncid, varid, (/start_pos/), (/1/), ORDER(i))

         !-- station index.. will repeat for every timesstep
         iret = nf_inq_varid(ncid,"parent_index", varid)
         iret = nf_put_vara_int(ncid, varid, (/start_pos/), (/1/), cnt)

          !--record number of previous record for same station
!obsolete format         prev_pos = cnt+(nstations*(output_count-1))
         prev_pos = cnt+(nobs*(output_count-2))
         if(output_count.ne.1) then !-- only write next set of records
           iret = nf_inq_varid(ncid,"prevChild", varid)
           iret = nf_put_vara_int(ncid, varid, (/start_pos/), (/1/), prev_pos)
         endif

         cnt=cnt+1  !--indices are 0 based
         rec_num_of_station(cnt) = start_pos-1  !-- save position for last child, 0-based!!


       endif
    enddo
!    close(999) 

    !-- output  only observation points
    cnt=0
    do i=1,NLINKS

      if(STRMFRXSTPTS(i) .ne. -9999) then 
         start_posO = (cnt+1)+(nobs * (output_count-1))
!Write frxst_pts to text file...
!yw          write(999,117) seconds_since,trim(date),cnt,chlon(i),chlat(i), &
           write(999,117) seconds_since,date(1:10),date(12:19),cnt,chlon(i),chlat(i), &
                    qlink(i,1), qlink(i,1)*35.315,hlink(i)
!yw 117 FORMAT(I8,1X,A25,1X,I7,1X,F10.5,1X,F8.5,1X,F9.3,1x,F12.3,1X,F6.3)
!yw 117 FORMAT(I8,1X,A10,1X,A8,1x,I7,1X,F10.5,1X,F8.5,1X,F9.3,1x,F12.3,1X,F6.3)
  117 FORMAT(I8,",",A10,1X,A8,",",I7,",",F10.5,",",F8.5,",",F9.3,",",F12.3,",",F6.3)

         !!--time in seconds since startdate
         iret = nf_inq_varid(ncid2,"time_observation", varid)
         iret = nf_put_vara_int(ncid2, varid, (/start_posO/), (/1/), seconds_since)

         iret = nf_inq_varid(ncid2,"streamflow", varid)
         iret = nf_put_vara_real(ncid2, varid, (/start_posO/), (/1/), qlink(i,1))

         iret = nf_inq_varid(ncid2,"head", varid)
         iret = nf_put_vara_real(ncid2, varid, (/start_posO/), (/1/), hlink(i))

         iret = nf_inq_varid(ncid,"order", varid)
         iret = nf_put_vara_int(ncid2, varid, (/start_posO/), (/1/), ORDER(i))

         !-- station index.. will repeat for every timesstep
         iret = nf_inq_varid(ncid2,"parent_index", varid)
         iret = nf_put_vara_int(ncid2, varid, (/start_posO/), (/1/), cnt)

          !--record number of previous record for same station
!obsolete format          prev_posO = cnt+(nobs*(output_count-1))
         prev_posO = cnt+(nobs*(output_count-2))
         if(output_count.ne.1) then !-- only write next set of records
           iret = nf_inq_varid(ncid2,"prevChild", varid)
           iret = nf_put_vara_int(ncid2, varid, (/start_posO/), (/1/), prev_posO)

!IF block to add -1 to last element of prevChild array to designate end of list...
!           if(cnt+1.eq.nobs.AND.output_count.eq.split_output_count) then
!             iret = nf_put_vara_int(ncid2, varid, (/start_posO/), (/1/), -1)
!           else
!             iret = nf_put_vara_int(ncid2, varid, (/start_posO/), (/1/), prev_posO)
!           endif 


         endif

         cnt=cnt+1  !--indices are 0 based
         rec_num_of_stationO(cnt) = start_posO - 1  !-- save position for last child, 0-based!!


      endif

    enddo
    close(999) 


      !-- lastChild variable gives the record number of the most recent report for the station
      iret = nf_inq_varid(ncid,"lastChild", varid)
      iret = nf_put_vara_int(ncid, varid, (/1/), (/nstations/), rec_num_of_station)

      !-- lastChild variable gives the record number of the most recent report for the station
      iret = nf_inq_varid(ncid2,"lastChild", varid)
      iret = nf_put_vara_int(ncid2, varid, (/1/), (/nobs/), rec_num_of_stationO)

      iret = nf_redef(ncid)
      date19(1:19) = "0000-00-00_00:00:00"
      date19(1:len_trim(date)) = date
      iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_end", 19, date19)

      iret = nf_redef(ncid2)
      iret = nf_put_att_text(ncid2, NF_GLOBAL, "time_coverage_end", 19, date19)

      iret = nf_enddef(ncid)
      iret = nf_sync(ncid)

      iret = nf_enddef(ncid2)
      iret = nf_sync(ncid2)

      if (output_count == split_output_count) then
        output_count = 0
        iret = nf_close(ncid)
        iret = nf_close(ncid2)
     endif

     deallocate(chanlat)
     deallocate(chanlon)
     deallocate(elevation)
     deallocate(station_id)
     deallocate(lOrder)
     deallocate(rec_num_of_station)
     deallocate(stname)

     deallocate(chanlatO)
     deallocate(chanlonO)
     deallocate(elevationO)
     deallocate(station_idO)
     deallocate(lOrderO)
     deallocate(rec_num_of_stationO)
     deallocate(stnameO)
#ifdef HYDRO_D
     print *, "Exited Subroutine output_chrt"
#endif
     close(16)

20 format(i8,',',f12.7,',',f10.7,',',f6.2,',',i3)

end subroutine output_chrt

#ifdef MPP_LAND
!-- output the channel route in an IDV 'station' compatible format
   subroutine mpp_output_chrt(gnlinks,map_l2g,igrid, &
        split_output_count, NLINKS, ORDER, &
        startdate, date, chlon, chlat, hlink,zelev,qlink,dtrt, &
        K,STRMFRXSTPTS,order_to_write)

       USE module_mpp_land

!!output the routing variables over just channel
     integer,                                  intent(in) :: igrid,K
     integer,                                  intent(in) :: split_output_count
     integer,                                  intent(in) :: NLINKS
     real, dimension(NLINKS),               intent(in) :: chlon,chlat
     real, dimension(NLINKS),                  intent(in) :: hlink,zelev

     integer, dimension(NLINKS),               intent(in) :: ORDER
     integer, dimension(NLINKS),               intent(inout) :: STRMFRXSTPTS

     real,                                     intent(in) :: dtrt
     real, dimension(NLINKS,2),                intent(in) :: qlink

     character(len=*),                         intent(in) :: startdate
     character(len=*),                         intent(in) :: date

     integer  :: gnlinks, map_l2g(nlinks), order_to_write
     !real, dimension(gNLINKS) :: g_chlon,g_chlat, g_hlink,g_zelev
     !real, dimension(gNLINKS,2) :: g_qlink
     !integer ,  dimension(gNLINKS) :: g_order,g_STRMFRXSTPTS
      real, allocatable,dimension(:) :: g_chlon,g_chlat, g_hlink,g_zelev
      integer, allocatable,dimension(:) :: g_order,g_STRMFRXSTPTS
      real,allocatable,dimension(:,:) :: g_qlink

     if(my_id .eq. io_id ) then
        allocate(g_chlon(gNLINKS))
        allocate(g_chlat(gNLINKS))
        allocate(g_hlink(gNLINKS))
        allocate(g_zelev(gNLINKS))
        allocate(g_qlink(gNLINKS,2))
        allocate(g_order(gNLINKS))
        allocate(g_STRMFRXSTPTS(gNLINKS))
     else
        allocate(g_chlon(1))
        allocate(g_chlat(1))
        allocate(g_hlink(1))
        allocate(g_zelev(1))
        allocate(g_qlink(1,2))
        allocate(g_order(1))
        allocate(g_STRMFRXSTPTS(1))
     endif

     call write_chanel_int(order,map_l2g,gnlinks,nlinks,g_order)
     call write_chanel_int(STRMFRXSTPTS,map_l2g,gnlinks,nlinks,g_STRMFRXSTPTS)
     call write_chanel_real(chlon,map_l2g,gnlinks,nlinks,g_chlon)
     call write_chanel_real(chlat,map_l2g,gnlinks,nlinks,g_chlat)
     call write_chanel_real(hlink,map_l2g,gnlinks,nlinks,g_hlink)
     call write_chanel_real(zelev,map_l2g,gnlinks,nlinks,g_zelev)
     call write_chanel_real(qlink(:,1),map_l2g,gnlinks,nlinks,g_qlink(:,1))
     call write_chanel_real(qlink(:,2),map_l2g,gnlinks,nlinks,g_qlink(:,2))

     if(my_id .eq. IO_id) then
       call output_chrt(igrid, split_output_count, GNLINKS, g_ORDER, &
          startdate, date, g_chlon, g_chlat, g_hlink,g_zelev,g_qlink,dtrt,K,&
          g_STRMFRXSTPTS,order_to_write)
    end if
    deallocate(g_order, g_STRMFRXSTPTS, g_chlon, g_chlat, g_hlink , g_zelev, g_qlink)

end subroutine mpp_output_chrt

!---------  lake netcdf output -----------------------------------------
!-- output the ilake info an IDV 'station' compatible format -----------
   subroutine mpp_output_lakes(lake_index,igrid, split_output_count, NLAKES, &
        startdate, date, latlake, lonlake, elevlake, &
        qlakei,qlakeo, resht,dtrt,K)

   USE module_mpp_land

!!output the routing variables over just channel
     integer,                                  intent(in) :: igrid, K
     integer,                                  intent(in) :: split_output_count
     integer,                                  intent(in) :: NLAKES
     real, dimension(NLAKES),                  intent(in) :: latlake,lonlake,elevlake,resht
     real, dimension(NLAKES),                  intent(in) :: qlakei,qlakeo  !-- inflow and outflow of lake
     real,                                     intent(in) :: dtrt

     character(len=*),                         intent(in) :: startdate
     character(len=*),                         intent(in) :: date
     integer lake_index(nlakes)

     call write_lake_real(latlake,lake_index,nlakes)
     call write_lake_real(lonlake,lake_index,nlakes)
     call write_lake_real(elevlake,lake_index,nlakes)
     call write_lake_real(resht,lake_index,nlakes)
     call write_lake_real(qlakei,lake_index,nlakes)
     call write_lake_real(qlakeo,lake_index,nlakes)
     if(my_id.eq. IO_id) then
        call output_lakes(igrid, split_output_count, NLAKES, &
           startdate, date, latlake, lonlake, elevlake, &
           qlakei,qlakeo, resht,dtrt,K)
     end if
     return
     end subroutine mpp_output_lakes

#endif

!----------------------------------- lake netcdf output
!-- output the ilake info an IDV 'station' compatible format
   subroutine output_lakes(igrid, split_output_count, NLAKES, &
        startdate, date, latlake, lonlake, elevlake, &
        qlakei,qlakeo, resht,dtrt,K)

!!output the routing variables over just channel
     integer,                                  intent(in) :: igrid, K
     integer,                                  intent(in) :: split_output_count
     integer,                                  intent(in) :: NLAKES
     real, dimension(NLAKES),                  intent(in) :: latlake,lonlake,elevlake,resht
     real, dimension(NLAKES),                  intent(in) :: qlakei,qlakeo  !-- inflow and outflow of lake
     real,                                     intent(in) :: dtrt

     character(len=*),                         intent(in) :: startdate
     character(len=*),                         intent(in) :: date

     integer, allocatable, DIMENSION(:)                   :: station_id
     integer, allocatable, DIMENSION(:)                   :: rec_num_of_lake

     integer, save  :: output_count
     integer, save :: ncid

     integer :: stationdim, dimdata, varid, charid, n
     integer :: iret,i, start_pos, prev_pos  !-- 
     integer :: previous_pos        !-- used for the station model
     character(len=256) :: output_flnm
     character(len=19)  :: date19, date19start
     character(len=34)  :: sec_since_date
     integer :: seconds_since,cnt
     character(len=32)  :: convention
     character(len=6),allocatable, DIMENSION(:)  :: stname

    !--- all this for writing the station id string
     INTEGER   TDIMS, TXLEN
     PARAMETER (TDIMS=2)    ! number of TX dimensions
     PARAMETER (TXLEN = 6) ! length of example string
     INTEGER  TIMEID        ! record dimension id
     INTEGER  TXID          ! variable ID
     INTEGER  TXDIMS(TDIMS) ! variable shape
     INTEGER  TSTART(TDIMS), TCOUNT(TDIMS)

!    sec_since_date = 'seconds since '//date(1:4)//'-'//date(6:7)//'-'//date(9:10)//' '//date(12:13)//':'//date(15:16)//' UTC'
!    seconds_since = int(dtrt)*output_count
     seconds_since = int(dtrt)*K

     allocate(station_id(NLAKES))
     allocate(rec_num_of_lake(NLAKES))
     allocate(stname(NLAKES))

     if (output_count == 0) then

!-- have moved sec_since_date from above here..
      sec_since_date = 'seconds since '//startdate(1:4)//'-'//startdate(6:7)//'-'//startdate(9:10) &
                  //' '//startdate(12:13)//':'//startdate(15:16)//' UTC'

      date19start(1:len_trim(startdate)) = startdate(1:4)//'-'//startdate(6:7)//'-'//startdate(9:10)//'_' &
                  //startdate(12:13)//':'//startdate(15:16)//':00'

      write(output_flnm, '(A12,".LAKEOUT_DOMAIN",I1)') date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16), igrid
#ifdef HYDRO_D
      print*, 'output_flnm = "'//trim(output_flnm)//'"'
#endif

#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
       iret = nf_create(trim(output_flnm), IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
#else
       iret = nf_create(trim(output_flnm), NF_CLOBBER, ncid)
#endif

#ifdef HYDRO_D
      if (iret /= 0) then
          print*, "Problem nf_create"
          call hydro_stop("output_lakes")
      endif
#endif

      do i=1,NLAKES
         station_id(i) = i
         write(stname(i),'(I6)') i
      enddo 

      iret = nf_def_dim(ncid, "recNum", NF_UNLIMITED, dimdata)  !--for linked list approach
      iret = nf_def_dim(ncid, "station", nlakes, stationdim)

      !- station location definition,  lat
      iret = nf_def_var(ncid,"latitude",NF_FLOAT, 1, (/stationdim/), varid)
      iret = nf_put_att_text(ncid,varid,'long_name',13,'Lake latitude')
      iret = nf_put_att_text(ncid,varid,'units',13,'degrees_north')

      !- station location definition,  long
      iret = nf_def_var(ncid,"longitude",NF_FLOAT, 1, (/stationdim/), varid)
      iret = nf_put_att_text(ncid,varid,'long_name',14,'Lake longitude')
      iret = nf_put_att_text(ncid,varid,'units',12,'degrees_east')

!     !-- lake's phyical elevation
      iret = nf_def_var(ncid,"altitude",NF_FLOAT, 1, (/stationdim/), varid)
      iret = nf_put_att_text(ncid,varid,'long_name',13,'Lake altitude')
      iret = nf_put_att_text(ncid,varid,'units',6,'meters')

     !-- parent index
      iret = nf_def_var(ncid,"parent_index",NF_INT,1,(/dimdata/), varid)
      iret = nf_put_att_text(ncid,varid,'long_name',33,'index of the lake for this record')

     !-- prevChild
      iret = nf_def_var(ncid,"prevChild",NF_INT,1,(/dimdata/), varid)
      iret = nf_put_att_text(ncid,varid,'long_name',54,'record number of the previous record for the same lake')
!ywtmp      iret = nf_put_att_int(ncid,varid,'_FillValue',NF_INT,2,-1)
      iret = nf_put_att_int(ncid,varid,'_FillValue',2,-1)

     !-- lastChild
      iret = nf_def_var(ncid,"lastChild",NF_INT,1,(/stationdim/), varid)
      iret = nf_put_att_text(ncid,varid,'long_name',27,'latest report for this lake')
!ywtmp      iret = nf_put_att_int(ncid,varid,'_FillValue',NF_INT,2,-1)
      iret = nf_put_att_int(ncid,varid,'_FillValue',2,-1)

!     !- water surface elevation
      iret = nf_def_var(ncid, "elevation", NF_FLOAT, 1, (/dimdata/), varid)
      iret = nf_put_att_text(ncid,varid,'units',6,'meters')
      iret = nf_put_att_text(ncid,varid,'long_name',14,'Lake Elevation')

!     !- inflow to lake
      iret = nf_def_var(ncid, "inflow", NF_FLOAT, 1, (/dimdata/), varid)
      iret = nf_put_att_text(ncid,varid,'units',13,'meter^3 / sec')

!     !- outflow to lake
      iret = nf_def_var(ncid, "outflow", NF_FLOAT, 1, (/dimdata/), varid)
      iret = nf_put_att_text(ncid,varid,'units',13,'meter^3 / sec')

     !-- station  id
     ! define character-position dimension for strings of max length 6
         iret = NF_DEF_DIM(ncid, "id_len", 6, charid)
         TXDIMS(1) = charid   ! define char-string variable and position dimension first
         TXDIMS(2) = stationdim
         iret = nf_def_var(ncid,"station_id",NF_CHAR, TDIMS, TXDIMS, varid)
         iret = nf_put_att_text(ncid,varid,'long_name',10,'Station id')

!     !- time definition, timeObs
         iret = nf_def_var(ncid,"time_observation",NF_INT, 1, (/dimdata/), varid)
         iret = nf_put_att_text(ncid,varid,'units',34,sec_since_date)
         iret = nf_put_att_text(ncid,varid,'long_name',19,'time of observation')

!       date19(1:19) = "0000-00-00_00:00:00"
!       date19(1:len_trim(startdate)) = startdate
!       iret = nf_put_att_text(ncid, NF_GLOBAL, "Conventions",32, convention)
!
        date19(1:19) = "0000-00-00_00:00:00"
        date19(1:len_trim(startdate)) = startdate
        convention(1:32) = "Unidata Observation Dataset v1.0"
        iret = nf_put_att_text(ncid, NF_GLOBAL, "Conventions",32, convention)
        iret = nf_put_att_text(ncid, NF_GLOBAL, "cdm_datatype",7, "Station")
        iret = nf_put_att_text(ncid, NF_GLOBAL, "geospatial_lat_max",4, "90.0")
        iret = nf_put_att_text(ncid, NF_GLOBAL, "geospatial_lat_min",5, "-90.0")
        iret = nf_put_att_text(ncid, NF_GLOBAL, "geospatial_lon_max",5, "180.0")
        iret = nf_put_att_text(ncid, NF_GLOBAL, "geospatial_lon_min",6, "-180.0")
        iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_start",19, date19start)
        iret = nf_put_att_text(ncid, NF_GLOBAL, "stationDimension",7, "station")
!!       iret = nf_put_att_text(ncid, NF_GLOBAL, "observationDimension",6, "recNum")
!!        iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coordinate",16,"time_observation")
        iret = nf_put_att_real(ncid, NF_GLOBAL, "missing_value", NF_FLOAT, 1, -9E15)
        iret = nf_enddef(ncid)

        !-- write latitudes
        iret = nf_inq_varid(ncid,"latitude", varid)
        iret = nf_put_vara_real(ncid, varid, (/1/), (/NLAKES/), LATLAKE)

        !-- write longitudes
        iret = nf_inq_varid(ncid,"longitude", varid)
        iret = nf_put_vara_real(ncid, varid, (/1/), (/NLAKES/), LONLAKE)

        !-- write physical height of lake
        iret = nf_inq_varid(ncid,"altitude", varid)
        iret = nf_put_vara_real(ncid, varid, (/1/), (/NLAKES/), elevlake)

        !-- write station id's 
         do i=1,nlakes
          TSTART(1) = 1
          TSTART(2) = i
          TCOUNT(1) = TXLEN
          TCOUNT(2) = 1
          iret = nf_inq_varid(ncid,"station_id", varid)
          iret = nf_put_vara_text(ncid, varid, TSTART, TCOUNT, stname(i))
         enddo

     endif

     output_count = output_count + 1

     cnt=0
     do i=1,NLAKES

         start_pos = (cnt+1)+(nlakes*(output_count-1))

         !!--time in seconds since startdate
         iret = nf_inq_varid(ncid,"time_observation", varid)
         iret = nf_put_vara_int(ncid, varid, (/start_pos/), (/1/), seconds_since)

         iret = nf_inq_varid(ncid,"elevation", varid)
         iret = nf_put_vara_real(ncid, varid, (/start_pos/), (/1/), resht(i))

         iret = nf_inq_varid(ncid,"inflow", varid)
         iret = nf_put_vara_real(ncid, varid, (/start_pos/), (/1/), qlakei(i))

         iret = nf_inq_varid(ncid,"outflow", varid)
         iret = nf_put_vara_real(ncid, varid, (/start_pos/), (/1/), qlakeo(i))

         !-- station index.. will repeat for every timesstep
         iret = nf_inq_varid(ncid,"parent_index", varid)
         iret = nf_put_vara_int(ncid, varid, (/start_pos/), (/1/), cnt)

          !--record number of previous record for same station
         prev_pos = cnt+(nlakes*(output_count-1))
         if(output_count.ne.1) then !-- only write next set of records
           iret = nf_inq_varid(ncid,"prevChild", varid)
           iret = nf_put_vara_int(ncid, varid, (/start_pos/), (/1/), prev_pos)
         endif

         cnt=cnt+1  !--indices are 0 based
         rec_num_of_lake(cnt) = start_pos-1  !-- save position for last child, 0-based!!

    enddo

      !-- lastChild variable gives the record number of the most recent report for the station
      iret = nf_inq_varid(ncid,"lastChild", varid)
      iret = nf_put_vara_int(ncid, varid, (/1/), (/nlakes/), rec_num_of_lake)

     !-- number of children reported for this station, OPTIONAL
     !--  iret = nf_inq_varid(ncid,"numChildren", varid)
     !--  iret = nf_put_vara_int(ncid, varid, (/1/), (/nlakes/), rec_num_of_lake)

    iret = nf_redef(ncid)
    date19(1:19) = "0000-00-00_00:00:00"
    date19(1:len_trim(date)) = date
    iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_end", 19, date19)
    iret = nf_enddef(ncid)

    iret = nf_sync(ncid)
     if (output_count == split_output_count) then
        output_count = 0
        iret = nf_close(ncid)
     endif

     deallocate(station_id)
     deallocate(rec_num_of_lake)
     deallocate(stname)
#ifdef HYDRO_D
     print *, "Exited Subroutine output_lakes"
#endif
     close(16)

 end subroutine output_lakes
!----------------------------------- lake netcdf output

#ifdef MPP_LAND

!-- output the channel route in an IDV 'grid' compatible format
   subroutine mpp_output_chrtgrd(igrid, split_output_count, ixrt,jxrt, &
        NLINKS,CH_NETLNK_in, startdate, date, &
        qlink, dt, geo_finegrid_flnm, gnlinks,map_l2g,g_ixrt,g_jxrt )

   USE module_mpp_land

     implicit none
#include <netcdf.inc>
     integer g_ixrt,g_jxrt
     integer,                                  intent(in) :: igrid
     integer,                                  intent(in) :: split_output_count
     integer,                                  intent(in) :: NLINKS,ixrt,jxrt
     real,                                     intent(in) :: dt
     real, dimension(NLINKS,2),                intent(in) :: qlink
     integer, dimension(IXRT,JXRT),            intent(in) :: CH_NETLNK_in
     character(len=*),          intent(in)     :: geo_finegrid_flnm
     character(len=*),                         intent(in) :: startdate
     character(len=*),                         intent(in) :: date

     integer::  gnlinks , map_l2g(nlinks)

     integer, allocatable,dimension(:,:)         :: CH_NETLNK
     real, allocatable,dimension(:,:)                :: g_qlink

     if(my_id .eq. io_id) then
        allocate(CH_NETLNK(g_IXRT,g_JXRT))
        allocate(g_qlink(gNLINKS,2) )
     else
        allocate(CH_NETLNK(1,1))
        allocate(g_qlink(1,2) )
     endif

     call write_chanel_real(qlink(:,1),map_l2g,gnlinks,nlinks,g_qlink(:,1))
     call write_chanel_real(qlink(:,2),map_l2g,gnlinks,nlinks,g_qlink(:,2))

     call write_IO_rt_int(CH_NETLNK_in, CH_NETLNK)

    if(my_id.eq.IO_id) then
        call  output_chrtgrd(igrid, split_output_count, g_ixrt,g_jxrt, &
           GNLINKS, CH_NETLNK, startdate, date, &
           g_qlink, dt, geo_finegrid_flnm)
    endif

     deallocate(g_qlink, CH_NETLNK)  
     return
     end subroutine mpp_output_chrtgrd
#endif

!-- output the channel route in an IDV 'grid' compatible format
   subroutine output_chrtgrd(igrid, split_output_count, ixrt,jxrt, &
        NLINKS, CH_NETLNK, startdate, date, &
        qlink, dt, geo_finegrid_flnm)

     integer,                                  intent(in) :: igrid
     integer,                                  intent(in) :: split_output_count
     integer,                                  intent(in) :: NLINKS,ixrt,jxrt
     real,                                     intent(in) :: dt
     real, dimension(NLINKS,2),                intent(in) :: qlink
     integer, dimension(IXRT,JXRT),            intent(in) :: CH_NETLNK
     character(len=*),          intent(in)     :: geo_finegrid_flnm
     character(len=*),                         intent(in) :: startdate
     character(len=*),                         intent(in) :: date
     character(len=32)  :: convention
     integer,save  :: output_count
     integer, save :: ncid,ncstatic
     real, dimension(IXRT,JXRT)          :: tmpflow
     real, dimension(IXRT)            :: xcoord
     real, dimension(JXRT)            :: ycoord
     real                                :: long_cm,lat_po,fe,fn
     real, dimension(2)                  :: sp

    integer :: varid, n
    integer :: jxlatdim,ixlondim,timedim !-- dimension ids

    integer :: iret,i,j
    character(len=256) :: output_flnm
    character(len=19)  :: date19
    character(len=34)  :: sec_since_date
 

    integer :: seconds_since




      tmpflow = -9E15

 
        write(output_flnm, '(A12,".CHRTOUT_GRID",I1)') date(1:4)//date(6:7)//date(9:10)//date(12:13)//date(15:16), igrid
#ifdef HYDRO_D
        print*, 'output_flnm = "'//trim(output_flnm)//'"'
#endif
 

!--- define dimension
#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
       iret = nf_create(trim(output_flnm), IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
#else
       iret = nf_create(trim(output_flnm), NF_CLOBBER, ncid)
#endif

#ifdef HYDRO_D
        if (iret /= 0) then
            print*, "Problem nf_create"
            call hydro_stop("output_chrtgrd")
        endif
#endif


        iret = nf_def_dim(ncid, "time", NF_UNLIMITED, timedim)
        iret = nf_def_dim(ncid, "x", ixrt, ixlondim)
        iret = nf_def_dim(ncid, "y", jxrt, jxlatdim)

!--- define variables
!     !- time definition, timeObs
         iret = nf_def_var(ncid,"time",NF_INT, 1, (/timedim/), varid)

       !- x-coordinate in cartesian system
!yw         iret = nf_def_var(ncid,"x",NF_DOUBLE, 1, (/ixlondim/), varid)
!yw         iret = nf_put_att_text(ncid,varid,'long_name',26,'x coordinate of projection')
!yw         iret = nf_put_att_text(ncid,varid,'standard_name',23,'projection_x_coordinate')
!yw         iret = nf_put_att_text(ncid,varid,'units',5,'Meter')

       !- y-coordinate in cartesian ssystem
!yw         iret = nf_def_var(ncid,"y",NF_DOUBLE, 1, (/jxlatdim/), varid)
!yw         iret = nf_put_att_text(ncid,varid,'long_name',26,'y coordinate of projection')
!yw         iret = nf_put_att_text(ncid,varid,'standard_name',23,'projection_y_coordinate')
!yw         iret = nf_put_att_text(ncid,varid,'units',5,'Meter')

!     !- flow definition, var
        iret = nf_def_var(ncid,"streamflow",NF_REAL, 3, (/ixlondim,jxlatdim,timedim/), varid)
        iret = nf_put_att_text(ncid,varid,'units',6,'m3 s-1')
        iret = nf_put_att_text(ncid,varid,'long_name',15,'water flow rate')
        iret = nf_put_att_text(ncid,varid,'coordinates',3,'x y')
        iret = nf_put_att_text(ncid,varid,'grid_mapping',23,'lambert_conformal_conic')
        iret = nf_put_att_real(ncid,varid,'missing_value',NF_REAL,1,-9E15)
        iret = nf_def_var(ncid,"index",NF_INT, 2, (/ixlondim,jxlatdim/), varid)


!-- place prjection information


      date19(1:19) = "0000-00-00_00:00:00"
      date19(1:len_trim(startdate)) = startdate
      convention(1:32) = "CF-1.0"
      iret = nf_put_att_real(ncid, NF_GLOBAL, "missing_value", NF_FLOAT, 1, -9E15)
      iret = nf_put_att_text(ncid, NF_GLOBAL, "Conventions",6, convention)
      iret = nf_put_att_text(ncid, NF_GLOBAL, "time_coverage_start",19, date19)

      iret = nf_enddef(ncid)

!!-- write latitude and longitude locations

!DJG inv    do j=jxrt,1,-1
    do j=1,jxrt
     do i=1,ixrt
       if(CH_NETLNK(i,j).GE.0) then
         tmpflow(i,j) = qlink(CH_NETLNK(i,j),1) 
       else
         tmpflow(i,j) = -9E15
       endif
     enddo
    enddo

!!time in seconds since startdate
    iret = nf_inq_varid(ncid,"index", varid)
    iret = nf_put_vara_int(ncid, varid, (/1,1/), (/ixrt,jxrt/),CH_NETLNK)

    iret = nf_inq_varid(ncid,"streamflow", varid)
    iret = nf_put_vara_real(ncid, varid, (/1,1,1/), (/ixrt,jxrt,1/),tmpflow)

        iret = nf_close(ncid)



 end subroutine output_chrtgrd


 subroutine read_chan_forcing( &
       indir,olddate,startdate,hgrid,&
       ixrt,jxrt,QSTRMVOLRT_ACC,QINFLOWBASE,QSUBRT)
! This subrouting is going to read channel forcing for
!  channel only simulations (ie when CHANRTSWCRT = 2)

   implicit none
#include <netcdf.inc>
   ! in variable
   character(len=*) :: olddate,hgrid,indir,startdate
   character(len=256) :: filename
   integer :: ixrt,jxrt
   real,dimension(ixrt,jxrt):: QSTRMVOLRT_ACC,QINFLOWBASE,QSUBRT
   ! tmp variable
   character(len=256) :: inflnm, product
   integer  :: i,j,mmflag
   character(len=256) :: units
   integer :: ierr
   integer :: ncid


!DJG Create filename...
        inflnm = trim(indir)//"/"//&
             olddate(1:4)//olddate(6:7)//olddate(9:10)//olddate(12:13)//&
             olddate(15:16)//".RTOUT_DOMAIN"//hgrid
#ifdef HYDRO_D
        print *, "Channel forcing file...",inflnm
#endif


!DJG Open NetCDF file...
#ifdef PARALLELIO
    iret = nf90_open_par(trim(inflnm), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#else
    ierr = nf_open(inflnm, NF_NOWRITE, ncid)
#endif
    if (ierr /= 0) then
#ifdef HYDRO_D
       write(*,'("READFORC_chan Problem opening netcdf file: ''", A, "''")') trim(inflnm)
       call hydro_stop("read_chan_forcing")
#endif
    endif

!DJG read data...
    call get_2d_netcdf("QSTRMVOLRT",  ncid, QSTRMVOLRT_ACC, units, ixrt, jxrt, .TRUE., ierr)
!DJG TBC    call get_2d_netcdf("T2D", ncid, t,     units, ixrt, jxrt, .TRUE., ierr)
!DJG TBC    call get_2d_netcdf("T2D", ncid, t,     units, ixrt, jxrt, .TRUE., ierr)

    ierr = nf_close(ncid)

 end subroutine read_chan_forcing




      subroutine get2d_int(var_name,out_buff,ix,jx,fileName)
          implicit none
#include <netcdf.inc>
          integer :: iret,varid,ncid,ix,jx
          integer out_buff(ix,jx)
          character(len=*), intent(in) :: var_name
          character(len=*), intent(in) :: fileName
#ifdef PARALLELIO
    iret = nf90_open_par(trim(fileName), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
#else
          iret = nf_open(trim(fileName), NF_NOWRITE, ncid)
#endif
          if (iret .ne. 0) then
#ifdef HYDRO_D
            print*,"aaa failed to open the netcdf file: ",trim(fileName)
            call hydro_stop("get2d_int")
#endif
          endif
          iret = nf_inq_varid(ncid,trim(var_name),  varid)
          if(iret .ne. 0) then
#ifdef HYDRO_D
            print*,"failed to read the variabe: ",trim(var_name)
            print*,"failed to read the netcdf file: ",trim(fileName)
#endif
          endif
          iret = nf_get_var_int(ncid, varid, out_buff)
          iret = nf_close(ncid)
         return
      end subroutine get2d_int


#ifdef MPP_LAND
       SUBROUTINE MPP_READ_ROUTEDIM(did,g_IXRT,g_JXRT, GCH_NETLNK,GNLINKS,IXRT,JXRT, &
            route_chan_f,route_link_f, &
            route_direction_f, route_lake_f,NLINKS, NLAKES, &
            CH_NETLNK, channel_option, geo_finegrid_flnm)


         USE module_mpp_land

         implicit none
#include <netcdf.inc>
        INTEGER                                      :: channel_option, did
        INTEGER                                      :: g_IXRT,g_JXRT
        INTEGER, INTENT(INOUT)                       :: NLINKS, NLAKES, GNLINKS
        INTEGER, INTENT(IN)                          :: IXRT,JXRT
        INTEGER                                      :: CHNID,cnt
        INTEGER, DIMENSION(IXRT,JXRT)                :: CH_NETRT   !- binary channel mask
        INTEGER, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: CH_NETLNK  !- each node gets unique id
        INTEGER, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: GCH_NETLNK  !- each node gets unique id based on global domain
        ! INTEGER, DIMENSION(g_IXRT,g_JXRT) :: g_CH_NETLNK  ! temp array
        INTEGER, allocatable,DIMENSION(:,:) :: g_CH_NETLNK  ! temp array
        INTEGER, DIMENSION(IXRT,JXRT)                :: DIRECTION  !- flow direction
        INTEGER, DIMENSION(IXRT,JXRT)                :: LAKE_MSKRT
        REAL, DIMENSION(IXRT,JXRT)                   :: LAT, LON
        integer:: i,j

        CHARACTER(len=256)       :: route_chan_f, route_link_f,route_direction_f,route_lake_f
        CHARACTER(len=256)       :: geo_finegrid_flnm
!       CHARACTER(len=*)       :: geo_finegrid_flnm

!       integer, allocatable, dimension(:) :: tmp_int
        integer :: ywcount

        if(my_id .eq. IO_id) then
           allocate(g_CH_NETLNK(g_IXRT,g_JXRT))
           g_CH_NETLNK = -9999
           CALL READ_ROUTEDIM(g_IXRT, g_JXRT, route_chan_f, route_link_f, &
              route_direction_f, route_lake_f, GNLINKS, NLAKES, &
              g_CH_NETLNK, channel_option,geo_finegrid_flnm)
        else
           allocate(g_CH_NETLNK(1,1))
        endif
        

        call mpp_land_bcast_int1(NLAKES)
        call mpp_land_bcast_int1(GNLINKS)


        call decompose_RT_int(g_CH_NETLNK,GCH_NETLNK,g_IXRT,g_JXRT,ixrt,jxrt)
        deallocate(g_CH_NETLNK)
        ywcount = 0 
        CH_NETLNK = -9999
        do j = 1, jxrt
           do i = 1, ixrt
                  if(GCH_NETLNK(i,j) .gt. 0) then
                       ywcount = ywcount + 1
                       CH_NETLNK(i,j) = ywcount
                  endif
           end do
        end do
        NLINKS = ywcount

!ywcheck
!        CH_NETLNK = GCH_NETLNK


        allocate(rt_domain(did)%map_l2g(NLINKS))

        rt_domain(did)%map_l2g = -1
        do j = 1, jxrt
           do i = 1, ixrt
              if(CH_NETLNK(i,j) .gt. 0) then
                  rt_domain(did)%map_l2g(CH_NETLNK(i,j)) = GCH_NETLNK(i,j)
              endif
           end do       
        end do       

        call mpp_chrt_nlinks_collect(NLINKS)
        return 
        end SUBROUTINE MPP_READ_ROUTEDIM




#endif
        
      SUBROUTINE READ_ROUTING_seq(IXRT,JXRT,ELRT,CH_NETRT,LKSATFAC,route_topo_f,    &
            route_chan_f, geo_finegrid_flnm,OVROUGHRTFAC,RETDEPRTFAC)


#include <netcdf.inc>
!        INTEGER, INTENT(IN)                          :: IXRT,JXRT
!        REAL,    INTENT(INOUT), DIMENSION(IXRT,JXRT) :: ELRT,LKSATFAC
!        INTEGER, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: CH_NETRT
!Dummy inverted grids
!        REAL,    INTENT(INOUT), DIMENSION(IXRT,JXRT) :: OVROUGHRTFAC
!        REAL,    INTENT(INOUT), DIMENSION(IXRT,JXRT) :: RETDEPRTFAC

        ! modified by Qian, 2016.06
        INTEGER,           INTENT(IN)                          :: IXRT, JXRT
        REAL,              INTENT(INOUT), DIMENSION(IXRT,JXRT) :: ELRT
        INTEGER, OPTIONAL, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: CH_NETRT
        REAL,    OPTIONAL, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: LKSATFAC
        CHARACTER(len=256), OPTIONAL, INTENT(IN)               :: route_topo_f
        CHARACTER(len=256), OPTIONAL, INTENT(IN)               :: route_chan_f
        CHARACTER(len=256),           INTENT(IN)               :: geo_finegrid_flnm
        REAL,    OPTIONAL, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: OVROUGHRTFAC
        REAL,    OPTIONAL, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: RETDEPRTFAC

        INTEGER                  :: I,J, iret, jj
        CHARACTER(len=256)       :: var_name

        var_name = "TOPOGRAPHY"
        call nreadRT2d_real(var_name,ELRT,ixrt,jxrt,&
                     trim(geo_finegrid_flnm))
#ifdef HYDRO_D
        write(6,*) "read ",var_name
#endif

!!!DY to be fixed ... 6/27/08
!        var_name = "BED_ELEVATION"
!        iret = get2d_real(var_name,ELRT,ixrt,jxrt,&
!                     trim(geo_finegrid_flnm))

        IF (PRESENT(CH_NETRT)) THEN
           var_name = "CHANNELGRID"
           call nreadRT2d_int(var_name,CH_NETRT,ixrt,jxrt,&
                              trim(geo_finegrid_flnm))
#ifdef HYDRO_D
           write(6,*) "read ",var_name
#endif
        ENDIF

        IF (PRESENT(LKSATFAC)) THEN
           var_name = "LKSATFAC"
           LKSATFAC = -9999.9
           call nreadRT2d_real(var_name,LKSATFAC,ixrt,jxrt,&
                               trim(geo_finegrid_flnm))
#ifdef HYDRO_D
           write(6,*) "read ",var_name
#endif
           where (LKSATFAC == -9999.9) LKSATFAC = 1000.0  !specify LKSAFAC if no term avail...
        ENDIF

!1.12.2012...Read in routing calibration factors...
        IF (PRESENT(RETDEPRTFAC)) THEN
           var_name = "RETDEPRTFAC"
           call nreadRT2d_real(var_name,RETDEPRTFAC,ixrt,jxrt,&
                               trim(geo_finegrid_flnm))
#ifdef HYDRO_D
           write(6,*) "read ", var_name
#endif
           where (RETDEPRTFAC < 0.) RETDEPRTFAC = 1.0  ! reset grid to = 1.0 if non-valid value exists
        ENDIF

     
        IF (PRESENT(OVROUGHRTFAC)) THEN
           var_name = "OVROUGHRTFAC"
           call nreadRT2d_real(var_name,OVROUGHRTFAC,ixrt,jxrt,&
                               trim(geo_finegrid_flnm))
           where (OVROUGHRTFAC <= 0.) OVROUGHRTFAC = 1.0 ! reset grid to = 1.0 if non-valid value exists
#ifdef HYDRO_D
           write(6,*) "read ", var_name
#endif
        ENDIF

#ifdef HYDRO_D
        write(6,*) "finish READ_ROUTING_seq"
#endif
 
        return

!DJG -----------------------------------------------------
   END SUBROUTINE READ_ROUTING_seq

!DJG _____________________________
   subroutine output_lsm(outFile,did)


   implicit none

   integer did

   character(len=*) outFile

    integer :: ncid,irt, dimid_ix, dimid_jx,  &
             dimid_ixrt, dimid_jxrt, varid, &
             dimid_links, dimid_basns, dimid_soil
    integer :: iret, n
    character(len=2) tmpStr 



#ifdef MPP_LAND
     if(IO_id.eq.my_id) &
#endif

#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
       iret = nf_create(trim(outFile), IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
#else
       iret = nf_create(trim(outFile), NF_CLOBBER, ncid)
#endif

#ifdef MPP_LAND
       call mpp_land_bcast_int1(iret)
#endif

#ifdef HYDRO_D
       if (iret /= 0) then
          print*, "Problem nf_create"
          call hydro_stop("output_lsm")
       endif
#endif


#ifdef MPP_LAND
     if(IO_id.eq.my_id) then
#endif
#ifdef HYDRO_D
         write(6,*) "output file ", outFile
#endif
! define dimension for variables 
          iret = nf_def_dim(ncid, "depth", nlst_rt(did)%nsoil, dimid_soil)  !-- 3-d soils
   
#ifdef MPP_LAND
          iret = nf_def_dim(ncid, "ix", global_nx, dimid_ix)  !-- make a decimated grid
          iret = nf_def_dim(ncid, "iy", global_ny, dimid_jx)
#else
          iret = nf_def_dim(ncid, "ix", rt_domain(did)%ix, dimid_ix)  !-- make a decimated grid
          iret = nf_def_dim(ncid, "iy", rt_domain(did)%jx, dimid_jx)
#endif
    
!define variables
          do n = 1, nlst_rt(did)%nsoil
             if( n .lt. 10) then
                write(tmpStr, '(i1)') n
             else
                write(tmpStr, '(i2)') n
             endif
             iret = nf_def_var(ncid,"stc"//trim(tmpStr),NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
             iret = nf_def_var(ncid,"smc"//trim(tmpStr),NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
             iret = nf_def_var(ncid,"sh2ox"//trim(tmpStr),NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          end do

          iret = nf_def_var(ncid,"smcmax1",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          iret = nf_def_var(ncid,"smcref1",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          iret = nf_def_var(ncid,"smcwlt1",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          iret = nf_def_var(ncid,"infxsrt",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          iret = nf_def_var(ncid,"sfcheadrt",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)

          iret = nf_enddef(ncid)

#ifdef MPP_LAND
    endif
#endif
        call w_rst_nc3(ncid,rt_domain(did)%ix,rt_domain(did)%jx,nlst_rt(did)%nsoil,rt_domain(did)%stc,"stc")
        call w_rst_nc3(ncid,rt_domain(did)%ix,rt_domain(did)%jx,nlst_rt(did)%nsoil,rt_domain(did)%smc,"smc")
        call w_rst_nc3(ncid,rt_domain(did)%ix,rt_domain(did)%jx,nlst_rt(did)%nsoil,rt_domain(did)%sh2ox,"sh2ox")
        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%SMCMAX1,"smcmax1") 
        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%SMCREF1,"smcref1" )
        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%SMCWLT1,"smcwlt1"  ) 
        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%INFXSRT,"infxsrt"  ) 
        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%SFCHEADRT,"sfcheadrt" )


#ifdef MPP_LAND
     if(IO_id.eq.my_id) then
#endif

        iret = nf_close(ncid)
#ifdef HYDRO_D
        write(6,*) "finish writing outFile : ", outFile
#endif

#ifdef MPP_LAND
    endif
#endif

        return
        end subroutine output_lsm


   subroutine RESTART_OUT_nc(outFile,did)
   implicit none

   integer :: did
   integer :: n
   character(len=2) :: tmpStr
   character(len=*) :: outFile

   integer :: ncid,irt,    dimid_ix,    dimid_jx,   &
              dimid_ixrt,  dimid_jxrt,  varid,      &
              dimid_links, dimid_basns, dimid_soil, &
              dimid_lakes
   integer :: iret


#ifdef MPP_LAND
   if(IO_id.eq.my_id) &
#endif

#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
     iret = nf_create(trim(outFile), IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
     write(6,*) "yyywww using large netcdf file definition. "
     call flush(6)
#else
     iret = nf_create(trim(outFile), NF_CLOBBER, ncid)
     write(6,*) "yyywww do not use large netcdf file definition. "
     call flush(6)
#endif


#ifdef MPP_LAND
     call mpp_land_bcast_int1(iret)
#endif

#ifdef HYDRO_D
     if (iret /= 0) then
        print*, "Problem nf_create"
        call hydro_stop("RESTART_OUT_nc")
     endif
#endif

#ifdef MPP_LAND
     if(IO_id.eq.my_id) then
#endif
       ! define dimension for variables 
       iret = nf_def_dim(ncid, "depth", nlst_rt(did)%nsoil, dimid_soil)  !-- 3-d soils
   
#ifdef MPP_LAND
       iret = nf_def_dim(ncid, "ix", global_nx, dimid_ix)                !-- make a decimated grid
       iret = nf_def_dim(ncid, "iy", global_ny, dimid_jx)
       iret = nf_def_dim(ncid, "ixrt", global_rt_nx , dimid_ixrt)        !-- make a decimated grid
       iret = nf_def_dim(ncid, "iyrt", global_rt_ny, dimid_jxrt)
#else
       iret = nf_def_dim(ncid, "ix", rt_domain(did)%ix, dimid_ix)        !-- make a decimated grid
       iret = nf_def_dim(ncid, "iy", rt_domain(did)%jx, dimid_jx)
       iret = nf_def_dim(ncid, "ixrt", rt_domain(did)%ixrt , dimid_ixrt) !-- make a decimated grid
       iret = nf_def_dim(ncid, "iyrt", rt_domain(did)%jxrt, dimid_jxrt)
#endif

       if(nlst_rt(did)%SUBRTSWCRT .EQ. 1 .OR. nlst_rt(did)%OVRTSWCRT .EQ.1 ) then  ! Qian, 2016.08
          iret = nf_def_dim(ncid, "links", rt_domain(did)%gnlinks, dimid_links)
          if(rt_domain(did)%nlakes .gt. 0) then
            iret = nf_def_dim(ncid, "lakes", rt_domain(did)%nlakes, dimid_lakes)
          endif
            iret = nf_def_dim(ncid, "basns", rt_domain(did)%numbasns, dimid_basns)
          endif

          ! define variables
          do n = 1, nlst_rt(did)%nsoil
             if( n .lt. 10) then
                write(tmpStr, '(i1)') n
             else
                write(tmpStr, '(i2)') n
             endif
             iret = nf_def_var(ncid,"stc"//trim(tmpStr),NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
             iret = nf_def_var(ncid,"smc"//trim(tmpStr),NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
             iret = nf_def_var(ncid,"sh2ox"//trim(tmpStr),NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          end do
    
          iret = nf_def_var(ncid,"smcmax1",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          iret = nf_def_var(ncid,"smcref1",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          iret = nf_def_var(ncid,"smcwlt1",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          iret = nf_def_var(ncid,"infxsrt",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          iret = nf_def_var(ncid,"soldrain",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)
          iret = nf_def_var(ncid,"sfcheadrt",NF_FLOAT,2,(/dimid_ix,dimid_jx/),varid)

          if ( nlst_rt(did)%SUBRTSWCRT .EQ.1 .OR. nlst_rt(did)%OVRTSWCRT  .EQ.1 .OR. &
               nlst_rt(did)%GWBASESWCRT.NE.0 .OR. nlst_rt(did)%OVCHRTSWCRT.EQ.1      &
             ) then

            if (nlst_rt(did)%SUBRTSWCRT.EQ.1 .OR. nlst_rt(did)%OVRTSWCRT.EQ.1) then  ! Qian, 2016.08 
               iret = nf_def_var(ncid,"QBDRYRT",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
               iret = nf_def_var(ncid,"qstrmvolrt",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
            endif

            iret = nf_def_var(ncid,"infxswgt",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
            iret = nf_def_var(ncid,"SFCHEADSUBRT",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)

            do n = 1, nlst_rt(did)%nsoil
               if ( n .lt. 10) then
                  write(tmpStr, '(i1)') n
               else
                  write(tmpStr, '(i2)') n
               endif
               iret = nf_def_var(ncid,"sh2owgt"//trim(tmpStr),NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
            end do
            iret = nf_def_var(ncid,"RETDEPRT",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)

            if(nlst_rt(did)%CHANRTSWCRT.EQ.1) then
              iret = nf_def_var(ncid,"hlink",NF_FLOAT,1,(/dimid_links/),varid)
              iret = nf_def_var(ncid,"qlink1",NF_FLOAT,1,(/dimid_links/),varid)
              iret = nf_def_var(ncid,"qlink2",NF_FLOAT,1,(/dimid_links/),varid)
              iret = nf_def_var(ncid,"cvol",NF_FLOAT,1,(/dimid_links/),varid)
              if(rt_domain(did)%nlakes .gt. 0) then
                  iret = nf_def_var(ncid,"resht",NF_FLOAT,1,(/dimid_lakes/),varid)
                  iret = nf_def_var(ncid,"qlakeo",NF_FLOAT,1,(/dimid_lakes/),varid)
              endif
              iret = nf_def_var(ncid,"lake_inflort",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
            end if

            if(nlst_rt(did)%GWBASESWCRT.EQ.1) then
              iret = nf_def_var(ncid,"z_gwsubbas",NF_FLOAT,1,(/dimid_basns/),varid)
              !yw test bucket model
!             iret = nf_def_var(ncid,"gwbas_pix_ct",NF_FLOAT,1,(/dimid_basns/),varid)
!             iret = nf_def_var(ncid,"gw_buck_exp",NF_FLOAT,1,(/dimid_basns/),varid)
!             iret = nf_def_var(ncid,"z_max",NF_FLOAT,1,(/dimid_basns/),varid)
!             iret = nf_def_var(ncid,"gw_buck_coeff",NF_FLOAT,1,(/dimid_basns/),varid)
!             iret = nf_def_var(ncid,"qin_gwsubbas",NF_FLOAT,1,(/dimid_basns/),varid)
!             iret = nf_def_var(ncid,"qinflowbase",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
!             iret = nf_def_var(ncid,"qout_gwsubbas",NF_FLOAT,1,(/dimid_basns/),varid)
            end if

            if(nlst_rt(did)%gwBaseSwCRT .eq. 3)then
              iret = nf_def_var(ncid,"HEAD",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
            end if 
 
          end if   

          if ((nlst_rt(did)%OVCHRTSWCRT .eq. 1) .OR. (nlst_rt(did)%INTERGWSFW .eq. 1)) then
             iret = nf_def_var(ncid,"HLAKE",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
             iret = nf_def_var(ncid,"HLAKEAV",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
             iret = nf_def_var(ncid,"SINFILAV",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
             if (nlst_rt(did)%INTERGWSFW .eq. 1) &
             iret = nf_def_var(ncid,"HEND",NF_FLOAT,2,(/dimid_ixrt,dimid_jxrt/),varid)
          endif
	      
!         output global attribute
          iret = nf_put_att_int(ncid,NF_GLOBAL,"his_out_counts",NF_INT, 1,rt_domain(did)%his_out_counts)
          iret = nf_put_att_text(ncid,NF_GLOBAL,"Restart_Time",19,nlst_rt(did)%olddate(1:19))
          iret = nf_put_att_text(ncid,NF_GLOBAL,"Since_Date",19,nlst_rt(did)%sincedate(1:19))
          iret = nf_put_att_real(ncid,NF_GLOBAL,"DTCT",NF_REAL, 1,nlst_rt(did)%DTCT)
          iret = nf_enddef(ncid)

#ifdef MPP_LAND
    endif
#endif
        call w_rst_nc3(ncid,rt_domain(did)%ix,rt_domain(did)%jx,nlst_rt(did)%nsoil,rt_domain(did)%stc,"stc")
        call w_rst_nc3(ncid,rt_domain(did)%ix,rt_domain(did)%jx,nlst_rt(did)%nsoil,rt_domain(did)%smc,"smc")
        call w_rst_nc3(ncid,rt_domain(did)%ix,rt_domain(did)%jx,nlst_rt(did)%nsoil,rt_domain(did)%sh2ox,"sh2ox")


        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%SMCMAX1,"smcmax1") 
        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%SMCREF1,"smcref1" )
        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%SMCWLT1,"smcwlt1"  ) 
        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%INFXSRT,"infxsrt"  ) 
        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%soldrain,"soldrain"  ) 
        call w_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%sfcheadrt,"sfcheadrt"  ) 


ifvar:  if (nlst_rt(did)%SUBRTSWCRT  .EQ. 1 .OR. nlst_rt(did)%OVRTSWCRT   .EQ.1 .OR. &
            nlst_rt(did)%GWBASESWCRT .ne. 0 .OR. nlst_rt(did)%OVCHRTSWCRT .EQ.1      &
           ) then  

            call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%INFXSWGT, "infxswgt" )
            call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%SFCHEADSUBRT, "SFCHEADSUBRT" )
            call w_rst_rt_nc3(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,nlst_rt(did)%nsoil,rt_domain(did)%SH2OWGT, "sh2owgt" )
            call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%RETDEPRT, "RETDEPRT" )
       
            if (nlst_rt(did)%SUBRTSWCRT.EQ.1 .OR. nlst_rt(did)%OVRTSWCRT.EQ.1  ) then  ! Qian, 2016.08
               call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%QBDRYRT, "QBDRYRT" )
               call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%QSTRMVOLRT, "qstrmvolrt" )
            endif

            if(nlst_rt(did)%CHANRTSWCRT.EQ.1) then
              call w_rst_crt_nc1(ncid,rt_domain(did)%nlinks,rt_domain(did)%HLINK,"hlink" &
#ifdef MPP_LAND
                 ,rt_domain(did)%map_l2g, rt_domain(did)%gnlinks  &
#endif
                  )

              call w_rst_crt_nc1(ncid,rt_domain(did)%nlinks,rt_domain(did)%QLINK(:,1),"qlink1" &
#ifdef MPP_LAND
                 ,rt_domain(did)%map_l2g, rt_domain(did)%gnlinks  &
#endif
                  )


              call w_rst_crt_nc1(ncid,rt_domain(did)%nlinks,rt_domain(did)%QLINK(:,2),"qlink2" &
#ifdef MPP_LAND
                 ,rt_domain(did)%map_l2g, rt_domain(did)%gnlinks  &
#endif
                  )


              call w_rst_crt_nc1(ncid,rt_domain(did)%nlinks,rt_domain(did)%cvol,"cvol" &
#ifdef MPP_LAND
                 ,rt_domain(did)%map_l2g, rt_domain(did)%gnlinks  &
#endif
                  )

!              call w_rst_crt_nc1(ncid,rt_domain(did)%nlinks,rt_domain(did)%resht,"resht" &
!#ifdef MPP_LAND
!                 ,rt_domain(did)%map_l2g, rt_domain(did)%gnlinks  &
!#endif
!                  )

              call w_rst_crt_nc1_lake(ncid,rt_domain(did)%nlakes,rt_domain(did)%resht,"resht" &
#ifdef MPP_LAND
                 ,rt_domain(did)%lake_index  &
#endif
                  )

              call w_rst_crt_nc1_lake(ncid,rt_domain(did)%nlakes,rt_domain(did)%qlakeo,"qlakeo" &
#ifdef MPP_LAND
                 ,rt_domain(did)%lake_index  &
#endif
                  )

              call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%LAKE_INFLORT,"lake_inflort")

            end if  ! if(nlst_rt(did)%CHANRTSWCRT.EQ.1)

            if (nlst_rt(did)%GWBASESWCRT.EQ.1) then
              call w_rst_crt_nc1g(ncid,rt_domain(did)%numbasns,rt_domain(did)%z_gwsubbas,"z_gwsubbas" )
!yw test bucket model
!             call w_rst_crt_nc1g(ncid,rt_domain(did)%numbasns,rt_domain(did)%gwbas_pix_ct,"gwbas_pix_ct" )
!             call w_rst_crt_nc1g(ncid,rt_domain(did)%numbasns,rt_domain(did)%gw_buck_exp,"gw_buck_exp" )
!             call w_rst_crt_nc1g(ncid,rt_domain(did)%numbasns,rt_domain(did)%z_max,"z_max" )
!             call w_rst_crt_nc1g(ncid,rt_domain(did)%numbasns,rt_domain(did)%gw_buck_coeff,"gw_buck_coeff" )
!             call w_rst_crt_nc1g(ncid,rt_domain(did)%numbasns,rt_domain(did)%qin_gwsubbas,"qin_gwsubbas" )
!             call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%qinflowbase,"qinflowbase")
!             call w_rst_crt_nc1g(ncid,rt_domain(did)%numbasns,rt_domain(did)%qout_gwsubbas,"qout_gwsubbas" )
            end if  ! if (nlst_rt(did)%GWBASESWCRT.EQ.1)

            if(nlst_rt(did)%GWBASESWCRT.EQ.3) then
	      call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,gw2d(did)%h, "HEAD" )
            end if

        end if ifvar

        ! added by Qian, 2016.11
        if ((nlst_rt(did)%OVCHRTSWCRT .eq. 1) .OR. (nlst_rt(did)%INTERGWSFW .eq. 1)) then
           call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt2d(did)%hlake, "HLAKE" )
           call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt2d(did)%hlakeav, "HLAKEAV" )
           call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt2d(did)%sinfilav, "SINFILAV" )
           if (nlst_rt(did)%INTERGWSFW .eq. 1) &
           call w_rst_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt2d(did)%hend, "HEND")
        endif


#ifdef MPP_LAND
        if(IO_id.eq.my_id) &
#endif
        iret = nf_close(ncid)

        return
        end subroutine RESTART_OUT_nc

#ifdef MPP_LAND

   subroutine RESTART_OUT_bi(outFile,did)


   implicit none

   integer did

   character(len=*) outFile

    integer :: iunit


    iunit = 80
    open(iunit, file = outFile, form="unformatted")
          write(iunit) rt_domain(did)%his_out_counts
          write(iunit) nlst_rt(did)%olddate(1:19)
          write(iunit) nlst_rt(did)%sincedate(1:19)
          write(iunit) nlst_rt(did)%DTCT 
          write(iunit) rt_domain(did)%stc
          write(iunit) rt_domain(did)%smc
          write(iunit) rt_domain(did)%sh2ox
          write(iunit) rt_domain(did)%SMCMAX1
          write(iunit) rt_domain(did)%SMCREF1
          write(iunit) rt_domain(did)%SMCWLT1
          write(iunit) rt_domain(did)%INFXSRT
          write(iunit) rt_domain(did)%soldrain
          write(iunit) rt_domain(did)%sfcheadrt

          if(nlst_rt(did)%SUBRTSWCRT.EQ.1.OR.nlst_rt(did)%OVRTSWCRT.EQ.1 .or. nlst_rt(did)%GWBASESWCRT .ne. 0) then
                write(iunit) rt_domain(did)%QBDRYRT
                write(iunit) rt_domain(did)%INFXSWGT
                write(iunit) rt_domain(did)%SFCHEADSUBRT
                write(iunit) rt_domain(did)%SH2OWGT
                write(iunit) rt_domain(did)%QSTRMVOLRT
                write(iunit) rt_domain(did)%RETDEPRT
                if(nlst_rt(did)%CHANRTSWCRT.EQ.1) then
                   write(iunit) rt_domain(did)%HLINK
                   write(iunit) rt_domain(did)%QLINK(:,1)
                   write(iunit) rt_domain(did)%QLINK(:,2)
                   write(iunit) rt_domain(did)%cvol
                   write(iunit) rt_domain(did)%resht
                   write(iunit) rt_domain(did)%qlakeo
                   write(iunit) rt_domain(did)%LAKE_INFLORT
                end if
                if(nlst_rt(did)%GWBASESWCRT.EQ.1) then
                   write(iunit) rt_domain(did)%z_gwsubbas
                end if
                if(nlst_rt(did)%GWBASESWCRT .eq. 3) then
                   write(iunit) gw2d(did)%h
                endif
          end if  

          ! added by Qian, 2016.11
          if ((nlst_rt(did)%OVCHRTSWCRT .eq. 1) .OR. (nlst_rt(did)%INTERGWSFW .eq. 1)) then
             write(iunit) rt2d(did)%hlake
             write(iunit) rt2d(did)%hlakeav
             write(iunit) rt2d(did)%sinfilav
             if (nlst_rt(did)%INTERGWSFW .eq. 1) &
             write(iunit) rt2d(did)%hend
          endif

        close(iunit)

        return
        end subroutine RESTART_OUT_bi

   subroutine RESTART_in_bi(inFileTmp,did)


   implicit none

   integer did

   character(len=*) inFileTmp
   character(len=256) inFile
   character(len=19) str_tmp

    integer :: iunit
    logical :: fexist


    iunit = 80

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
             endif

    inFile = trim(inFileTmp)//"."//str_tmp 

    inquire (file=trim(inFile), exist=fexist)
    if(.not. fexist) then
        call hydro_stop("could not find restart file")
    endif

    open(iunit, file = inFile, form="unformatted")
          read(iunit) rt_domain(did)%his_out_counts
          read(iunit) nlst_rt(did)%olddate(1:19)
          read(iunit) nlst_rt(did)%sincedate(1:19)
          read(iunit) nlst_rt(did)%DTCT 
          read(iunit) rt_domain(did)%stc
          read(iunit) rt_domain(did)%smc
          read(iunit) rt_domain(did)%sh2ox
          read(iunit) rt_domain(did)%SMCMAX1
          read(iunit) rt_domain(did)%SMCREF1
          read(iunit) rt_domain(did)%SMCWLT1
          read(iunit) rt_domain(did)%INFXSRT
          read(iunit) rt_domain(did)%soldrain
          read(iunit) rt_domain(did)%sfcheadrt

          if(nlst_rt(did)%SUBRTSWCRT.EQ.1.OR.nlst_rt(did)%OVRTSWCRT.EQ.1 .or. nlst_rt(did)%GWBASESWCRT .ne. 0) then
                read(iunit) rt_domain(did)%QBDRYRT
                read(iunit) rt_domain(did)%INFXSWGT
                read(iunit) rt_domain(did)%SFCHEADSUBRT
                read(iunit) rt_domain(did)%SH2OWGT
                read(iunit) rt_domain(did)%QSTRMVOLRT
                read(iunit) rt_domain(did)%RETDEPRT
                if(nlst_rt(did)%CHANRTSWCRT.EQ.1) then
                   read(iunit) rt_domain(did)%HLINK
                   read(iunit) rt_domain(did)%QLINK(:,1)
                   read(iunit) rt_domain(did)%QLINK(:,2)
                   read(iunit) rt_domain(did)%cvol
                   read(iunit) rt_domain(did)%resht
                   read(iunit) rt_domain(did)%qlakeo
                   read(iunit) rt_domain(did)%LAKE_INFLORT
                end if
                if(nlst_rt(did)%GWBASESWCRT.EQ.1) then
                   read(iunit) rt_domain(did)%z_gwsubbas
                end if
                if(nlst_rt(did)%GWBASESWCRT .eq. 3) then
                   read(iunit) gw2d(did)%ho
                endif
          end if  

          ! added by Qian, 2016.11
          if ((nlst_rt(did)%OVCHRTSWCRT .eq. 1) .OR. (nlst_rt(did)%INTERGWSFW .eq. 1)) then
             read(iunit) rt2d(did)%hlake
             if (nlst_rt(did)%INTERGWSFW .eq. 1) &
             read(iunit) rt2d(did)%hbeg
          endif

        close(iunit)

        return
        end subroutine RESTART_in_bi
#endif

        subroutine w_rst_rt_nc2(ncid,ix,jx,inVar,varName)
           implicit none
           integer:: ncid,ix,jx,varid , iret
           character(len=*) varName
           real, dimension(ix,jx):: inVar
#ifdef MPP_LAND
           real, dimension(global_rt_nx, global_rt_ny):: varTmp
           call write_IO_rt_real(inVar,varTmp) 
           if(my_id .eq. IO_id) then
              iret = nf_inq_varid(ncid,varName, varid)
              iret = nf_put_vara_real(ncid, varid, (/1,1/), (/global_rt_nx,global_rt_ny/),varTmp)
           endif
#else
           iret = nf_inq_varid(ncid,varName, varid)
           if(iret .eq. 0) then
              iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ix,jx/),inVar)
           else
              write(6,*) "Error : variable not defined in rst file before write: ", varName
           endif
#endif
           
           return
        end subroutine w_rst_rt_nc2

        subroutine w_rst_rt_nc3(ncid,ix,jx,NSOIL,inVar, varName)
           implicit none
           integer:: ncid,ix,jx,varid , iret, nsoil
           character(len=*) varName
           real,dimension(ix,jx,nsoil):: inVar
           character(len=2) tmpStr
           integer k
#ifdef MPP_LAND
           real varTmp(global_rt_nx,global_rt_ny)
           do k = 1, nsoil
              call write_IO_rt_real(inVar(:,:,k),varTmp(:,:)) 
              if(my_id .eq. IO_id) then
                 if( k .lt. 10) then
                    write(tmpStr, '(i1)') k
                 else
                    write(tmpStr, '(i2)') k
                 endif
                 iret = nf_inq_varid(ncid,varName//trim(tmpStr), varid)
                 iret = nf_put_vara_real(ncid, varid, (/1,1/), (/global_rt_nx,global_rt_ny/),varTmp)
              endif
           end do
#else
           do k = 1, nsoil
                 if( k .lt. 10) then
                    write(tmpStr, '(i1)') k
                 else
                    write(tmpStr, '(i2)') k
                 endif
              iret = nf_inq_varid(ncid,varName//trim(tmpStr), varid)
              iret = nf_put_vara_real(ncid, varid, (/1,1/),(/ix,jx/),inVar(:,:,k)) 
           end do 
#endif
           return
        end subroutine w_rst_rt_nc3

        subroutine w_rst_nc2(ncid,ix,jx,inVar,varName)
           implicit none
           integer:: ncid,ix,jx,varid , iret
           character(len=*) varName
           real inVar(ix,jx)

#ifdef MPP_LAND
           real varTmp(global_nx,global_ny)
           call write_IO_real(inVar,varTmp) 
           if(my_id .eq. IO_id) then
              iret = nf_inq_varid(ncid,varName, varid)
              iret = nf_put_vara_real(ncid, varid, (/1,1/), (/global_nx,global_ny/),varTmp)
           endif
#else
           iret = nf_inq_varid(ncid,varName, varid)
           iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ix,jx/),invar)
#endif
           
           return
        end subroutine w_rst_nc2

        subroutine w_rst_nc3(ncid,ix,jx,NSOIL,inVar, varName)
           implicit none
           integer:: ncid,ix,jx,varid , iret, nsoil
           character(len=*) varName
           real inVar(ix,jx,nsoil)
           integer k
           character(len=2) tmpStr
           
#ifdef MPP_LAND
           real varTmp(global_nx,global_ny)
           do k = 1, nsoil
              call write_IO_real(inVar(:,:,k),varTmp(:,:)) 
              if(my_id .eq. IO_id) then
                 if( k .lt. 10) then
                    write(tmpStr, '(i1)') k
                 else
                    write(tmpStr, '(i2)') k
                 endif
                iret = nf_inq_varid(ncid,varName//trim(tmpStr), varid)
                iret = nf_put_vara_real(ncid, varid, (/1,1/), (/global_nx,global_ny/),varTmp)
              endif
           end do
#else
           do k = 1, nsoil
                 if( k .lt. 10) then
                    write(tmpStr, '(i1)') k
                 else
                    write(tmpStr, '(i2)') k
                 endif
             iret = nf_inq_varid(ncid,varName//trim(tmpStr), varid)
             iret = nf_put_vara_real(ncid, varid, (/1,1/), (/ix,jx/),inVar(:,:,k) )
           end do 
#endif
           return
        end subroutine w_rst_nc3

        subroutine w_rst_crt_nc1_lake(ncid,n,inVar,varName &
#ifdef MPP_LAND
                 ,nodelist     &
#endif
                  )
           implicit none
           integer:: ncid,n,varid , iret
           character(len=*) varName
           real inVar(n)
#ifdef MPP_LAND
           integer:: nodelist(n)
           if(n .eq. 0) return

           call write_lake_real(inVar,nodelist,n)          
           if(my_id .eq. IO_id) then
#endif
              iret = nf_inq_varid(ncid,varName, varid)
              iret = nf_put_vara_real(ncid, varid, (/1/), (/n/),inVar)
#ifdef MPP_LAND
           endif
#endif
           return
        end subroutine w_rst_crt_nc1_lake

        subroutine w_rst_crt_nc1(ncid,n,inVar,varName &
#ifdef MPP_LAND
                 ,map_l2g, gnlinks&
#endif
                  )
           implicit none
           integer:: ncid,n,varid , iret
           character(len=*) varName
           real inVar(n)
#ifdef MPP_LAND
           integer:: gnlinks, map_l2g(n)
           real g_var(gnlinks)
           call write_chanel_real(inVar,map_l2g,gnlinks,n,g_var)          
           if(my_id .eq. IO_id) then
              iret = nf_inq_varid(ncid,varName, varid)
              iret = nf_put_vara_real(ncid, varid, (/1/), (/gnlinks/),g_var)
#else
              iret = nf_inq_varid(ncid,varName, varid)
              iret = nf_put_vara_real(ncid, varid, (/1/), (/n/),inVar)
#endif
#ifdef MPP_LAND
           endif
#endif
           return
        end subroutine w_rst_crt_nc1

        subroutine w_rst_crt_nc1g(ncid,n,inVar,varName)
           implicit none
           integer:: ncid,n,varid , iret
           character(len=*) varName
           real inVar(n)
#ifdef MPP_LAND
           if(my_id .eq. IO_id) then
#endif
              iret = nf_inq_varid(ncid,varName, varid)
              iret = nf_put_vara_real(ncid, varid, (/1/), (/n/),inVar)
#ifdef MPP_LAND
           endif
#endif
           return
        end subroutine w_rst_crt_nc1g

   subroutine RESTART_IN_NC(inFile,did)


   implicit none
   character(len=*) inFile
   integer  :: ierr, iret,ncid, did

    integer :: i, j


#ifdef MPP_LAND
     if(IO_id .eq. my_id) then
#endif
!open a netcdf file 
    iret = nf_open(trim(inFile), NF_NOWRITE, ncid)
#ifdef MPP_LAND
    endif
    call mpp_land_bcast_int1(iret)
#endif
    if (iret /= 0) then
#ifdef HYDRO_D
       write(*,'("Problem opening file: ''", A, "''")') &
            trim(inFile)
       call hydro_stop("RESTART_IN_NC") 
#endif
    endif

#ifdef MPP_LAND
     if(IO_id .eq. my_id) then
#endif
        iret = NF_GET_ATT_INT(ncid, NF_GLOBAL, 'his_out_counts', rt_domain(did)%his_out_counts) 
        iret = NF_GET_ATT_REAL(ncid, NF_GLOBAL, 'DTCT', nlst_rt(did)%DTCT)
        iret = nf_get_att_text(ncid,NF_GLOBAL,"Since_Date",nlst_rt(did)%sincedate(1:19))
        if(iret /= 0) nlst_rt(did)%sincedate = nlst_rt(did)%startdate
        if(nlst_rt(did)%DTCT .gt. 0) then
           nlst_rt(did)%DTCT = min(nlst_rt(did)%DTCT, nlst_rt(did)%DTRT)
        else
           nlst_rt(did)%DTCT = nlst_rt(did)%DTRT
        endif
#ifdef MPP_LAND
    endif
    call mpp_land_bcast_int1(rt_domain(did)%out_counts)
    call mpp_land_bcast_real1(nlst_rt(did)%DTCT)
#endif

#ifdef HYDRO_D
     write(6,*) "nlst_rt(did)%nsoil=",nlst_rt(did)%nsoil
#endif

     if(nlst_rt(did)%rst_typ .eq. 1 ) then 
        call read_rst_nc3(ncid,rt_domain(did)%ix,rt_domain(did)%jx,nlst_rt(did)%nsoil,rt_domain(did)%stc,"stc")
        call read_rst_nc3(ncid,rt_domain(did)%ix,rt_domain(did)%jx,nlst_rt(did)%nsoil,rt_domain(did)%smc,"smc")
        call read_rst_nc3(ncid,rt_domain(did)%ix,rt_domain(did)%jx,nlst_rt(did)%nsoil,rt_domain(did)%sh2ox,"sh2ox")
        call read_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%INFXSRT,"infxsrt")
        call read_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%sfcheadrt,"sfcheadrt")
        call read_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%soldrain,"soldrain")
     endif

!yw check
 
        call read_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%SMCMAX1,"smcmax1")
        call read_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%SMCREF1,"smcref1")
        call read_rst_nc2(ncid,rt_domain(did)%ix,rt_domain(did)%jx,rt_domain(did)%SMCWLT1,"smcwlt1")


        if(nlst_rt(did)%SUBRTSWCRT.EQ.1 .OR. nlst_rt(did)%OVRTSWCRT.EQ.1 .OR. nlst_rt(did)%GWBASESWCRT .ne. 0) then
            call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%INFXSWGT,"infxswgt")
            call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%SFCHEADSUBRT,"SFCHEADSUBRT")
            call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%QBDRYRT,"QBDRYRT")
            call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%QSTRMVOLRT,"qstrmvolrt")
            call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%RETDEPRT,"RETDEPRT")
            call read_rst_rt_nc3(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,nlst_rt(did)%nsoil,rt_domain(did)%SH2OWGT,"sh2owgt")


            if(nlst_rt(did)%CHANRTSWCRT.EQ.1) then
              call read_rst_crt_stream_nc(ncid,rt_domain(did)%HLINK,rt_domain(did)%NLINKS,"hlink",rt_domain(did)%GNLINKS,rt_domain(did)%map_l2g)
              call read_rst_crt_stream_nc(ncid,rt_domain(did)%QLINK(:,1),rt_domain(did)%NLINKS,"qlink1",rt_domain(did)%GNLINKS,rt_domain(did)%map_l2g)
              call read_rst_crt_stream_nc(ncid,rt_domain(did)%QLINK(:,2),rt_domain(did)%NLINKS,"qlink2",rt_domain(did)%GNLINKS,rt_domain(did)%map_l2g)
              call read_rst_crt_stream_nc(ncid,rt_domain(did)%CVOL,rt_domain(did)%NLINKS,"cvol",rt_domain(did)%GNLINKS,rt_domain(did)%map_l2g)
              if(rt_domain(did)%NLAKES .gt. 0) then
                 call read_rst_crt_nc(ncid,rt_domain(did)%RESHT,rt_domain(did)%NLAKES,"resht")
                 call read_rst_crt_nc(ncid,rt_domain(did)%QLAKEO,rt_domain(did)%NLAKES,"qlakeo")
              endif
              call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt_domain(did)%LAKE_INFLORT,"lake_inflort")
            end if

            if(nlst_rt(did)%GWBASESWCRT.EQ.1.AND.nlst_rt(did)%GW_RESTART.NE.0 .and. rt_domain(did)%numbasns .gt. 0) then
               call read_rst_crt_nc(ncid,rt_domain(did)%z_gwsubbas,rt_domain(did)%numbasns,"z_gwsubbas")
            end if

            if(nlst_rt(did)%GWBASESWCRT.EQ.3) then
              call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,gw2d(did)%ho,"HEAD")
	    end if
        end if

        ! added by Qian, 2016.11
        if (nlst_rt(did)%OVCHRTSWCRT.EQ.1 .OR. nlst_rt(did)%INTERGWSFW .EQ.1) then
           call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt2d(did)%hlake,"HLAKE")
           call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt2d(did)%hlakeav,"HLAKEAV")
           call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt2d(did)%sinfilav,"SINFILAV")
           if (nlst_rt(did)%INTERGWSFW .EQ.1) &
           call read_rt_nc2(ncid,rt_domain(did)%ixrt,rt_domain(did)%jxrt,rt2d(did)%hbeg, "HEND")
        endif

          if(nlst_rt(did)%rstrt_swc.eq.1) then  !Switch for rest of restart accum vars...
#ifdef HYDRO_D
            print *, "1 Resetting RESTART Accumulation Variables to 0...",nlst_rt(did)%rstrt_swc
#endif
            rt_domain(did)%INFXSRT=0.
            rt_domain(did)%LAKE_INFLORT=0.
            rt_domain(did)%QSTRMVOLRT=0.
          end if

      
#ifdef MPP_LAND
        if(my_id .eq. IO_id) &
#endif
        iret =  nf_close(ncid) 
#ifdef HYDRO_D
        write(6,*) "end of RESTART_IN"
        call flush(6)
#endif

         !call check_channel(81,rt_domain(did)%QLINK(:,1),1,rt_domain(did)%NLINKS)
         !call check_channel(83,rt_domain(did)%QLINK(:,2),1,rt_domain(did)%NLINKS)
         !call check_channel(84,rt_domain(did)%HLINK,1,rt_domain(did)%NLINKS)
         !call check_channel(85,rt_domain(did)%CVOL,1,rt_domain(did)%NLINKS)
         !call hydro_stop("666666666666")

        return
        end subroutine RESTART_IN_nc

      subroutine read_rst_nc3(ncid,ix,jx,NSOIL,var,varStr)
         implicit none 
         integer ::  ix,jx,nsoil, ireg, ncid, varid, iret
         real,dimension(ix,jx,nsoil) ::  var
         character(len=*) :: varStr
         character(len=2) :: tmpStr
         integer :: n
         integer i
#ifdef MPP_LAND
         real,dimension(global_nx,global_ny) :: xtmp
#endif

         do i = 1, nsoil
#ifdef MPP_LAND
         if(my_id .eq. IO_id) then
#endif
                 if( i .lt. 10) then
                    write(tmpStr, '(i1)') i
                 else
                    write(tmpStr, '(i2)') i
                 endif
           iret = nf_inq_varid(ncid,  trim(varStr)//trim(tmpStr),  varid)
#ifdef MPP_LAND
         endif
         call mpp_land_bcast_int1(iret)
#endif

         if (iret /= 0) then
#ifdef HYDRO_D
            print*, 'variable not found: name = "', trim(varStr)//'"'
#endif
            return
         endif
#ifdef HYDRO_D
         print*, "read restart variable ", varStr//trim(tmpStr)
#endif
#ifdef MPP_LAND
         if(my_id .eq. IO_id) & 
            iret = nf_get_var_real(ncid, varid, xtmp)

            call decompose_data_real(xtmp(:,:), var(:,:,i))
#else
            iret = nf_get_var_real(ncid, varid, var(:,:,i))
#endif
         end do

         return
      end subroutine read_rst_nc3

      subroutine read_rst_nc2(ncid,ix,jx,var,varStr)
         implicit none
         integer ::  ix,jx,ireg, ncid, varid, iret
         real,dimension(ix,jx) ::  var
         character(len=*) :: varStr
#ifdef MPP_LAND
         real,dimension(global_nx,global_ny) :: xtmp 
         if(my_id .eq. IO_id) & 
#endif
           iret = nf_inq_varid(ncid,  trim(varStr),  varid)

#ifdef MPP_LAND
         call mpp_land_bcast_int1(iret)
#endif

         if (iret /= 0) then
#ifdef HYDRO_D
            print*, 'variable not found: name = "', trim(varStr)//'"'
#endif
            return
         endif
#ifdef HYDRO_D
         print*, "read restart variable ", varStr
#endif
#ifdef MPP_LAND
         if(my_id .eq. IO_id) & 
            iret = nf_get_var_real(ncid, varid, xtmp)

         call decompose_data_real(xtmp, var)
#else
            var = 0.0
            iret = nf_get_var_real(ncid, varid, var)
#endif
         return
      end subroutine read_rst_nc2

      subroutine read_rst_rt_nc3(ncid,ix,jx,NSOIL,var,varStr)
         implicit none
         integer ::  ix,jx,nsoil, ireg, ncid, varid, iret
         real,dimension(ix,jx,nsoil) ::  var
         character(len=*) :: varStr
         character(len=2) :: tmpStr
         integer i
#ifdef MPP_LAND
         real,dimension(global_rt_nx,global_rt_ny) :: xtmp
#endif
         do i = 1, nsoil
                 if( i .lt. 10) then
                    write(tmpStr, '(i1)') i
                 else
                    write(tmpStr, '(i2)') i
                 endif
#ifdef MPP_LAND
         if(my_id .eq. IO_id) & 
#endif
            iret = nf_inq_varid(ncid,  trim(varStr)//trim(tmpStr),  varid)
#ifdef MPP_LAND
         call mpp_land_bcast_int1(iret)
#endif
         if (iret /= 0) then
#ifdef HYDRO_D
            print*, 'variable not found: name = "', trim(varStr)//'"'
#endif
            return
         endif
#ifdef HYDRO_D
         print*, "read restart variable ", varStr//trim(tmpStr)
#endif
#ifdef MPP_LAND
         iret = nf_get_var_real(ncid, varid, xtmp)
            call decompose_RT_real(xtmp(:,:),var(:,:,i),global_rt_nx,global_rt_ny,ix,jx)
#else
         iret = nf_get_var_real(ncid, varid, var(:,:,i))
#endif
         end do
         return
      end subroutine read_rst_rt_nc3

      subroutine read_rst_rt_nc2(ncid,ix,jx,var,varStr)
         implicit none
         integer ::  ix,jx,ireg, ncid, varid, iret
         real,dimension(ix,jx) ::  var
         character(len=*) :: varStr
#ifdef MPP_LAND
         real,dimension(global_rt_nx,global_rt_ny) :: xtmp 
#endif
         iret = nf_inq_varid(ncid,  trim(varStr),  varid)
#ifdef MPP_LAND
         call mpp_land_bcast_int1(iret)
#endif
         if (iret /= 0) then
#ifdef HYDRO_D
            print*, 'variable not found: name = "', trim(varStr)//'"'
#endif
            return
         endif
#ifdef HYDRO_D
         print*, "read restart variable ", varStr
#endif
#ifdef MPP_LAND
         if(my_id .eq. IO_id) &   
             iret = nf_get_var_real(ncid, varid, xtmp)
         call decompose_RT_real(xtmp,var,global_rt_nx,global_rt_ny,ix,jx)
#else
            iret = nf_get_var_real(ncid, varid, var)
#endif
         return
      end subroutine read_rst_rt_nc2

      subroutine read_rt_nc2(ncid,ix,jx,var,varStr)
         implicit none
         integer ::  ix,jx, ncid, varid, iret
         real,dimension(ix,jx) ::  var
         character(len=*) :: varStr

#ifdef MPP_LAND
         real,dimension(global_rt_nx,global_rt_ny) :: xtmp
         xtmp = 0.0
#endif
            iret = nf_inq_varid(ncid,  trim(varStr),  varid)
#ifdef MPP_LAND
         call mpp_land_bcast_int1(iret)
#endif
            if (iret /= 0) then
#ifdef HYDRO_D
               print*, 'variable not found: name = "', trim(varStr)//'"'
#endif
               return
            endif
#ifdef HYDRO_D
         print*, "read restart variable ", varStr
#endif
#ifdef MPP_LAND
         if(my_id .eq. IO_id) then
            iret = nf_get_var_real(ncid, varid, xtmp)
         endif
         call decompose_RT_real(xtmp,var,global_rt_nx,global_rt_ny,ix,jx)
#else
            iret = nf_get_var_real(ncid, varid, var)
#endif
         return
      end subroutine read_rt_nc2

      subroutine read_rst_crt_nc(ncid,var,n,varStr)
         implicit none
         integer ::  ireg, ncid, varid, n, iret
         real,dimension(n) ::  var
         character(len=*) :: varStr
        
         if( n .le. 0)  return
#ifdef MPP_LAND
         if(my_id .eq. IO_id) & 
#endif
            iret = nf_inq_varid(ncid,  trim(varStr),  varid)
#ifdef MPP_LAND
         call mpp_land_bcast_int1(iret)
#endif
            if (iret /= 0) then
#ifdef HYDRO_D
               print*, 'variable not found: name = "', trim(varStr)//'"'
#endif
               return
            endif
#ifdef HYDRO_D
         print*, "read restart variable ", varStr
#endif
#ifdef MPP_LAND
         if(my_id .eq. IO_id) then
#endif
            iret = nf_get_var_real(ncid, varid, var)
#ifdef MPP_LAND
         endif
         call mpp_land_bcast_real(n,var)
#endif
         return
      end subroutine read_rst_crt_nc 

      subroutine read_rst_crt_stream_nc(ncid,var_out,n,varStr,gnlinks,map_l2g)
         implicit none
         integer ::  ncid, varid, n, iret, gnlinks
         integer, intent(in), dimension(:) :: map_l2g
         character(len=*) :: varStr
         integer :: l, g
         real,intent(out) , dimension(:) ::  var_out
#ifdef MPP_LAND
         real,dimension(gnlinks) ::  var
#else
         real,dimension(n) ::  var
#endif


#ifdef MPP_LAND
         if(my_id .eq. IO_id) & 
#endif
            iret = nf_inq_varid(ncid,  trim(varStr),  varid)
#ifdef MPP_LAND
         call mpp_land_bcast_int1(iret)
#endif
            if (iret /= 0) then
#ifdef HYDRO_D
               print*, 'variable not found: name = "', trim(varStr)//'"'
#endif
               return
            endif
#ifdef HYDRO_D
         print*, "read restart variable ", varStr
#endif
#ifdef MPP_LAND
         if(my_id .eq. IO_id) then
#endif
            var = 0.0
            iret = nf_get_var_real(ncid, varid, var)
#ifdef MPP_LAND
         endif
         call mpp_land_bcast_real(gnlinks,var)
        
         if(n .le. 0) return
         var_out = 0

         do l = 1, n
            g = map_l2g(l)
            var_out(l) = var(g)
         end do
#else
         var_out = var
#endif
         return
      end subroutine read_rst_crt_stream_nc 

      subroutine hrldas_out()
      end subroutine hrldas_out




       SUBROUTINE READ_CHROUTING1(IXRT,JXRT,fgDEM,CH_NETRT, LAKE_MSKRT, &
            FROM_NODE, TO_NODE, TYPEL, ORDER, MAXORDER, NLINKS, &
            NLAKES, MUSK, MUSX, QLINK, CHANLEN, MannN, So, ChSSlp, Bw, &
            HRZAREA, LAKEMAXH, WEIRC, WEIRL, ORIFICEC, ORIFICEA, &
            ORIFICEE, LATLAKE, LONLAKE, ELEVLAKE, &
            route_link_f, &
            route_lake_f, route_direction_f, route_order_f, &
            CHANRTSWCRT,dist, ZELEV, LAKENODE, CH_NETLNK, &
            CHANXI, CHANYJ, CHLAT, CHLON,  &
            channel_option,LATVAL,LONVAL, &
            STRMFRXSTPTS,geo_finegrid_flnm &
#ifdef MPP_LAND
               ,Link_Location &
#endif
            )
#ifdef MPP_LAND
        use module_mpp_land, only:  my_id, io_id
#endif
#include <netcdf.inc>
        INTEGER, INTENT(IN)                          :: IXRT,JXRT
        INTEGER                                      :: CHANRTSWCRT, NLINKS, NLAKES
        REAL, INTENT(IN), DIMENSION(IXRT,JXRT)       :: fgDEM
        INTEGER, DIMENSION(IXRT,JXRT)                :: DIRECTION
        INTEGER, DIMENSION(IXRT,JXRT)                :: GSTRMFRXSTPTS
        INTEGER, INTENT(IN), DIMENSION(IXRT,JXRT) :: CH_NETRT
        INTEGER, INTENT(INOUT), DIMENSION(IXRT,JXRT) :: LAKE_MSKRT
        INTEGER,                DIMENSION(IXRT,JXRT) :: GORDER  !-- gridded stream orderk
#ifdef MPP_LAND
        INTEGER,                DIMENSION(IXRT,JXRT) :: Link_Location !-- gridded stream orderk
#endif
        INTEGER                                      :: I,J,channel_option
        REAL, INTENT(OUT), DIMENSION(IXRT,JXRT)      :: LATVAL, LONVAL
        CHARACTER(len=28)                            :: dir
!Dummy inverted grids from arc


!----DJG,DNY New variables for channel and lake routing
        CHARACTER(len=155)	 :: header
        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: FROM_NODE
        REAL, INTENT(INOUT),  DIMENSION(NLINKS)      :: ZELEV
        REAL, INTENT(INOUT),  DIMENSION(NLINKS)      :: CHLAT,CHLON

        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: TYPEL
        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: TO_NODE,ORDER
        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: STRMFRXSTPTS

        INTEGER, INTENT(INOUT)                       :: MAXORDER
        REAL, INTENT(INOUT),  DIMENSION(NLINKS)      :: MUSK, MUSX !muskingum
        REAL, INTENT(INOUT),  DIMENSION(NLINKS,2)    :: QLINK  !channel flow
        REAL, INTENT(INOUT),  DIMENSION(NLINKS)      :: CHANLEN   !channel length
        REAL, INTENT(INOUT),  DIMENSION(NLINKS)      :: MannN, So !mannings N
        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: LAKENODE  ! identifies which nodes pour into which lakes
        REAL, INTENT(IN)                             :: dist(ixrt,jxrt,9)

        INTEGER, INTENT(IN), DIMENSION(IXRT,JXRT)    :: CH_NETLNK
        REAL,  DIMENSION(IXRT,JXRT)                  ::  ChSSlpG,BwG,MannNG  !channel properties on Grid
        REAL,  DIMENSION(IXRT,JXRT)                  ::  chanDepth, elrt


!-- store the location x,y location of the channel element
         INTEGER, INTENT(INOUT), DIMENSION(NLINKS)   :: CHANXI, CHANYJ

!--reservoir/lake attributes
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: HRZAREA
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: LAKEMAXH
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: WEIRC
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: WEIRL
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: ORIFICEC
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: ORIFICEA
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: ORIFICEE
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: LATLAKE,LONLAKE,ELEVLAKE
        REAL, INTENT(INOUT), DIMENSION(NLINKS)       :: ChSSlp, Bw

        CHARACTER(len=256)                           :: route_link_f
        CHARACTER(len=256)                           :: route_lake_f
        CHARACTER(len=256)                           :: route_direction_f
        CHARACTER(len=256)                           :: route_order_f
        CHARACTER(len=256)                           :: geo_finegrid_flnm
        CHARACTER(len=256)                           :: var_name

        INTEGER                                      :: tmp, cnt, ncid, iret, jj,ct
        real                                         :: gc,n
        integer :: did

        did = 1

!---------------------------------------------------------
! End Declarations
!---------------------------------------------------------
        MAXORDER = -9999
!initialize GSTRM
        GSTRMFRXSTPTS = -9999

!yw initialize the array.
        to_node =   MAXORDER
        from_node = MAXORDER
#ifdef MPP_LAND
        Link_location = MAXORDER
#endif

#ifdef HYDRO_D
        print *, "reading routing initialization files..."
        print *, "route direction", route_direction_f
        print *, "route order", route_order_f
        print *, "route linke",route_link_f
        print *, "route lake",route_lake_f
#endif

!DJG Edited code here to retrieve data from hires netcdf file....

   IF((CHANRTSWCRT.eq.1.or.CHANRTSWCRT.eq.2).AND.channel_option.eq.3) then

        var_name = "LATITUDE"
        call nreadRT2d_real  (   &
             var_name,LATVAL,ixrt,jxrt,trim(geo_finegrid_flnm))


        var_name = "LONGITUDE"
        call nreadRT2d_real(   &
             var_name,LONVAL,ixrt,jxrt,trim(geo_finegrid_flnm))

    END IF


       IF(CHANRTSWCRT.eq.1.or.CHANRTSWCRT.eq.2) then
!DJG change filename to LAKEPARM.TBL        open(unit=79,file=trim(route_link_f),          &
#ifdef MPP_LAND
         if(my_id .eq. IO_id) then 
#endif
            open(unit=79,file='LAKEPARM.TBL',          &
                  form='formatted',status='old')
#ifdef MPP_LAND
         endif
#endif
       END IF


       var_name = "LAKEGRID"
       call nreadRT2d_int(&
             var_name,LAKE_MSKRT,ixrt,jxrt,trim(geo_finegrid_flnm))

       var_name = "FLOWDIRECTION"
       call nreadRT2d_int(& 
              var_name,DIRECTION,ixrt,jxrt,trim(geo_finegrid_flnm))

       var_name = "STREAMORDER"
       call nreadRT2d_int(&
               var_name,GORDER,ixrt,jxrt,trim(geo_finegrid_flnm))


       var_name = "frxst_pts"
       call nreadRT2d_int(&
             var_name,GSTRMFRXSTPTS,ixrt,jxrt,trim(geo_finegrid_flnm))

!!!Flip y-dimension of highres grids from exported Arc files...

        var_name = "CHAN_DEPTH"
        call nreadRT2d_real(   &
             var_name,chanDepth,ixrt,jxrt,trim(geo_finegrid_flnm))

       if(nlst_rt(did)%GWBASESWCRT .eq. 3) then
             elrt = fgDEM - chanDepth
       else
             elrt = fgDEM     !ywtmp
       endif
       
       ct = 0
     
! temp fix for buggy Arc export...
        do j=1,jxrt
          do i=1,ixrt
            if(DIRECTION(i,j).eq.-128) DIRECTION(i,j)=128
          end do
        end do

      cnt =0
      if ((CHANRTSWCRT.eq.1.or.CHANRTSWCRT.eq.2).AND.channel_option .ne. 3) then ! not routing on grid, read from file
      
#ifdef MPP_LAND
       if(my_id .eq. IO_id) then
#endif
          read(79,*)  header
#ifdef MPP_LAND
       endif
#endif
       call hydro_stop("Possible Error for this code")
 
#ifdef MPP_LAND
     if(my_id .eq. IO_id) then
#endif
       do i=1,NLINKS
          read (79,*) tmp, FROM_NODE(i), TO_NODE(i), TYPEL(i),&
                   ORDER(i), QLINK(i,1), MUSK(i), MUSX(i), CHANLEN(i), &
                   MannN(i), So(i), ChSSlp(i), Bw(i), HRZAREA(i),&
                   LAKEMAXH(i), WEIRC(i), WEIRL(i), ORIFICEC(i), &
                   ORIFICEA(i),ORIFICEE(i)

           !-- hardwire QLINK
          QLINK(i,1) = 1.0
          QLINK(i,2) = QLINK(i,1)

          if (So(i).lt.0.005) So(i) = 0.005  !-- impose a minimum slope requireement
 
          if (ORDER(i) .gt. MAXORDER) then
              MAXORDER = ORDER(i)
          endif

        end do

#ifdef MPP_LAND
       endif
        call mpp_land_bcast_int(NLINKS,FROM_NODE)
        call mpp_land_bcast_int(NLINKS,TO_NODE)
        call mpp_land_bcast_int(NLINKS,TYPEL  )
        call mpp_land_bcast_int(NLINKS,ORDER  )
        call mpp_land_bcast_real(NLINKS,QLINK  )
        call mpp_land_bcast_real(NLINKS,MUSK   )
        call mpp_land_bcast_real(NLINKS,MUSX   )
        call mpp_land_bcast_real(NLINKS,CHANLEN)
        call mpp_land_bcast_real(NLINKS,MannN  )
        call mpp_land_bcast_real(NLINKS,So     )
        call mpp_land_bcast_real(NLINKS,ChSSlp )
        call mpp_land_bcast_real(NLINKS,Bw     )
        call mpp_land_bcast_real(NLINKS,HRZAREA)
        call mpp_land_bcast_real(NLINKS,LAKEMAXH)
        call mpp_land_bcast_real(NLINKS,WEIRC  )
        call mpp_land_bcast_real(NLINKS,WEIRL  )
        call mpp_land_bcast_real(NLINKS,ORIFICEC)
        call mpp_land_bcast_real(NLINKS,ORIFICEA)
        call mpp_land_bcast_real(NLINKS,ORIFICEE)
        call mpp_land_bcast_int1(MAXORDER)
        
#endif

      elseif ((CHANRTSWCRT.eq.1.or.CHANRTSWCRT.eq.2).AND.channel_option.eq.3) then  !-- handle setting up topology on the grid for diffusion scheme

#ifdef MPP_LAND
       if(my_id .eq. IO_id) then
#endif
       read(79,*)  header  !-- read the lake file
#ifdef HYDRO_D
       write(*,*) "output message: reading lake file ", header
       write(6,*) "output message: error check read file ",route_link_f
#endif
#ifdef MPP_LAND
       endif
#endif


      if (NLAKES.gt.0) then !read in only if there are lakes

#ifdef MPP_LAND
       if(my_id .eq. IO_id) then
#endif

       do i=1, NLAKES
        read (79,*) tmp, HRZAREA(i),LAKEMAXH(i), &
          WEIRC(i), WEIRL(i), ORIFICEC(i), ORIFICEA(i), ORIFICEE(i),&
           LATLAKE(i), LONLAKE(i),ELEVLAKE(i)
#ifdef HYDRO_D
        write (*,*) tmp, HRZAREA(i),LAKEMAXH(i), LATLAKE(i), LONLAKE(i),ELEVLAKE(i),NLAKES
#endif
       enddo

#ifdef MPP_LAND
       endif
          call mpp_land_bcast_real(NLAKES,HRZAREA)
          call mpp_land_bcast_real(NLAKES,LAKEMAXH)
          call mpp_land_bcast_real(NLAKES,WEIRC  )
          call mpp_land_bcast_real(NLAKES,WEIRL  )
          call mpp_land_bcast_real(NLAKES,ORIFICEC)
          call mpp_land_bcast_real(NLAKES,ORIFICEA)
          call mpp_land_bcast_real(NLAKES,ORIFICEE)
          call mpp_land_bcast_real(NLAKES,LATLAKE )
          call mpp_land_bcast_real(NLAKES,LONLAKE )
          call mpp_land_bcast_real(NLAKES,ELEVLAKE)
#endif

      end if   !end if for NLAKES >0 check

       cnt = 0 

      
     BwG = 0.0
     ChSSlpG = 0.0
     MannNG = 0.0
     TYPEL = 0
     MannN = 0.0
     Bw = 0.0
     ChSSlp = 0.0

!DJG inv       DO j = JXRT,1,-1  !rows
       DO j = 1,JXRT  !rows
        DO i = 1 ,IXRT   !colsumns
         If (CH_NETRT(i, j) .ge. 0) then !get its direction and assign its elevation and order
          If ((DIRECTION(i, j) .EQ. 64) .AND. (j + 1 .LE. JXRT) .AND. &
               (CH_NETRT(i,j+1).ge.0) ) then !North
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
             ORDER(cnt) = GORDER(i,j)
             STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
             ZELEV(cnt) = ELRT(i,j)
             MannN(cnt) = MannNG(i,j)
             ChSSlp(cnt) = ChSSlpG(i,j)
             Bw(cnt) = BwG(i,j)
             CHLAT(cnt) = LATVAL(i,j)
             CHLON(cnt) = LONVAL(i,j)
             FROM_NODE(cnt) = CH_NETLNK(i, j)
             TO_NODE(cnt) = CH_NETLNK(i, j + 1)
             CHANLEN(cnt) = dist(i,j,1)
             CHANXI(cnt) = i
             CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
          else if ((DIRECTION(i, j) .EQ. 128) .AND. (i + 1 .LE. IXRT) &
                    .AND. (j + 1 .LE. JXRT) .AND. (CH_NETRT(i+1,j+1).ge.0) ) then !North East
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
             ORDER(cnt) = GORDER(i,j)
             STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
             ZELEV(cnt) = ELRT(i,j)
             MannN(cnt) = MannNG(i,j)
             ChSSlp(cnt) = ChSSlpG(i,j)
             Bw(cnt) = BwG(i,j)
             CHLAT(cnt) = LATVAL(i,j)
             CHLON(cnt) = LONVAL(i,j)
             FROM_NODE(cnt) = CH_NETLNK(i, j)
             TO_NODE(cnt) = CH_NETLNK(i + 1, j + 1)
             CHANLEN(cnt) = dist(i,j,2)
             CHANXI(cnt) = i
             CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
          else if ((DIRECTION(i, j) .EQ. 1) .AND. (i + 1 .LE. IXRT) &
                    .AND. (CH_NETRT(i+1,j).ge.0) ) then !East
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
             ORDER(cnt) = GORDER(i,j)
             STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
             ZELEV(cnt) = ELRT(i,j)
             MannN(cnt) = MannNG(i,j)
             ChSSlp(cnt) = ChSSlpG(i,j)
             Bw(cnt) = BwG(i,j)
             CHLAT(cnt) = LATVAL(i,j)
             CHLON(cnt) = LONVAL(i,j)
             FROM_NODE(cnt) = CH_NETLNK(i, j)
             TO_NODE(cnt) = CH_NETLNK(i + 1, j)
             CHANLEN(cnt) = dist(i,j,3)
             CHANXI(cnt) = i
             CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
          else if ((DIRECTION(i, j) .EQ. 2) .AND. (i + 1 .LE. IXRT) &
                    .AND. (j - 1 .NE. 0) .AND. (CH_NETRT(i+1,j-1).ge.0) ) then !south east
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
             ORDER(cnt) = GORDER(i,j)
             STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
             ZELEV(cnt) = ELRT(i,j)
             MannN(cnt) = MannNG(i,j)
             ChSSlp(cnt) = ChSSlpG(i,j)
             Bw(cnt) = BwG(i,j)
             CHLAT(cnt) = LATVAL(i,j)
             CHLON(cnt) = LONVAL(i,j)
             FROM_NODE(cnt) = CH_NETLNK(i, j)
             TO_NODE(cnt) = CH_NETLNK(i + 1, j - 1)
             CHANLEN(cnt) = dist(i,j,4)
             CHANXI(cnt) = i
             CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
          else if ((DIRECTION(i, j) .EQ. 4) .AND. (j - 1 .NE. 0).AND.(CH_NETRT(i,j-1).ge.0) ) then !due south
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
             ORDER(cnt) = GORDER(i,j)
             STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
             ZELEV(cnt) = ELRT(i,j)
             MannN(cnt) = MannNG(i,j)
             ChSSlp(cnt) = ChSSlpG(i,j)
             Bw(cnt) = BwG(i,j)
             CHLAT(cnt) = LATVAL(i,j)
             CHLON(cnt) = LONVAL(i,j)
             FROM_NODE(cnt) = CH_NETLNK(i, j)
             TO_NODE(cnt) = CH_NETLNK(i, j - 1)
             CHANLEN(cnt) = dist(i,j,5)
             CHANXI(cnt) = i
             CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
          else if ((DIRECTION(i, j) .EQ. 8) .AND. (i - 1 .GT. 0) &
               .AND. (j - 1 .NE. 0) .AND. (CH_NETRT(i-1,j-1).ge.0)) then !south west
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
             ORDER(cnt) = GORDER(i,j)
             STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
             ZELEV(cnt) = ELRT(i,j)
             MannN(cnt) = MannNG(i,j)
             ChSSlp(cnt) = ChSSlpG(i,j)
             Bw(cnt) = BwG(i,j)
             CHLAT(cnt) = LATVAL(i,j)
             CHLON(cnt) = LONVAL(i,j)
             FROM_NODE(cnt) = CH_NETLNK(i,j)
             TO_NODE(cnt) = CH_NETLNK(i - 1, j - 1)
             CHANLEN(cnt) = dist(i,j,6)
             CHANXI(cnt) = i
             CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
          else if ((DIRECTION(i, j) .EQ. 16) .AND. (i - 1 .GT. 0).AND.(CH_NETRT(i-1,j).ge.0) ) then !West
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
             FROM_NODE(cnt) = CH_NETLNK(i, j)
             ORDER(cnt) = GORDER(i,j)
             STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
             ZELEV(cnt) = ELRT(i,j)
             MannN(cnt) = MannNG(i,j)
             ChSSlp(cnt) = ChSSlpG(i,j)
             Bw(cnt) = BwG(i,j)
             CHLAT(cnt) = LATVAL(i,j)
             CHLON(cnt) = LONVAL(i,j)
             TO_NODE(cnt) = CH_NETLNK(i - 1, j)
             CHANLEN(cnt) = dist(i,j,7)
             CHANXI(cnt) = i
             CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
          else if ((DIRECTION(i, j) .EQ. 32) .AND. (i - 1 .GT. 0) &
                    .AND. (j + 1 .LE. JXRT) .AND. (CH_NETRT(i-1,j+1).ge.0) ) then !North West
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
             ORDER(cnt) = GORDER(i,j)
             STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
             ZELEV(cnt) = ELRT(i,j)
             MannN(cnt) = MannNG(i,j)
             ChSSlp(cnt) = ChSSlpG(i,j)
             Bw(cnt) = BwG(i,j)
             CHLAT(cnt) = LATVAL(i,j)
             CHLON(cnt) = LONVAL(i,j)
             FROM_NODE(cnt) = CH_NETLNK(i, j)
             TO_NODE(cnt) = CH_NETLNK(i - 1, j + 1)
             CHANLEN(cnt) = dist(i,j,8)
             CHANXI(cnt) = i
             CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
          else 
#ifdef HYDRO_D
             print *, "NO MATCH", i,j,CH_NETLNK(i,j),DIRECTION(i,j),i + 1,j - 1 !south east
#endif
          End If

         End If !CH_NETRT check for this node

        END DO
       END DO 

#ifdef HYDRO_D
       print *, "after exiting the channel, this many nodes", cnt
       write(*,*) " " 
#endif


!Find out if the boundaries are on an edge
!DJG inv       DO j = JXRT,1,-1
       DO j = 1,JXRT
         DO i = 1 ,IXRT
          If (CH_NETRT(i, j) .ge. 0) then !get its direction

           If (((DIRECTION(i, j).EQ. 64) .AND. (j + 1 .GT. JXRT)) .OR. &        !-- 64's can only flow north
              ((DIRECTION(i, j) .EQ. 64) .and. (j < jxrt) .AND. (CH_NETRT(i,j+1) .lt. 0))) then !North
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
              ORDER(cnt) = GORDER(i,j)
              STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
              ZELEV(cnt) = ELRT(i,j)
              MannN(cnt) = MannNG(i,j)
              ChSSlp(cnt) = ChSSlpG(i,j)
              Bw(cnt) = BwG(i,j)

              CHLAT(cnt) = LATVAL(i,j)
              CHLON(cnt) = LONVAL(i,j)
              if(j+1 .GT. JXRT) then !-- an edge
               TYPEL(cnt) = 1
              elseif(LAKE_MSKRT(i,j+1).gt.0) then 
               TYPEL(cnt) = 2
               LAKENODE(cnt) = LAKE_MSKRT(i,j+1)
              else
               TYPEL(cnt) = 1 
              endif
              FROM_NODE(cnt) = CH_NETLNK(i, j)
              CHANLEN(cnt) = dist(i,j,1)
              CHANXI(cnt) = i
              CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
#ifdef HYDRO_D
              print *, "Pour Point N", TYPEL(cnt), LAKENODE(cnt), CHANLEN(cnt), cnt
#endif

           else if ( ((DIRECTION(i, j) .EQ. 128) .AND. (i + 1 .GT. IXRT))  & !-- 128's can flow out of the North or East edge
               .OR.  ((DIRECTION(i, j) .EQ. 128) .AND. (j + 1 .GT. JXRT))  & !   this is due north edge
               .OR.  ((DIRECTION(i, j) .EQ. 128) .AND. (i<ixrt .and. j<jxrt) .AND.(CH_NETRT(i + 1, j + 1).lt.0))) then !North East
!#endif
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
              ORDER(cnt) = GORDER(i,j)
              STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
              ZELEV(cnt) = ELRT(i,j)
              MannN(cnt) = MannNG(i,j)
              ChSSlp(cnt) = ChSSlpG(i,j)
              Bw(cnt) = BwG(i,j)

              CHLAT(cnt) = LATVAL(i,j)
              CHLON(cnt) = LONVAL(i,j)
              if((i+1 .GT. IXRT) .OR. (j+1 .GT. JXRT))  then ! an edge
                TYPEL(cnt) = 1
              elseif(LAKE_MSKRT(i+1,j+1).gt.0) then 
                TYPEL(cnt) = 2
                LAKENODE(cnt) = LAKE_MSKRT(i+1,j+1)
              else
                TYPEL(cnt) = 1
              endif
              FROM_NODE(cnt) = CH_NETLNK(i, j)
              CHANLEN(cnt) = dist(i,j,2)  
              CHANXI(cnt) = i
              CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
#ifdef HYDRO_D
              print *, "Pour Point NE", TYPEL(cnt), LAKENODE(cnt), CHANLEN(cnt), cnt
#endif

           else if (((DIRECTION(i, j) .EQ. 1) .AND. (i + 1 .GT. IXRT)) .OR. &    !-- 1's can only flow due east
               ((DIRECTION(i, j) .EQ. 1) .and. (i<ixrt) .AND. (CH_NETRT(i + 1, j) .lt. 0))) then !East
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
              ORDER(cnt) = GORDER(i,j)
              STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
              ZELEV(cnt) = ELRT(i,j)
              MannN(cnt) = MannNG(i,j)
              ChSSlp(cnt) = ChSSlpG(i,j)
              Bw(cnt) = BwG(i,j)
              CHLAT(cnt) = LATVAL(i,j)
              CHLON(cnt) = LONVAL(i,j)
              if(i+1 .GT. IXRT) then  !an edge
                TYPEL(cnt) = 1
              elseif(LAKE_MSKRT(i+1,j).gt.0) then 
                TYPEL(cnt) = 2
                LAKENODE(cnt) = LAKE_MSKRT(i+1,j)
              else
                TYPEL(cnt) = 1
              endif
              FROM_NODE(cnt) = CH_NETLNK(i, j)
              CHANLEN(cnt) = dist(i,j,3)
              CHANXI(cnt) = i
              CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
#ifdef HYDRO_D
              print *, "Pour Point E", TYPEL(cnt), LAKENODE(cnt), CHANLEN(cnt), cnt
#endif
           else if ( ((DIRECTION(i, j) .EQ. 2) .AND. (i + 1 .GT. IXRT))    &      !-- 2's can flow out of east or south edge
                .OR. ((DIRECTION(i, j) .EQ. 2) .AND. (j - 1 .EQ. 0))       &      !-- this is the south edge
                .OR. ((DIRECTION(i, j) .EQ. 2) .and. (i<ixrt .and. j>1) .AND.(CH_NETRT(i + 1, j - 1) .lt.0))) then !south east
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
              ORDER(cnt) = GORDER(i,j)
              STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
              ZELEV(cnt) = ELRT(i,j)
              MannN(cnt) = MannNG(i,j)
              ChSSlp(cnt) = ChSSlpG(i,j)
              Bw(cnt) = BwG(i,j)
              CHLAT(cnt) = LATVAL(i,j)
              CHLON(cnt) = LONVAL(i,j)
              if((i+1 .GT. IXRT) .OR. (j-1 .EQ. 0)) then !an edge 
                TYPEL(cnt) = 1
              elseif(LAKE_MSKRT(i+1,j-1).gt.0) then 
                TYPEL(cnt) = 2
                LAKENODE(cnt) = LAKE_MSKRT(i+1,j-1)
              else
                TYPEL(cnt) = 1
              endif
              FROM_NODE(cnt) = CH_NETLNK(i, j)
              CHANLEN(cnt) = dist(i,j,4)
              CHANXI(cnt) = i
              CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
#ifdef HYDRO_D
              print *, "Pour Point SE", TYPEL(cnt), LAKENODE(cnt), CHANLEN(cnt), cnt
#endif

           else if (((DIRECTION(i, j) .EQ. 4) .AND. (j - 1 .EQ. 0)) .OR. &       !-- 4's can only flow due south
               ((DIRECTION(i, j) .EQ. 4) .and. (j>1) .AND.(CH_NETRT(i, j - 1) .lt. 0))) then !due south
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
              ORDER(cnt) = GORDER(i,j)
              STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
              ZELEV(cnt) = ELRT(i,j)
              MannN(cnt) = MannNG(i,j)
              ChSSlp(cnt) = ChSSlpG(i,j)
              Bw(cnt) = BwG(i,j)
              CHLAT(cnt) = LATVAL(i,j)
              CHLON(cnt) = LONVAL(i,j)
              if(j-1 .EQ. 0) then !- an edge
                TYPEL(cnt) =1
              elseif(LAKE_MSKRT(i,j-1).gt.0) then 
                TYPEL(cnt) = 2
                LAKENODE(cnt) = LAKE_MSKRT(i,j-1)
              else
                TYPEL(cnt) = 1
              endif
              FROM_NODE(cnt) = CH_NETLNK(i, j)
              CHANLEN(cnt) = dist(i,j,5)
              CHANXI(cnt) = i
              CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
#ifdef HYDRO_D
              print *, "Pour Point S", TYPEL(cnt), LAKENODE(cnt), CHANLEN(cnt), cnt
#endif

          else if ( ((DIRECTION(i, j) .EQ. 8) .AND. (i - 1 .LE. 0))      &      !-- 8's can flow south or west
               .OR.  ((DIRECTION(i, j) .EQ. 8) .AND. (j - 1 .EQ. 0))      &      !-- this is the south edge
               .OR.  ((DIRECTION(i, j).EQ.8).and. (i>1 .and. j>1) .AND.(CH_NETRT(i - 1, j - 1).lt.0))) then !south west
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
              ORDER(cnt) = GORDER(i,j)
              STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
              ZELEV(cnt) = ELRT(i,j)
              MannN(cnt) = MannNG(i,j)
              ChSSlp(cnt) = ChSSlpG(i,j)
              Bw(cnt) = BwG(i,j)
              CHLAT(cnt) = LATVAL(i,j)
              CHLON(cnt) = LONVAL(i,j)
              if( (i-1 .EQ. 0) .OR. (j-1 .EQ. 0) ) then !- an edge
               TYPEL(cnt) = 1
              elseif(LAKE_MSKRT(i-1,j-1).gt.0) then 
                TYPEL(cnt) = 2
                LAKENODE(cnt) = LAKE_MSKRT(i-1,j-1)
              else
                TYPEL(cnt) = 1
              endif
              FROM_NODE(cnt) = CH_NETLNK(i, j)
              CHANLEN(cnt) = dist(i,j,6) 
              CHANXI(cnt) = i
              CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
#ifdef HYDRO_D
              print *, "Pour Point SW", TYPEL(cnt), LAKENODE(cnt), CHANLEN(cnt), cnt
#endif
           else if (((DIRECTION(i, j) .EQ. 16) .AND. (i - 1 .LE.0) ) &                  !16's can only flow due west
               .OR.((DIRECTION(i, j).EQ.16) .and. (i>1) .AND.(CH_NETRT(i - 1, j).lt.0))) then !West
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
              ORDER(cnt) = GORDER(i,j)
              STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
              ZELEV(cnt) = ELRT(i,j)
              MannN(cnt) = MannNG(i,j)
              ChSSlp(cnt) = ChSSlpG(i,j)
              Bw(cnt) = BwG(i,j)
              CHLAT(cnt) = LATVAL(i,j)
              CHLON(cnt) = LONVAL(i,j)
              if(i-1 .EQ. 0) then !-- an edge
                TYPEL(cnt) = 1
              elseif(LAKE_MSKRT(i-1,j).gt.0) then 
                TYPEL(cnt) = 2
                LAKENODE(cnt) = LAKE_MSKRT(i-1,j)
              else
                TYPEL(cnt) = 1
              endif
              FROM_NODE(cnt) = CH_NETLNK(i, j)
              CHANLEN(cnt) = dist(i,j,7)
              CHANXI(cnt) = i
              CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
#ifdef HYDRO_D
              print *, "Pour Point W", TYPEL(cnt), LAKENODE(cnt), CHANLEN(cnt), cnt
#endif

           else if ( ((DIRECTION(i, j) .EQ. 32) .AND. (i - 1 .LE. 0))      &      !-- 32's can flow either west or north
               .OR.  ((DIRECTION(i, j) .EQ. 32) .AND. (j + 1 .GT. JXRT))   &      !-- this is the north edge
!#endif
               .OR.  ((DIRECTION(i, j).EQ.32) .and. (i>1 .and. j<jxrt) .AND.(CH_NETRT(i - 1, j + 1).lt.0))) then !North West
#ifdef MPP_LAND
             cnt = CH_NETLNK(i,j)
#else
             cnt = cnt + 1
#endif
              ORDER(cnt) = GORDER(i,j)
              STRMFRXSTPTS(cnt) = GSTRMFRXSTPTS(i,j)
              ZELEV(cnt) = ELRT(i,j)
              MannN(cnt) = MannNG(i,j)
              ChSSlp(cnt) = ChSSlpG(i,j)
              Bw(cnt) = BwG(i,j)
              CHLAT(cnt) = LATVAL(i,j)
              CHLON(cnt) = LONVAL(i,j)
              if( (i-1 .EQ. 0) .OR. (j+1 .GT. JXRT)) then !-- an edge
               TYPEL(cnt) = 1
              elseif(LAKE_MSKRT(i-1,j+1).gt.0) then 
                TYPEL(cnt) = 2
                LAKENODE(cnt) = LAKE_MSKRT(i-1,j+1)
              else
                TYPEL(cnt) = 1
              endif
              FROM_NODE(cnt) = CH_NETLNK(i, j)
              CHANLEN(cnt) = dist(i,j,8)   
              CHANXI(cnt) = i
              CHANYJ(cnt) = j
#ifdef MPP_LAND
             Link_Location(i,j) = cnt
#endif
#ifdef HYDRO_D
              print *, "Pour Point NW", TYPEL(cnt), LAKENODE(cnt), CHANLEN(cnt), cnt
#endif
           endif
          endif !CH_NETRT check for this node
         END DO
       END DO 
      endif


#ifdef MPP_LAND
     if(my_id .eq. IO_id)  &
#endif
     close(79)

!       if (cnt .ne. NLINKS) then 
!#ifdef HYDRO_D
!         print *, "Apparent error in network topology", cnt, NLINKS
!         print* , "ixrt =", ixrt, "jxrt =", jxrt
!#endif
!         call hydro_stop("READ_CHROUTING")
!       endif

#ifdef MPP_LAND
#ifdef HYDRO_D
     print*, "my_id=",my_id, "cnt = ", cnt 
#endif
#endif


     do i=1,NLINKS
#ifdef HYDRO_D
      if (STRMFRXSTPTS(i).ne.-9999) print *,"Frxst_pt: ",STRMFRXSTPTS(i)
#endif
     end do
 
#ifdef MPP_LAND
     Link_location = CH_NETLNK
     call MPP_CHANNEL_COM_INT(Link_location,ixrt,jxrt,TYPEL,NLINKS,99) 
#endif
!    call check_channel(71,MannN,1,nlinks)
!    call check_channel(72,TYPEL*1.0,1,nlinks)
!    call check_channel(73,Bw*1.0,1,nlinks)
!    call check_channel(74,ChSSlp*1.0,1,nlinks)
!    call check_channel(75,CHANLEN*1.0,1,nlinks)
!    call hydro_stop("55555555")


     return

!DJG -----------------------------------------------------
   END SUBROUTINE READ_CHROUTING1




#ifdef MPP_LAND

!yw continue

       SUBROUTINE MPP_READ_CHROUTING_new(IXRT,JXRT,ELRT,CH_NETRT, LAKE_MSKRT, &
            FROM_NODE, TO_NODE, TYPEL, ORDER, MAXORDER, NLINKS, &
            NLAKES, MUSK, MUSX, QLINK, CHANLEN, MannN, So, ChSSlp, Bw, &
            HRZAREA, LAKEMAXH, WEIRC, WEIRL, ORIFICEC, ORIFICEA, &
            ORIFICEE, LATLAKE, LONLAKE, ELEVLAKE, &
            route_link_f, &
            route_lake_f, route_direction_f, route_order_f, &
            CHANRTSWCRT,dist, ZELEV, LAKENODE, CH_NETLNK, &
            CHANXI, CHANYJ, CHLAT, CHLON,  &
            channel_option,LATVAL,&
            LONVAL,STRMFRXSTPTS,geo_finegrid_flnm,g_ixrt,g_jxrt, &
            gnlinks,GCH_NETLNK, map_l2g, link_location,yw_mpp_nlinks, &
            lake_index, nlinks_index)
        implicit none
        INTEGER, INTENT(IN)                          :: IXRT,JXRT,g_IXRT,g_JXRT, GNLINKS
        INTEGER                                      :: CHANRTSWCRT, NLINKS, NLAKES
        INTEGER                                      :: I,J,channel_option
        CHARACTER(len=28)                            :: dir

!----DJG,DNY New variables for channel and lake routing
        CHARACTER(len=155)	 :: header
        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: FROM_NODE
        REAL, INTENT(INOUT),  DIMENSION(NLINKS)      :: ZELEV
        REAL, INTENT(INOUT),  DIMENSION(NLINKS)      :: CHLAT,CHLON

        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: TYPEL
        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: TO_NODE,ORDER
        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: STRMFRXSTPTS

        INTEGER, INTENT(INOUT)                       :: MAXORDER
        REAL, INTENT(INOUT),  DIMENSION(NLINKS)      :: MUSK, MUSX !muskingum
        REAL, INTENT(INOUT),  DIMENSION(NLINKS,2)    :: QLINK  !channel flow
        REAL, INTENT(INOUT),  DIMENSION(NLINKS)      :: CHANLEN   !channel length
        REAL, INTENT(INOUT),  DIMENSION(NLINKS)      :: MannN, So !mannings N
        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: LAKENODE  ! identifies which nodes pour into which lakes
        REAL, INTENT(IN)                             :: dist(ixrt,jxrt,9)
        INTEGER, INTENT(INOUT),  DIMENSION(NLINKS)   :: map_l2g


!-- store the location x,y location of the channel element
         INTEGER, INTENT(INOUT), DIMENSION(NLINKS)   :: CHANXI, CHANYJ

!--reservoir/lake attributes
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: HRZAREA
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: LAKEMAXH
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: WEIRC
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: WEIRL
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: ORIFICEC
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: ORIFICEA
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: ORIFICEE
        REAL, INTENT(INOUT),  DIMENSION(NLAKES)      :: LATLAKE,LONLAKE,ELEVLAKE
        REAL, INTENT(INOUT), DIMENSION(NLINKS)       :: ChSSlp, Bw

        CHARACTER(len=256)                           :: route_link_f
        CHARACTER(len=256)                           :: route_lake_f
        CHARACTER(len=256)                           :: route_direction_f
        CHARACTER(len=256)                           :: route_order_f
        CHARACTER(len=256)                           :: geo_finegrid_flnm
        CHARACTER(len=256)                           :: var_name

        INTEGER                                      :: tmp, cnt, ncid
        real                                         :: gc,n

        INTEGER, INTENT(IN), DIMENSION(IXRT,JXRT)    :: CH_NETLNK,GCH_NETLNK
        REAL, INTENT(IN), DIMENSION(IXRT,JXRT)       :: ELRT
        INTEGER, INTENT(IN), DIMENSION(IXRT,JXRT) :: CH_NETRT
        INTEGER, INTENT(OUT), DIMENSION(IXRT,JXRT) :: LAKE_MSKRT, link_location
        REAL, INTENT(OUT), DIMENSION(IXRT,JXRT)    :: latval,lonval
        integer :: k
        INTEGER, DIMENSION(nlinks)            :: node_table, nlinks_index
        INTEGER, DIMENSION(nlakes)            :: lake_index
        integer :: yw_mpp_nlinks , l, mpp_nlinks


!        link_location = -999


          call READ_CHROUTING1(IXRT,JXRT,ELRT,CH_NETRT, LAKE_MSKRT, &
            FROM_NODE, TO_NODE, TYPEL, ORDER, MAXORDER, NLINKS, &
            NLAKES, MUSK, MUSX, QLINK, CHANLEN, MannN, So, ChSSlp, Bw, &
            HRZAREA, LAKEMAXH, WEIRC, WEIRL, ORIFICEC, ORIFICEA, &
            ORIFICEE, LATLAKE, LONLAKE, ELEVLAKE, &
            route_link_f, &
            route_lake_f, route_direction_f, route_order_f, &
            CHANRTSWCRT,dist, ZELEV, LAKENODE, CH_NETLNK, &
            CHANXI, CHANYJ, CHLAT, CHLON,  &
            channel_option,LATVAL,LONVAL, &
            STRMFRXSTPTS,geo_finegrid_flnm &
#ifdef MPP_LAND
           ,Link_Location  &
#endif
           )


      call mpp_land_max_int1(MAXORDER)


      if(MAXORDER .eq. 0)  MAXORDER = -9999

       lake_index = -99
       do j = 1, jxrt
          do i = 1, ixrt
            if (LAKE_MSKRT(i,j) .gt. 0) then
               lake_index(LAKE_MSKRT(i,j)) = LAKE_MSKRT(i,j)
            endif
          enddo
       enddo


      CHANXI = 0
      CHANYj = 0
      do j = 1, jxrt
          do i = 1, ixrt
             if(CH_NETLNK(i,j) .gt. 0) then
               CHANXI(CH_NETLNK(i,j)) = i
               CHANYJ(CH_NETLNK(i,j)) = j
             endif
          end do
      end do

      node_table = 0
      yw_mpp_nlinks = 0
      do j = 1, jxrt
          do i = 1, ixrt
            if(CH_NETLNK(i,j) .ge. 0) then
               if( (i.eq.1) .and. (left_id .ge. 0) ) then
                    continue
               elseif ( (i.eq. ixrt) .and. (right_id .ge. 0) ) then
                    continue
               elseif ( (j.eq. 1) .and. (down_id .ge. 0) ) then
                    continue
               elseif ( (j.eq. jxrt) .and. (up_id .ge. 0) ) then
                    continue
               else
                        l = CH_NETLNK(i,j)
                        ! if(from_node(l) .gt. 0 .and. to_node(l) .gt. 0) then
                             yw_mpp_nlinks = yw_mpp_nlinks + 1
                             nlinks_index(yw_mpp_nlinks) = l
                        ! endif
               endif
            endif
          end do
      end do

       write(6,*) "nlinks = ", nlinks, " yw_mpp_nlinks = ", yw_mpp_nlinks
       call flush(6)

        call mpp_land_bcast_real(NLAKES,HRZAREA)
        call mpp_land_bcast_real(NLAKES,LAKEMAXH)
        call mpp_land_bcast_real(NLAKES,WEIRC)
        call mpp_land_bcast_real(NLAKES,WEIRL)
        call mpp_land_bcast_real(NLAKES,ORIFICEC)
        call mpp_land_bcast_real(NLAKES,ORIFICEA)
        call mpp_land_bcast_real(NLAKES,ORIFICEE)
        call mpp_land_bcast_real(NLAKES,LATLAKE)
        call mpp_land_bcast_real(NLAKES,LONLAKE)
        call mpp_land_bcast_real(NLAKES,ELEVLAKE)


        link_location = CH_NETLNK

        return 
        end SUBROUTINE MPP_READ_CHROUTING_new
#endif


#ifdef MPP_LAND
       subroutine out_day_crt(dayMean,outFile)
           implicit none
           integer :: did
           real ::  dayMean(:)
           character(len=*) :: outFile
           integer:: ywflag
           ywflag = -999
           did = 1
           if((nlst_rt(did)%olddate(12:13) .eq. "00") .and. (nlst_rt(did)%olddate(15:16) .eq. "00") ) ywflag = 99
           call mpp_land_bcast_int1(ywflag)
           if(ywflag <0) return
           ! output daily
           call out_obs_crt(did,dayMean,outFile)
       end subroutine out_day_crt

       subroutine out_obs_crt(did,dayMean,outFile)
           implicit none
           integer did, i, cnt
           real ::  dayMean(:)
           character(len=*) :: outFile
           real,dimension(rt_domain(did)%gnlinks) :: g_dayMean, chlat, chlon
           integer,dimension(rt_domain(did)%gnlinks) :: STRMFRXSTPTS
           
           g_dayMean = -999
           chlat = -999
           chlon = -999
           STRMFRXSTPTS = 0

           call write_chanel_int(RT_DOMAIN(did)%STRMFRXSTPTS,rt_domain(did)%map_l2g,rt_domain(did)%gnlinks,rt_domain(did)%nlinks,STRMFRXSTPTS)

           call write_chanel_real(dayMean,rt_domain(did)%map_l2g,rt_domain(did)%gnlinks,rt_domain(did)%nlinks,g_dayMean)

           call write_chanel_real(RT_DOMAIN(did)%CHLON,rt_domain(did)%map_l2g,rt_domain(did)%gnlinks,rt_domain(did)%nlinks,chlon)

           call write_chanel_real(RT_DOMAIN(did)%CHLAT,rt_domain(did)%map_l2g,rt_domain(did)%gnlinks,rt_domain(did)%nlinks,chlat)


           open (unit=95,file=outFile,status='unknown',position='append')
           cnt = 0
           do i = 1, rt_domain(did)%gnlinks
              if(STRMFRXSTPTS(i) .gt. 0) then
                   write(95,114) nlst_rt(did)%olddate(1:4),nlst_rt(did)%olddate(6:7),nlst_rt(did)%olddate(9:10), nlst_rt(did)%olddate(12:13), &
                         cnt,chlon(i),chlat(i),g_dayMean(i) 
                   cnt = cnt + 1
              endif
           end do
           close(95)
114 FORMAT(1x,A4,A2,A2,A2,",",I7,", ",F10.5,",",F10.5,",",F12.3)
       end subroutine out_obs_crt
#endif
       
    subroutine outPutChanInfo(fromNode,toNode,chlon,chlat)
        implicit none
        integer, dimension(:) :: fromNode,toNode
        real, dimension(:) :: chlat,chlon
        integer :: iret, nodes, i, ncid, dimid_n, varid

        nodes = size(chlon,1)         
#ifdef WRFIO_NCD_LARGE_FILE_SUPPORT
       iret = nf_create("nodeInfor.nc", IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
#else
       iret = nf_create("nodeInfor.nc", NF_CLOBBER, ncid)
#endif
       iret = nf_def_dim(ncid, "node", nodes, dimid_n)  !-- make a decimated grid
!  define the varialbes
       iret = nf_def_var(ncid,"fromNode",NF_INT,1,(/dimid_n/),varid)
       iret = nf_def_var(ncid,"toNode",NF_INT,1,(/dimid_n/),varid)
       iret = nf_def_var(ncid,"chlat",NF_FLOAT,1,(/dimid_n/),varid)
          iret = nf_put_att_text(ncid,varid,'long_name',13,'node latitude')
       iret = nf_def_var(ncid,"chlon",NF_FLOAT,1,(/dimid_n/),varid)
          iret = nf_put_att_text(ncid,varid,'long_name',14,'node longitude')
       iret = nf_enddef(ncid)
!write to the file
           iret = nf_inq_varid(ncid,"fromNode", varid)
           iret = nf_put_vara_int(ncid, varid, (/1/), (/nodes/), fromNode)
           iret = nf_inq_varid(ncid,"toNode", varid)
           iret = nf_put_vara_int(ncid, varid, (/1/), (/nodes/), toNode)
           iret = nf_inq_varid(ncid,"chlat", varid)
           iret = nf_put_vara_real(ncid, varid, (/1/), (/nodes/), chlat)
           iret = nf_inq_varid(ncid,"chlon", varid)
           iret = nf_put_vara_real(ncid, varid, (/1/), (/nodes/), chlon)
          iret = nf_close(ncid)
    end subroutine outPutChanInfo
end module module_HYDRO_io
