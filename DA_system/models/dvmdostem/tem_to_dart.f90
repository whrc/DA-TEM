program tem_to_dart

!-------------------------------------------------------------------------------
! purpose: Convert dvmdostem parameter file into suitable netcdf input for DART
!
! USAGE:  Edit namelist: tem_to_dart.nml 
!         ./tem_to_dart
!
! author: Chu-Chun Chang Sep 2023
!         
!-------------------------------------------------------------------------------

! This program is built with DART module for best consistency, needed to be compiled together with DART


use        types_mod, only : r8       

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_MSG, E_ERR

use time_manager_mod, only : time_type,set_date, set_time, set_calendar_type,&
                             get_date

use netcdf_utilities_mod
!use netcdf_utilities_mod, only : nc_check, nc_define_dimension, &
!                                 nc_define_real_scalar, & 
!                                 nc_add_attribute_to_variable, &
!                                 nc_put_variable, & 
!                                 nc_create_file, &
!                                 nc_close_file


use netcdf 
    

implicit none
 
character(len=*), parameter :: source = 'tem_to_dart.f90'

!--------------------------------------------------------------------------
! Default Dimensions of dvmdostem restart file
integer :: ngrid     = -1  ! number of grid cells
integer :: nx        = -1  ! X
integer :: ny        = -1  ! Y
integer :: npft      = 10  ! number of pfts
integer :: nlevsoi   = -1  ! Number of 'soil' layers

integer :: pftpart   =  3
integer :: snowlayer =  6  ! snow layer
integer :: rootlayer = 10  ! root layer
integer :: rocklayer =  5
integer ::    fronts = 10
integer ::   prevten = 10
integer :: prevtwelve = 12


!--------------------------------------------------------------
! Global variables
!---------------------------------------------------------------
 real(r8), allocatable,dimension(:,:,:)    :: leaf_C, stem_C, root_C
 real(r8), allocatable,dimension(:,:,:,:)  :: vegc

 character(len=6) :: str
 integer :: n, i, j, k, iunit, io
 logical :: exists
 

 !-------------------------------------------
 ! parameters that can move to namelist later
 !-------------------------------------------
 logical :: debug     = .true.
 logical :: add_time  = .false.
 character(len=256) :: nc_input ='restart-eq.nc'         ! input netcdf file name
 character(len=256) :: nc_output = 'model_bkgd.nc'       ! output netcdf file name
 integer            :: input_date
 type(time_type)    :: dart_date

 namelist /tem_to_dart_nml/ nc_input, nc_output, debug, add_time, input_date 


!==================================================================
! Main program
!==================================================================

 call initialize_utilities(progname='tem_to_dart')


 call find_namelist_in_file("model_to_dart.nml","tem_to_dart_nml", iunit)
 read(iunit, nml = tem_to_dart_nml, iostat = io)
 call check_namelist_read(iunit, io, "model_to_dart.nml")

 if(debug)then
   print*,'infile = ',nc_input
   print*,'npft = ',npft
 endif

!-----------------------------------------------------
! allocate arrays
!----------------------------------------------------
 
 call read_restart_vegc(nc_input,'vegc',vegc,leaf_C, stem_C, root_C)

!---------------------------------------------------
! output as netcdf file
!---------------------------------------------------
 call add_variable_to_nc(trim(nc_output),'leaf_C', leaf_C)
 call add_variable_to_nc(trim(nc_output),'stem_C', stem_C)
 call add_variable_to_nc(trim(nc_output),'root_C', root_C)

 ! Notes for future work: pararell computing output at rank=0
 write(*,*) "Output netcdf file : ",trim(nc_output)

!-----------------------------------------------------
! add date time if neccessary
!----------------------------------------------------
 if(add_time)then
     write(*,*) "Add dvmdostem_time : ",input_date
     call YYYYMMDD_to_darttime(input_date, dart_date)
     !write(*,*) " convert to dart time :", dart_date   
     call add_time_to_nc(trim(nc_output),'dvmdostem_date', input_date)    
 endif 

!----------------------------------------------------
! release memory space
!---------------------------------------------------
deallocate(vegc)
deallocate(leaf_C, stem_C, root_C)

 
 call finalize_utilities('tem_to_dart')

 write(*,*) " [tem_to_dart] Complete!! You nailed it!" 


contains


subroutine read_restart_vegc(ncfile,varname, vegc, leaf_C, stem_C, root_C)
!================================================================
! read and extract model restart VEGC for the DA timestep
! return selected data and date time
! split VEGC into leaf_C, stem_C, and root_C 
!================================================================
 character(len=*), intent(in)     :: ncfile        ! model transition var file name
 character(len=*), intent(in)     :: varname       ! var name 
 real(r8), allocatable, intent(out) :: vegc(:,:,:,:)
 real(r8), allocatable, intent(out) :: leaf_C(:,:,:) 
 real(r8), allocatable, intent(out) :: stem_C(:,:,:)
 real(r8), allocatable, intent(out) :: root_C(:,:,:)

 integer :: ncid
 character(len=*), parameter :: routine = 'read_restart_vegc'
 logical, parameter :: debug = .true.



! open model nc file
 ncid = nc_open_file_readonly(trim(ncfile), routine)

! get dimension names & lengths
 nx = nc_get_dimension_size(ncid,'X',routine)
 ny = nc_get_dimension_size(ncid,'Y',routine)
 npft = nc_get_dimension_size(ncid,'pft',routine)
 pftpart = nc_get_dimension_size(ncid,'pftpart',routine)

! read vegc
 allocate(vegc(npft, pftpart, nx, ny))
 allocate(leaf_C(npft, nx, ny))
 allocate(stem_C(npft, nx, ny))
 allocate(root_C(npft, nx, ny))

 call nc_get_variable(ncid, trim(varname), vegc, routine)

 leaf_C = vegc(:,1,:,:)
 stem_C = vegc(:,2,:,:)
 root_C = vegc(:,3,:,:)

 if(debug)then
   write(*,*)"Debugging from ",routine  
   write(*,*)"vegc(1,1,1,1) = ",vegc(1,1,2,1)
   write(*,*)"leaf_C(1,1,1) = ",leaf_C(1,2,1)
   write(*,*)"vegc(1,2,1,1) = ",vegc(1,2,2,1)
   write(*,*)"stem_C(1,1,1) =",stem_C(1,2,1)    
 endif

call nc_close_file(ncid)


end subroutine read_restart_vegc        

subroutine add_time_to_nc(ncfile,varname, YYYYMMDD)
!====================================================
! Add new variable to existing netcdf file
! design for leafc, stemc, rootc
! variable dimension = (npft, nx, ny)

character(len=*), intent(in) :: ncfile
character(len=*), intent(in) :: varname
integer, intent(in)          :: YYYYMMDD

! local variables
integer   :: ncid, varid, retval
integer   :: FillValue

! Attribute for the new variable
!character(len=*) :: attr_name = 'units'
!character(len=*) :: attr_value = 'unitless'


! Open the existing NetCDF file
retval = nf90_open(trim(ncfile), NF90_WRITE, ncid)
if (retval /= NF90_NOERR) then
    write(*,*) 'Error opening file : ',trim(ncfile)
    stop
endif

! Define the new variable
retval = nf90_redef(ncid)
if (retval /= NF90_NOERR) call handle_error(retval)

retval = nf90_def_var(ncid, trim(varname), NF90_INT, varid)
if (retval /= NF90_NOERR) then
    print *, 'Error defining variable :',trim(varname)
    call handle_error(retval)
endif

! Add attribute to the new variable
FillValue = -9999
retval = nf90_put_att(ncid, varid, '_FillValue', FillValue)

! End define mode
call nc_end_define_mode(ncid)


! Write data to the new variable
!call nc_put_variable(ncid, trim(varname), YYYYMMDD)
retval = nf90_put_var(ncid, varid, YYYYMMDD)
if (retval /= NF90_NOERR) then
     print *, 'Error writing variable :', trim(varname)
     call handle_error(retval)
endif

! Close the NetCDF file
call nc_close_file(ncid)

end subroutine add_time_to_nc



subroutine add_variable_to_nc(ncfile,varname, var)
!====================================================
! Add new variable to existing netcdf file 
! design for leafc, stemc, rootc
! variable dimension = (npft, nx, ny)

character(len=*), intent(in) :: ncfile
character(len=*), intent(in) :: varname
real(r8), intent(in)         :: var(npft, nx, ny)

! local variables
integer   :: ncid, varid, dimid_x, dimid_y, dimid_pft, retval
integer   :: dimids(3)
real(r8)  :: FillValue

! Attribute for the new variable
!character(len=*) :: attr_name = 'units'
!character(len=*) :: attr_value = 'unitless'


! Open the existing NetCDF file
retval = nf90_open(trim(ncfile), NF90_WRITE, ncid)
if (retval /= NF90_NOERR) then
    write(*,*) 'Error opening file : ',trim(ncfile)
    stop
endif

! Get the existing dimensions
retval = nf90_inq_dimid(ncid, 'X', dimid_x)
if (retval /= NF90_NOERR) then
    print *, 'Error getting dimension x'
    call handle_error(retval)
endif

retval = nf90_inq_dimid(ncid, 'Y', dimid_y)
if (retval /= NF90_NOERR) then
    print *, 'Error getting dimension y'
    call handle_error(retval)
endif

retval = nf90_inq_dimid(ncid, 'pft', dimid_pft)
if (retval /= NF90_NOERR) then
    print *, 'Error getting dimension pft'
    call handle_error(retval)
endif

call nc_get_attribute_from_variable(ncid,'vegc','_FillValue',FillValue)


! Define the new variable
retval = nf90_redef(ncid)
if (retval /= NF90_NOERR) call handle_error(retval)
dimids = (/dimid_pft, dimid_x, dimid_y/)
retval = nf90_def_var(ncid, trim(varname), NF90_DOUBLE, dimids, varid)
if (retval /= NF90_NOERR) then
    print *, 'Error defining variable :',trim(varname)
    call handle_error(retval)
endif

! define variable
!call nc_define_real_variable(ncid, trim(cmax_state(ns)%varname), (/"location"/) )


!call nc_define_real_variable(ncid, trim(varname), dimids)

! Add attribute to the new variable
retval = nf90_put_att(ncid, varid, '_FillValue', FillValue)
if (retval /= NF90_NOERR) then
    print *, 'Error adding attribute to variable :',trim(varname)
    call handle_error(retval)
endif

! End define mode
call nc_end_define_mode(ncid)


! Write data to the new variable
call nc_put_variable(ncid, trim(varname), var)
!retval = nf90_put_var(ncid, varid, var)
!if (retval /= NF90_NOERR) then
!     print *, 'Error writing variable :', trim(varname)
!     call handle_error(retval)
!endif

! Close the NetCDF file
!retval = nf90_close(ncid)
call nc_close_file(ncid)

end subroutine add_variable_to_nc        


subroutine handle_error(retval)
   integer, intent(in) :: retval
   print *, nf90_strerror(retval)
   stop
end subroutine handle_error


subroutine YYYYMMDD_to_darttime(YYYYMMDD, dart_time)
integer, intent(in)          :: YYYYMMDD
type(time_type), intent(out) :: dart_time
integer                      :: year, month, day 
integer                      :: hour, minute, second


call set_calendar_type( 'NOLEAP' )

! Extract year, month, and day using integer arithmetic
 year = YYYYMMDD / 10000
 month = (YYYYMMDD / 100) - (year * 100)
 day = YYYYMMDD - (year * 10000) - (month * 100)
 
 minute = 0
 hour = 0
 second = 0
 write(*,*) " file year = ",year
 write(*,*) "     month = ",month
 write(*,*) " file day  = ",day
 write(*,*) " file hour = ",hour
 write(*,*) " file mins = ",minute
 write(*,*) " file secs = ",second
! convert to dart time
 dart_time = set_date(year, month, day)

end subroutine YYYYMMDD_to_darttime




subroutine read_run_mask(filename,ndim,ngrid,lon,lat)
!=============================================================
! read run-mask.nc and define dimension
! only use in 2d/3d application
!=============================================================
 character(len=128), intent(in)     :: filename  ! run-mask.nc
 integer, intent(out)               :: ndim  ! dimension (1 or 2 or 3) 
 integer, intent(out)               :: ngrid ! number of grid cell
 real(r8), allocatable, intent(out) :: lon(:), lat(:)

 integer :: i, j, k
 integer :: ncid, nx,ny
 integer(r8), allocatable :: run(:,:)
 real(r8), allocatable :: londata(:,:), latdata(:,:)
 character(len=*), parameter :: routine = 'read_run_mask'

 logical, parameter :: debug = .true.


! open run-mask nc file
 ncid = nc_open_file_readonly(trim(filename), routine)
 
! get dimension names & lengths
 nx = nc_get_dimension_size(ncid,'X',routine)
 ny = nc_get_dimension_size(ncid,'Y',routine)

! allocate data
 allocate(run(ny,nx))
 allocate(londata(ny,nx))
 allocate(latdata(ny,nx))

! The nc utility in DART is associated with the dimention during compilation
! The compiler would complain if the rank of input data is not consistent with location in quickbuild.sh  

! read in data
! call nc_get_variable(ncid, 'lat', latdata, routine)
! call nc_get_variable(ncid, 'lon', londata, routine)
! call nc_get_variable(ncid, 'run', run, routine) 

! find which cell is used in dvmdostem (0= no use, 1 = use)
! ngrid = 0
! do i = 1,ny
!    do j=1, nx
!      if(run(i,j) > 0)then
!        ngrid = ngrid +1    
!        xindx(ngrid)=i
!        yindx(ngrid)=j  
!      endif
!    enddo
! enddo

! if (ngrid == 0)then
!    write(*,*)"========= ERROR : No cell is used in run-mask.nc ======="
!    return
! else
!    write(*,*)"[read_run_msk] grid cells used in run-mask : ",ngrid     
! endif


! allocate output lon lat data
! allocate(lon(ngrid))
! allocate(lat(ngrid))

! save lon, lat 
! do k=1,ngrid
!   i=xindx(k)
!   j=yindx(k)
!   lon(k)=londata(i,j)
!   lat(k)=latdata(i,j) 
! enddo

!! define ndim (1d or 2d) 
! if(ngrid == 1) then
!     ndim = 1   ! single cell (site) 
! elseif(ngrid > 1)then
!     ndim = 2   ! multiple cells   
! else
!     write(*,*)"========= [read_run_mask] Dimension errors!  ======="    
!     return
! endif

! if(debug)then
!    print*,"lon = ",lon(1:ngrid)
!    print*,"lat = ",lat(1:ngrid)    
! endif

 ! close nc file
 call nc_close_file(ncid)

end subroutine read_run_mask



end program tem_to_dart







