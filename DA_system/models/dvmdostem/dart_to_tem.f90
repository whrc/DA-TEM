program dart_to_tem

!-------------------------------------------------------------------------------
! purpose: Convert dvmdostem parameter file into suitable netcdf input for DART
!
! USAGE:  Edit namelist: dart.nml 
!         ./dart_to_tem
!
! author: Chu-Chun Chang May 2024
!         
!-------------------------------------------------------------------------------

! This program is built with DART module for best consistency, needed to be compiled together with DART


use        types_mod, only : r8       

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_MSG, E_ERR

use netcdf_utilities_mod
!use netcdf_utilities_mod, only : nc_check, nc_define_dimension, &
!                                 nc_define_real_scalar, & 
!                                 nc_add_attribute_to_variable, &
!                                 nc_put_variable, & 
!                                 nc_create_file, &
!                                 nc_close_file


use netcdf
 
    
implicit none
 
character(len=*), parameter :: source = 'dart_to_tem.f90'

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

real(r8):: FillValue = -9999.

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
 logical :: debug = .true.
 character(len=256) :: dart_posterior = 'dart_posterior.nc'        ! input netcdf file name
 character(len=256) :: tem_restart    = 'restart-eq.nc'       ! output netcdf file name

 namelist /dart_to_tem_nml/ dart_posterior, tem_restart, debug 


!==================================================================
! Main program
!==================================================================

 call initialize_utilities(progname='dart_to_tem')

 call find_namelist_in_file("model_to_dart.nml","dart_to_tem_nml", iunit)
 read(iunit, nml = dart_to_tem_nml, iostat = io)
 call check_namelist_read(iunit, io, "model_to_dart.nml")

 if(debug)then
   print*,'dart_posterior = ',dart_posterior
   print*,'tem_restart = ',tem_restart
 endif

!-----------------------------------------------------
! allocate arrays
!----------------------------------------------------
 
 call read_anal_vegc(trim(dart_posterior), vegc)

!---------------------------------------------------
! output as netcdf file
!---------------------------------------------------
 call add_vegc_to_nc(trim(tem_restart),'vegc', vegc)

 write(*,*) "Output netcdf file : ",trim(tem_restart)

!----------------------------------------------------
! release memory space
!---------------------------------------------------
 deallocate(vegc)
 
 call finalize_utilities('dart_to_tem')

 write(*,*) " [dart_to_tem] Complete!! You nailed it!" 


contains


subroutine read_anal_vegc(ncfile, vegc)
!================================================================
! read leaf_C, stem_C, root_C and organize to vegc
! vegc(Y, X, pftpart, pft)
! pftpart: leaf, stem, root
! leaf_C, stem_C, and root_C --> (Y, X, pft)
!================================================================
 character(len=*), intent(in)       :: ncfile        ! model transition var file name
 real(r8), allocatable, intent(out) :: vegc(:,:,:,:)

 integer :: ncid
 character(len=*), parameter :: routine = 'read_anal_vegc'
 logical, parameter :: debug = .true.


! open model nc file
 ncid = nc_open_file_readonly(trim(ncfile), routine)

! get dimension names & lengths
 nx = nc_get_dimension_size(ncid,'X',routine)
 ny = nc_get_dimension_size(ncid,'Y',routine)
 npft = nc_get_dimension_size(ncid,'pft',routine)

! read vegc
 allocate(vegc(npft, pftpart, nx, ny))
 allocate(leaf_C(npft, nx, ny))
 allocate(stem_C(npft, nx, ny))
 allocate(root_C(npft, nx, ny))

 ! read leafc stemc and rootc
 call nc_get_variable(ncid, 'leaf_C', leaf_C, routine)
 call nc_get_variable(ncid, 'stem_C', stem_C, routine)
 call nc_get_variable(ncid, 'root_C', root_C, routine)

 do i = 1,ny
   do j= 1,nx
      do k = 1, npft
        vegc(k,1,j,i)=leaf_C(k,j,i)
        vegc(k,2,j,i)=stem_C(k,j,i)
        vegc(k,3,j,i)=root_C(k,j,i)
      enddo
   enddo
 enddo

 if(debug)then
   write(*,*)"Debugging from ",routine  
   write(*,*)"vegc(1,1,1,1) = ",vegc(1,1,1,1)
   write(*,*)"leaf_C(1,1,1) = ",leaf_C(1,1,1)
   write(*,*)"vegc(1,2,1,1) = ",vegc(1,2,1,1)
   write(*,*)"stem_C(1,1,1) =",stem_C(1,1,1)    
 endif

call nc_close_file(ncid)

deallocate(leaf_C, stem_C, root_C)

end subroutine read_anal_vegc        


subroutine add_vegc_to_nc(ncfile,varname, var)
!====================================================
! Add new variable to existing netcdf file 
! design for leafc, stemc, rootc
! variable dimension = (npft, nx, ny)

character(len=*), intent(in) :: ncfile
character(len=*), intent(in) :: varname
real(r8), intent(in)         :: var(npft, pftpart, nx, ny)

! local variables
integer   :: dimid_x, dimid_y, dimid_pft, dimid_pftpart
integer   :: ncid, varid, retval
integer   :: dimids(4)



! Open the existing NetCDF file
retval = nf90_open(trim(ncfile), NF90_WRITE, ncid)
if (retval /= NF90_NOERR) then
    write(*,*) 'Error opening file : ',trim(ncfile)
    stop
endif

! Get the existing dimensions
retval = nf90_inq_dimid(ncid, 'X', dimid_x)
if (retval /= NF90_NOERR) call handle_error(retval)


retval = nf90_inq_dimid(ncid, 'Y', dimid_y)
if (retval /= NF90_NOERR) call handle_error(retval)

retval = nf90_inq_dimid(ncid, 'pft', dimid_pft)
if (retval /= NF90_NOERR) call handle_error(retval)


retval = nf90_redef(ncid) ! re-define mode

! define 
retval = nf90_def_dim(ncid, 'pftpart', pftpart, dimid_pftpart)
if (retval /= NF90_NOERR) call handle_error(retval)

!call nc_get_attribute_from_variable(ncid,'vegc','_FillValue',FillValue)

! Define the new variable
dimids = (/dimid_pft, dimid_pftpart, dimid_x, dimid_y/)

retval = nf90_def_var(ncid, trim(varname), NF90_DOUBLE, dimids, varid)
if (retval /= NF90_NOERR) call handle_error(retval)

! Add attribute to the new variable
retval = nf90_put_att(ncid, varid, '_FillValue', FillValue)
if (retval /= NF90_NOERR) call handle_error(retval)

call nc_end_define_mode(ncid) ! End define mode


! Write data to the new variable
call nc_put_variable(ncid, trim(varname), var)

! Close the NetCDF file
call nc_close_file(ncid)

end subroutine add_vegc_to_nc        



subroutine handle_error(retval)
   integer, intent(in) :: retval
   print *, nf90_strerror(retval)
   stop
end subroutine handle_error



end program dart_to_tem







