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

!--------------------------------------------------------------
! Global variables
!---------------------------------------------------------------

 integer :: cmt
 real(r8), allocatable,dimension(:)  :: cmax, nmax, kra, frg     ! dimension(PFT)
 real(r8), allocatable, dimension(:) :: cfall0, cfall1, cfall2    
 real(r8), allocatable, dimension(:) :: nfall0, nfall1, nfall2
 real(r8), allocatable, dimension(:) :: krb0, krb1, krb2 
  
 ! soil calibrated parameters
 real(r8) :: micbnup, kdcrawc, kdcsoma, kdcsompr, kdcsomcr


 character(len=6) :: str
 integer :: n, i, j, k, iunit, io
 logical :: exists
 

 !-------------------------------------------
 ! parameters that can move to namelist later
 !-------------------------------------------
 integer :: target_cmt = 13  ! Target CMT
 integer,parameter :: cmt_num = 14     ! total CMT number
 integer :: npft = 10        ! number of PFT
 integer :: hlines = 4                     ! number of headlines 
 logical :: debug = .true.
 character(len=256) :: para_infile ='cmt_calparbgc.txt'  ! parameter file for cmax
 character(len=256) :: para_outfile = 'cmax_output.nc'                 ! output netcdf file name
 character(len=256) :: run_mask_file = 'run-mask.nc'     ! run-mask.nc  

 namelist /tem_to_dart_nml/ para_infile, para_outfile, target_cmt, npft, run_mask_file,&
                            debug 


!==================================================================
! Main program
!==================================================================

 call initialize_utilities(progname='tem_to_dart')


 call find_namelist_in_file("model_to_dart.nml","tem_to_dart_nml", iunit)
 read(iunit, nml = tem_to_dart_nml, iostat = io)
 call check_namelist_read(iunit, io, "model_to_dart.nml")

 if(debug)then
   print*,'para_infile = ',para_infile
   print*,'target_cmt = ',target_cmt
   print*,'npft = ',npft
 endif

!-----------------------------------------------------
! allocate arrays
!----------------------------------------------------
 allocate(cmax(npft), nmax(npft))
 allocate(kra(npft), frg(npft))
 allocate(cfall0(npft), cfall1(npft),cfall2(npft))
 allocate(nfall0(npft), nfall1(npft), nfall2(npft))
 allocate(krb0(npft), krb1(npft), krb2(npft))

 call read_calparbgc()

!---------------------------------------------------
! output as netcdf file
!---------------------------------------------------
 call write_nc_state(trim(para_outfile))

 ! Notes for future work: pararell computing output at rank=0
 write(*,*) "Output netcdf file : ",trim(para_outfile)

!----------------------------------------------------
! release memory space
!---------------------------------------------------
deallocate(cmax, nmax)
deallocate(kra, frg, cfall0, cfall1, cfall2)
deallocate(nfall0, nfall1, nfall2, krb0, krb1, krb2)

 
 call finalize_utilities('tem_to_dart')

 write(*,*) " [tem_to_dart] Complete!! You nailed it!" 


contains

subroutine read_calparbgc()        
!==========================================================
! read parameter file 'cmt_calparbgc.txt' for target CMT
!==========================================================
 !character(len=256), intent(in) :: para_file
 

 ! check parameter file existence
   inquire (FILE = trim(para_infile), EXIST = exists)

   if (.not. exists) then
        write (*,'(2A/)') ' >> Cannot find file ', para_infile
        return
        !exit
   endif

 ! open fortran file
   open(10,file=para_infile, status='old')

 ! skip first 4 headlines   
 do i=1,hlines
   read(10,*)
 enddo

 ! start reading parameters
 ! will read until cmt_num (doesn't need to read them all for now) 
 
do n = 1,cmt_num
   read(10,*) ! divided line
   read(10,'(A6,I2)') str,cmt
   if(cmt .eq. 0)then
      read(10,*)       ! CMT00 has an extra headline     
   endif
   read(10,*)          ! PFT headline

   ! read target_cmt only
   if (cmt .eq. target_cmt)then
      read(10,*) cmax
      read(10,*) nmax
      read(10,*) cfall0
      read(10,*) cfall1
      read(10,*) cfall2
      read(10,*) nfall0
      read(10,*) nfall1
      read(10,*) nfall2
      read(10,*) kra
      read(10,*) krb0
      read(10,*) krb1
      read(10,*) krb2
      read(10,*) frg
   
      ! soil calibrated parameters
      read(10,*) str         ! soil headline
      read(10,*) micbnup
      read(10,*) kdcrawc
      read(10,*) kdcsoma
      read(10,*) kdcsompr
      read(10,*) kdcsomcr
   else
      do k =1,19
        read(10,*) str
      enddo     
   
   endif

 enddo

 if(debug)then
   print*,'Target CMT =',target_cmt   
   print*,'cmax(PFT 1:10) = ',cmax(:) 
   print*,'kdcsomcr = ',kdcsomcr 
 endif   

 close(10)

 end subroutine read_calparbgc


 subroutine write_nc_state(filename)
!===========================================================
! purpose: 
! write out traget states in netcdf
! only support for calparbgc for now
!
! NOTE: the netcdf here is using DART's netcdf interface
!       not the regular netcdf-fortran function 
!
!===========================================================
 character(len = 128), intent(in) :: filename ! output netcdf file name
 integer :: ncid, varid, iunit
 integer :: n, nsize, ns

 ! grid dimension and coordinate variables
 integer, parameter :: ndim = 1       ! dimension
 integer, parameter :: nx =1, ny = 1  ! nx*ny data 
 character(len = 128) :: varLongName

 ! self-defined array type
 type para_1D
  character(len=128) :: varname
  real(r8) :: state
 end type para_1D

 type(para_1D), allocatable :: cmax_state(:) !maximum for 10 PFTs


 ! coordinate variables 
 real :: xlon, ylat

!--------------------------------------------------------------
! Setting up
!--------------------------------------------------------------

! Check if we have avaiable cmax data
nsize = 0
do n = 1,npft
  ! loop over all PFTs but only save those with cmax > 0
  if(cmax(n) > 0)then
     nsize = nsize+1
   endif
enddo

if (nsize == 0)then
  write(*,*) "======================================================"
  write(*,*) "  WARNING : no avaiable state vector detected!!!"
  write(*,*) "  PROGRAM STOP"
  write(*,*) "======================================================"
  return
endif

!allocate data array
allocate(cmax_state(nsize))


! TODO: read geo location from run-mask.nc --> lon, lat
! for prototype, assume a number for lon and lat
 xlon = -147.85 
 ylat = 65.126
 
!------------------------------------------------------ 
! create a new nc file 
!------------------------------------------------------
 ncid = nc_create_file(filename)
 call nc_define_dimension(ncid, "location", nx)


!TODO: set time coordinate
! Setting a time frame is important for DART if you want to use 4D-EnKF
! TEM is very different from weather/ocean model ...
! currently I don't have a good idea on the time coordinate for TEM


 call nc_define_real_scalar(ncid, "lon")
 call nc_add_attribute_to_variable(ncid, 'lon', "units", 'degrees_east')
 call nc_define_real_scalar(ncid, "lat")
 call nc_add_attribute_to_variable(ncid, 'lat', "units", 'degrees_north')

 ns = 0
! Define state variables
 do n = 1, npft
 
   ! loop over all PFTs but only save those with cmax > 0
    if(cmax(n) > 0)then
       ns = ns + 1
       cmax_state(ns)%state = cmax(n)

       if(n < 10)then
           write(cmax_state(ns)%varname,'("cmax",I1)') n
           write(varLongName,'("cmax for PFT = ",I1)') n
       else
           write(cmax_state(ns)%varname,'("cmax",I2)') n
           write(varLongName,'("cmax for PFT = ",I2)') n
       endif
       
       ! define variable
       call nc_define_real_variable(ncid, trim(cmax_state(ns)%varname), (/"location"/) )
       call nc_add_attribute_to_variable(ncid, trim(cmax_state(ns)%varname), "units", 'cmax')
       call nc_add_attribute_to_variable(ncid, trim(cmax_state(ns)%varname), "long name", varLongName)
    endif

 enddo

 call nc_end_define_mode(ncid)

!--------------------------------------------------------
! put coordinate data (xlon & ylat) in
!---------------------------------------------------------
 call nc_put_variable(ncid, "lon", xlon)
 call nc_put_variable(ncid, "lat", ylat)

!-------------------------------------------------------
! putt state variables (cmax)
!--------------------------------------------------------

 do i = 1,nsize
    call nc_put_variable(ncid, trim(cmax_state(i)%varname), cmax_state(i)%state)
  
    if(debug)then
        print*,"write varname = "//trim(cmax_state(i)%varname)
        print*,"state (cmax) = ",cmax_state(i)%state    
    endif
 enddo

! close nc file
 call nc_close_file(ncid)

! release memory
 deallocate(cmax_state)


end subroutine write_nc_state

subroutine read_run_mask_nc()
!=============================================================
! TODO: read run-mask.nc
!=============================================================
  
 
end subroutine read_run_mask_nc



end program tem_to_dart







