! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. The public interfaces listed
! must all be supported with the argument lists as indicated. Many of the interfaces
! are not required for minimal implementation (see the discussion of each
! interface and look for NULL INTERFACE). 

use types_mod,             only : r8, i8, i4, MISSING_R8, &
                                  metadatalength, obstypelength

use time_manager_mod,      only : time_type, set_time, set_calendar_type

use location_mod,          only : location_type, set_location, get_location,    &
                                  get_close_obs, get_close_state,               &
                                  convert_vertical_obs, convert_vertical_state, &
                                  LocationDims

use utilities_mod,         only : error_handler,register_module, do_nml_file,   &
                                  do_nml_term,                                  &
                                  nmlfileunit, find_namelist_in_file,           &
                                  check_namelist_read

use location_io_mod,      only :  nc_write_location_atts, nc_write_location

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file,  &
                                 nc_add_global_creation_time, nc_check,         &
                                 nc_begin_define_mode, nc_end_define_mode,      &
                                 nc_open_file_readonly, nc_close_file,          &
                                 nc_get_dimension_size, nc_get_variable_size,   &
                                 nc_get_variable, nc_put_variable,              &
                                 nc_get_variable_dimension_names,               &


use         obs_kind_mod,  only : QTY_STATE_VARIABLE,             & 
                                  QTY_SOIL_TEMPERATURE,           &
                                  QTY_SOIL_MOISTURE,              &
                                  QTY_PAR_DIRECT,                 &
                                  QTY_SOLAR_INDUCED_FLUORESCENCE, &
                                  QTY_LEAF_AREA_INDEX,            &
                                  QTY_GROSS_PRIMARY_PROD_FLUX,    &
                                  get_index_for_quantity,         &
                                  get_name_for_quantity


use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain,                         &
                                  get_index_start, get_index_end,     &
                                  get_num_domains, get_num_variables, &
                                  get_variable_name 

use default_model_mod,     only : end_model, pert_model_copies, nc_write_model_vars

use dart_time_io_mod,      only : read_model_time, write_model_time

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: get_model_size,       &
          get_state_meta_data,  &
          model_interpolate,    &
          shortest_time_between_assimilations, &
          static_init_model,    &
          init_conditions,      &
          init_time,            &
          adv_1step,            &
          nc_write_model_atts

! required by DART but passed through from another module. 
! To write model specific versions of these routines
! remove the routine from  use statement above and add your code to
! this the file.
public :: pert_model_copies,      &
          nc_write_model_vars,    &
          get_close_obs,          &
          get_close_state,        &
          end_model,              &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time, &
          write_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = "dvmdostem/model_mod.f90"

! message strings
character(len=512) :: string1, string2, string3

! static_init_model
logical, save :: module_initialized = .false.
integer     :: dom_id ! used to access the state structure
integer(i8) :: model_size
type(time_type) :: time_step

! Variable Table 
integer :: nfields
integer, parameter :: max_state_variables = 10
integer, parameter :: num_state_table_columns = 5
character(len=obstypelength) :: variable_table(max_state_variables, num_state_table_columns)

! identifiers for variable_table
integer, parameter :: VAR_NAME_INDEX = 1
integer, parameter :: VAR_QTY_INDEX = 2
integer, parameter :: VAR_UPDATE_INDEX = 3


! Things should be defined in the namelist "model_mod.nml"
character(len=32)  :: calendar = 'NOLEAP'    !!dvmdostem uses 365-year calendar (no leap year)
character(len=256) :: model_input_filename = 'model_bkgd.nc'
character(len=256) :: model_para_filename  = 'cpara_bkgd.nc'
character(len=256) :: model_grid_filename  = 'run-mask.nc'
integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
integer            :: assimilation_date = 201000101     ! DA time YYYYMMDD
!logical            :: debug = .true.        ! true = print out debugging message    
logical            :: para_1d = .true.      ! true = 1d parameter estimation
character(len=metadatalength) :: tem_variables(max_state_variables * num_state_table_columns ) = ' '

namelist /model_nml/ model_input_filename, model_para_filename, model_grid_filename, & 
                     calendar, assimilation_date, tem_variables, para_1d

!--------------------------------------------------------------------------
! Things from dvmdostem simulation file





!-------------------------------------------------------------------------
! For 2d/3d application
!! Things from dvmdostem run_mask (grid) file
!integer :: ngrid    = -1 
!integer :: nlon     = -1
!integer :: nlat     = -1
!integer :: npfts    = -1
!integer :: nlevgrnd = -1 ! Number of 'ground' levels
!integer :: nlevsoi  = -1 ! Number of 'soil' levels

!real(r8), allocatable ::  lon(:) ! used grid longitude
!real(r8), allocatable ::  lat(:) ! used grid latitude
!integer,  allocatable ::  xindx(:), yindx(:) !grid indx for used cell   

! For 1d site application
type(location_type) :: state_loc(model_size)



contains


!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.
! Can be a NULL INTERFACE for the simplest models.

subroutine static_init_model()
!-----------------------------------------------------------------
real(r8) :: x_loc
integer  :: i
integer :: iunit, io

integer :: qty_list(max_state_variables)
logical :: update_list(max_state_variables)

integer :: nvars
character(len=obstypelength) :: var_names(max_state_variables)
real(r8) :: var_ranges(max_state_variables,2)
logical  :: var_update(max_state_variables)
integer  :: var_qtys(  max_state_variables)


if ( module_initialized ) return ! only need to do this once.


module_initialized = .true.


! Read namelist information
!----------------------------------------------
!call register_module(source) !Print module information

! Read the namelist
call find_namelist_in_file("input.nml", "model_nml", iunit)
!call find_namelist_in_file("model_mod.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Output the namelist values if requested
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)


! Define time step
!-----------------------------------------------
! TODO: read time step from transistion file??

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this isn't settable at runtime
! feel free to hardcode it and not add it to a namelist.
call set_calendar_type( calendar ) 

time_step = set_time(assimilation_period_seconds, assimilation_period_days)


! Define & check the variable table
!---------------------------------------------------
 call parse_variable_table(tem_variables, nfields, variable_table, qty_list, update_list) 

! Define domain index for model variables

dom_id = add_domain(model_input_filename, nfields, &
                    var_names = variable_table(1:nfields, VAR_NAME_INDEX), &
                    kind_list = qty_list(1:nfields), &
                    update_list = update_list(1:nfields))

model_size = get_domain_size(dom_id)

! Define dimension & location
!---------------------------------------------------
! For 1d application
if (para_1d == .true.)then
  do i = 1, model_size
     x_loc = (i - 1.0_r8) / model_size
     state_loc(i) =  set_location(x_loc)
  end do
endif

! If 2d/3d application, we need to get lon,lat
!call get_grid_dims() !TODO: for future 2d/3d application, additional subroutine is 
                      !      needed for reading lon/lat/depth from run_mask.nc

! Print out all the initialization info
!--------------------------------------------------


end subroutine static_init_model

!------------------------------------------------------------------
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to 
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the 
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a 
! NULL INTERFACE.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

end subroutine adv_1step


!------------------------------------------------------------------
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

x = MISSING_R8

end subroutine init_conditions

!------------------------------------------------------------------
! Companion interface to init_conditions. Returns a time that is somehow 
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

subroutine init_time(time)

type(time_type), intent(out) :: time

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time


!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 
! This interface is required for all applications.

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size



!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.
! This interface is required for all applications.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations



!------------------------------------------------------------------
! Given a state handle, a location, and a model state variable quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is an integer that specifies the quantity of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

subroutine model_interpolate(state_handle, ens_size, location, iqty, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: iqty
real(r8),            intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,             intent(out) :: istatus(ens_size)

! model_interpolate will interpolate any variable in the DART vector to the given location.
! The first variable matching the quantity of interest will be used for the interpolation.

! The return code for successful return should be 0.
! Any positive number is an error.

! istatus =  0 ... success
! istatus =  1 ... No avaiable quantity
! istatus =  3 ... quantity not in the DART vector

! Local variables
character(len=*), parameter  :: routine = 'model_interpolate'
character(len=256) :: qty_string



if ( .not. module_initialized ) call static_init_model

! initialize status var
istatus(:) = 0
expected_obs(:) = MISSING_R8

! check quantity existence
varid = get_varid_from_kind(dom_id, iqty)

if (varid < 1) then
   istatus = 1
   return
endif


! select interpolation method for variables
select case( iqty )

    case( QTY_STATE_VARIABLE ) ! 1D testing case
       call compute_grid_value(state_handle, ens_size, location, iqty, &
                                    expected_obs, istatus)
   
  
    case( QTY_LEAF_AREA_INDEX, QTY_GROSS_PRIMARY_PROD_FLUX, &
          QTY_PAR_DIRECT) ! 2D on grid cell 
   
        call compute_gridcell_value(state_handle, ens_size, location, iqty, &
                                    expected_obs, istatus) 

!case( QTY_SOIL_MOISTURE, QTY_SOIL_TEMPERATURE) ! 3D variables
     
         
    case default

      qty_string = get_name_for_quantity(iqty)

      write(string1,*)'not written for (integer) kind ',iqty
      write(string2,*)'AKA '//trim(qty_string)
      call error_handler(E_ERR,routine,string1,source,text2=string2)
      expected_obs = MISSING_R8
      istatus = 3

end select



end subroutine model_interpolate



!------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument quantity
! (qty) can be returned if the model has more than one type of field
! (for instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

subroutine get_state_meta_data(indx, location, qty_type)

integer(i8),         intent(in)  :: indx
type(location_type), intent(out) :: location
integer,             intent(out), optional :: qty_type

!real(r8) :: lat, lon
integer :: lon_indx, lat_indx, level, local_qty


if ( .not. module_initialized ) call static_init_model


! For 2d/3d application
! from the dart index get the local variables indices
!===========================================================
!if (present(qty_type)) then
!
!   qty_type = MISSING_I
!
!   call get_model_variable_indices(indx, lon_indx, lat_indx, level, &
!            var_id=var_id, dom_id=dom_id, kind_index=qty_type)
!
!   if( qty_type == MISSING_I ) then
!      write(string1,*) '[get_state_meta_data] Cannot find DART QTY  for indx ', indx
!      write(string2,*) 'variable "'//trim(get_variable_name(dom_id, var_id))//'"'
!      call error_handler(E_ERR,'get_state_meta_data',string1, source, text2=string2)
!   endif
!endif
! set_location subroutine varies according to model dimension: "oned" "twod" "threed"
!------------------------------------------------------------------------------------------

! 3d (future use)
!location = set_location(lon(lon_indx), lat(lat_indx), real(level,r8), VERTISLEVEL)



! For 1d site application (testing only)
!==============================
location = state_loc(indx)
if (present(var_type)) var_type = QTY_STATE_VARIABLE    ! default variable quantity (quick test)

end subroutine get_state_meta_data


!------------------------------------------------------------------
! Do any initialization/setup, including reading the
! namelist values.

subroutine initialize()

integer :: iunit, io

! Print module information
call register_module(source)

! Read the namelist 
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Output the namelist values if requested
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

end subroutine initialize


!------------------------------------------------------------------
! Writes model-specific attributes to a netCDF file

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid
integer, intent(in) :: domain_id

! put file into define mode.

integer :: msize

msize = int(model_size, i4)

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )

call nc_add_global_attribute(ncid, "model", "template")

call nc_write_location_atts(ncid, msize)
call nc_end_define_mode(ncid)
call nc_write_location(ncid, state_loc, msize)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts


subroutine get_grid_dims()
!=============================================================
! read run-mask.nc and define dimension
! CCC: only used in 2d/3d application
!=============================================================
 integer :: i, j, n, nx, ny
 integer :: ncid, ndim
 integer(r8), allocatable :: run(:,:)
 real(r8), allocatable :: londata(:,:), latdata(:,:)
 character(len=*), parameter :: routine = 'read_run_mask'

 logical, parameter :: debug = .true.


! open run-mask nc file
 ncid = nc_open_file_readonly(model_grid_filename, routine)

! get dimension names & lengths
 nx = nc_get_dimension_size(ncid,'X',routine)
 ny = nc_get_dimension_size(ncid,'Y',routine)

! allocate data
 allocate(run(ny,nx))
 allocate(londata(ny,nx))
 allocate(latdata(ny,nx))

! read in data
 call nc_get_variable(ncid, 'lat', latdata, routine)
 call nc_get_variable(ncid, 'lon', londata, routine)
 call nc_get_variable(ncid, 'run', run, routine)

 ! find how many cells are used in dvmdostem (0= no use, 1 = use)
 ngrid = 0
 do i = 1,ny
    do j=1, nx
      if(run(i,j) > 0)then
        ngrid = ngrid +1
      endif
    enddo
 enddo

 if (ngrid == 0)then
    write(*,*)"========= ERROR : No cell is used in run-mask.nc ======="
    return
 else
    write(*,*)"[read_run_msk] grid cells used in run-mask : ",ngrid
 endif


! allocate output lon lat data
 allocate(lon(ngrid))
 allocate(lat(ngrid))
 allocate(xindx(ngrid))
 allocate(yindx(ngrid))


! save lon, lat, and grid index
n=0
do i = 1,ny
   do j=1, nx
      if(run(i,j) > 0)then
         n=n+1      
         xindx(n)=i
         yindx(n)=j
         lon(n)=londata(i,j)
         lat(n)=latdata(i,j)
      endif
   enddo
enddo

! DART uses longitude [0 360]
where (lon <   0.0_r8) lon = lon + 360.0_r8
where (lon > 360.0_r8) lon = lon - 360.0_r8


if(debug)then
    print*,"lon = ",lon(1:ngrid)
    print*,"lat = ",lat(1:ngrid)
endif

 ! close nc file
 call nc_close_file(ncid)

end subroutine get_grid_dims


subroutine read_model_var(ncfile,varname,epochtime,outvar)
!================================================================
! read and extract model var for the DA timestep
! return selected data and date time
!
! NOTE:
! In cmax-GPP case, the var=simulated GPP and is assumed to be Hx
! read-in transition data
!================================================================
 character(len=128), intent(in)     :: ncfile        ! model transition var file name
 character(len=128), intent(in)     :: varname       ! var name
 integer, intent(in)                :: epochtime     ! target time
 real(r8), allocatable, intent(out) :: outvar(:,:,:) ! selected data at epochtime

 integer  :: nx, ny, npft, nt
 real(r8), allocatable  :: var(:,:,:,:), time(:)

 integer :: i, tnum, indx
 integer :: ncid
 character(len=*), parameter :: routine = 'read_model_var'
 logical, parameter :: debug = .true.

! open model nc file
 ncid = nc_open_file_readonly(trim(ncfile), routine)

! get dimension names & lengths
 nx = nc_get_dimension_size(ncid,'x',routine)
 ny = nc_get_dimension_size(ncid,'y',routine)
 nt = nc_get_dimension_size(ncid,'time',routine)
 npft = nc_get_dimension_size(ncid,'pft',routine)

 ! allocate input data
 allocate(time(nt))
 allocate(var(npft,nt,ny,nx))

! allocate output data
 allocate(outvar(npft,ny,nx))


! read in data
 call nc_get_variable(ncid, 'time', time, routine)
 call nc_get_variable(ncid, trim(varname), var, routine)


! identify the target time
 do i=1, nt
   if (time(i) == epochtime) then
      outvar = var(:,i,:,:)
   endif
 enddo

! close nc file
 call nc_close_file(ncid)

end subroutine read_model_var


subroutine epochtime_to_date(epochtime, run_type, year, month, day)
!==============================================
! convert dvmdostem  epoch time to date time
!
! NOTE: Current dvm-dos-tem uses 365 days calendar for 1 year
!       NO LEAPS YEAR!! 
!----------------------------------------------
 integer, intent(in)           :: epochtime
 character(len=8),intent(in)   :: run_type   ! transient('tr') or scenario('sc')          
 integer, intent(out)          :: year, month, day
 

 character(len=20) :: date_string
 integer :: base_year, days_remaining
 integer, dimension(12) :: days_in_month

 ! Define days in months for a common year
  days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


 ! Define the based date 
 if(trim(run_type) == 'tr')then
    base_year = 1901
 elseif (trim(run_type) == 'sc') then
    base_year = 2016
 else
    write(*,*) "[epochtime_to_date] Can't regonize the run_type!"     
    return
 endif

 ! Convert epochtime to date
  ! year 
  if (epochtime > 365) then
     year = base_year + days / 365
  else
     year = base_year      
  endif

  days_remaining = MOD(epochtime, 365)
  
  month = 1
  do while (days_remaining >= days_in_month(month) .and. month <= 12)
    days_remaining = days_remaining - days_in_month(month)
    month = month + 1
  enddo
  
  day = days_remaining + 1
  
  WRITE (date_string, '(I4.4,"-",I2.2,"-",I2.2)') year, month, day

end subroutine epochtime_to_date


subroutine date_to_epochtime(year, month, day, run_type, epochtime)
!==============================================
! convert dvmdostem  epoch time to date time
!
! NOTE: Current dvm-dos-tem uses 365 days calendar for 1 year
!       NO LEAPS YEAR!!
!----------------------------------------------
 integer, intent(in)           :: year, month, day
 character(len=8),intent(in)   :: run_type   ! transient('tr') or scenario('sc')
 integer, intent(out)          :: epochtime

 character(len=20) :: date_string
 integer :: base_year, base_month, base_day, total_days
 integer, dimension(12) :: days_in_month

 ! Define days in months for a common year
  days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

 ! Define the based date
 if(trim(run_type) == 'tr')then
    base_year = 1901
 elseif (trim(run_type) == 'sc') then
    base_year = 2016
 else
    write(*,*) "[date_to_epochtime] Can't regonize the run_type!"
    return
 endif

 ! Calculate total days from based year 01/01 to the given date
  total_days = (year - base_year) * 365

  if(total_days < 0) then
     write(*,*) "[date_to_epochtime] The date need to be later than ",base_year
     return
  endif

  if (month > 1) then
     do i = 1, month - 1
         total_days = total_days + days_in_month(i)
     enddo
  endif
  total_days = total_days + day - 1  ! subtract 1 to get the complete days

  epochtime = total_days

end subroutine date_to_epochtime


subroutine parse_variable_table(state_variables, ngood, table, qty_list, update_var)
! ===================================================================================
! Checks if the namelist was filled in correctly, and the user input against the variables available in the
! input netcdf file to see if it is possible to construct the DART state vector
! 
! It will return a table with each variable:
!
! cols     1    |    2      |   3   |   4    |  5 
! -------------------------------------------------------------------- 
! var1  varname | DART KIND | range | range  | update or not      
!  :                          (min)   (max)    
!  :
! varN
!
! 6. Two options for update column:
!   'UPDATE'       => update the variable in the restart file
!   'NO_COPY_BACK' => do not copy the variable back to the restart file
!
!
!
! NOTE: Currently column 5 & 6 are not available
!
! =====================================================================================

character(len=*),  intent(inout) :: state_variables(:)
integer,           intent(out) :: ngood
character(len=*),  intent(out) :: table(:,:)
integer,           intent(out) :: qty_list(:)   ! kind number
logical,           intent(out) :: update_var(:) ! logical update

character(len=*), parameter :: routine = 'parse_variable_table'

integer :: nrows, ncols, i
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: update       ! column 5

!character(len=*) :: string1, string2, string3



if ( .not. module_initialized ) call static_init_model

nrows = size(table,1) ! how many variables
ncols = 5             ! = size(table,2)

ngood = 0

! loop over all variables and assign values for the table
varLoop : do i = 1, nrows

   varname      = trim(state_variables(ncols*i - 4))
   dartstr      = trim(state_variables(ncols*i - 3))
   minvalstring = trim(state_variables(ncols*i - 2))
   maxvalstring = trim(state_variables(ncols*i - 1))
   update       = trim(state_variables(ncols*i    ))

   call to_upper(update) ! convert a string to upper case
 
   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)
   table(i,3) = trim(minvalstring)
   table(i,4) = trim(maxvalstring)
   table(i,5) = trim(update)

   ! If the first element is empty, we have found the end of the list.
   if ( table(i,1) == ' ' ) exit varLoop
 
   ! Any other condition is an error.
   if ( any(table(i,:) == ' ') ) then
      string1 = 'input.nml &model_nml:tem_variables not fully specified'
      string2 = 'must be 5 entries per variable. Last known variable name is'
      string3 = '['//trim(table(i,1))//'] ... (without the [], naturally)'
      call error_handler(E_ERR,routine,string1,source,text2=string2,text3=string3)
   endif

  ! Make sure DART qty is valid
   qty_list(i) = get_index_for_quantity(dartstr)
   if( qty_list(i)  < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1)
   endif

   ! Make sure the update variable has a valid name
   select case (update)
      case ('UPDATE')
         update_var(i) = .true.
      case ('NO_COPY_BACK')
         update_var(i) = .false.
      case default
         write(string1,'(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
         write(string2,'(6A)') 'you provided : ', trim(varname), ', ', trim(dartstr), ', ', trim(update)
         call error_handler(E_ERR,'verify_state_variables',string1, text2=string2)
   end select
 
  ngood = ngood + 1
enddo varLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_ALLMSG,routine,string1,source,text2=string2)
endif

end subroutine parse_variable_table


subroutine compute_grid_value(state_handle, ens_size, location, qty_index,&
                              interp_val, istatus)
!----------------------------------------------------
! For single grid point to do model_interpolation
!----------------------------------------------------

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location  ! For 1d case, set in static_init_model
integer,              intent(in) :: itype
real(r8),            intent(out) :: interp_val(ens_size)
integer,             intent(out) :: istatus(ens_size)

! local variables
integer(i8) :: lower_index, upper_index
real(r8) :: lctn, lctnfrac
logical :: debug =.true.

if ( .not. module_initialized ) call static_init_model

! initialize flags
interp_val = MISSING_R8  ! the DART bad value flag
istatus    = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
! In static_init_model, we set location = (1-i)/model_size
lctn = model_size * lctn

! Set lower & upper index window to include non-integer location setting
lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > model_size) lower_index = lower_index - model_size
if(upper_index > model_size) upper_index = upper_index - model_size


lctnfrac = lctn - int(lctn) !In most of case of integer location, this term will be zero
interp_val(:) = (1.0_r8 - lctnfrac) * get_state(lower_index, state_handle) + &
                            lctnfrac  * get_state(upper_index, state_handle)

if(debug)then
    write(string1, *)'[compute_grid_value] Variable index = ',lower_index, upper_index
    write(string2, *)'at grid location (lctn, lctnfrac) = (',lctn ,',',lctnfrac,')'
    write(string3, *)'interp_val = ',interp_val   
endif


end subroutine compute_grid_value 


!---------------------------------------------------------------------------------
! Stolen from DART-CLM... modified to dvm-dos-tem version
! Jan 2024 Note: will be used in PFT patches, not for 1d single-point testing  
!
!> Each gridcell may contain values for several land units, each land unit may contain
!> several columns, each column may contain several pft's.
!> aggregates across multiple pft's. So, each gridcell value
!> is an area-weighted value of an unknown number of column-based quantities.

subroutine compute_gridcell_value(state_handle, ens_size, location, qty_index, &
                                  interp_val, istatus)

! Purpose: model_interpolation for PFTs 1D/2D variables(no vertical interpolation)

! Passed variables
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
integer,             intent(in)  :: qty_index    ! QTY in DART state needed for interpolation
real(r8),            intent(out) :: interp_val(ens_size)   ! area-weighted result
integer,             intent(out) :: istatus(ens_size)      ! error code (0 == good)

character(len=*), parameter :: routine = 'compute_gridcell_value:'

! Local storage
integer  :: varid, counter(ens_size)
integer(i8) :: index1, indexi, indexN
integer  :: gridloni,gridlatj
real(r8) :: loc_lat, loc_lon
real(r8) :: state(ens_size)
real(r8) :: total(ens_size)
real(r8) :: total_area(ens_size)
real(r8), dimension(1) :: loninds,latinds
real(r8), dimension(LocationDims) :: loc
integer :: imem
character(len=obstypelength) :: varname

if ( .not. module_initialized ) call static_init_model

interp_val = MISSING_R8  ! the DART bad value flag
istatus    = 0

loc        = get_location(location)  ! loc is in DEGREES
loc_lon    = loc(1)
loc_lat    = loc(2)
! No vertical component is considered here


!! Check if variables are in state
! initialization
varname = "missing_varname"
varid = -1

varid = get_varid_from_kind(dom_id, qty_index)
if (varid < 0) then
   istatus = 11
   return
else
   varname = get_variable_name(dom_id,varid)
endif


! TODO: below are from CLM, need to modify to TEM coordinate

! determine the grid cell for the location
latinds  = minloc(abs(LAT - loc_lat))   ! these return 'arrays' ...
loninds  = minloc(abs(LON - loc_lon))   ! these return 'arrays' ...
gridlatj = latinds(1)
gridloni = loninds(1)

if (debug > 0) then
   write(string1,*)'Working on "',trim(varname),'"'
   write(string2,*)'targetlon, lon, lon index is ',loc_lon,LON(gridloni),gridloni
   write(string3,*)'targetlat, lat, lat index is ',loc_lat,LAT(gridlatj),gridlatj
   call error_handler(E_ALLMSG,routine,string1,source,text2=string2,text3=string3)
endif


index1 = get_index_start(dom_id, varid)
indexN = get_index_end(  dom_id, varid)

counter    = 0
total      = 0.0_r8      ! temp storage for state vector
total_area = 0.0_r8      ! temp storage for area
ELEMENTS : do indexi = index1, indexN

   if (   lonixy(indexi) /=  gridloni ) cycle ELEMENTS
   if (   latjxy(indexi) /=  gridlatj ) cycle ELEMENTS
   if ( landarea(indexi) ==   0.0_r8  ) cycle ELEMENTS

   state = get_state(indexi, state_handle)

   MEMBERS: do imem = 1, ens_size

      if(state(imem) == MISSING_R8) cycle MEMBERS

      counter(imem)    = counter(imem)    + 1
      total(imem)      = total(imem)      + state(imem)*landarea(indexi)
      total_area(imem) = total_area(imem) +             landarea(indexi)


   enddo MEMBERS
enddo ELEMENTS

do imem = 1,ens_size
   if (total_area(imem) > 0.0_r8 .and. istatus(imem) == 0) then
      interp_val(imem) = total(imem) / total_area(imem)
   else
      interp_val(imem) = MISSING_R8
      istatus(imem)    = 32
   endif
enddo

end subroutine compute_gridcell_value




!===================================================================
! End of model_mod
!===================================================================
end module model_mod

