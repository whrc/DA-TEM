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

use types_mod,             only : r8, i8, i4, MISSING_R8, MISSING_I,            &
                                  metadatalength, obstypelength

use time_manager_mod,      only : time_type,set_date, set_time,                 &
                                  set_calendar_type, get_date

use location_mod,          only : location_type, set_location, get_location,    &
                                  get_close_obs, get_close_state,               &
                                  convert_vertical_obs, convert_vertical_state, &
                                  LocationDims, VERTISLEVEL, VERTISHEIGHT

use utilities_mod,         only : error_handler,register_module, do_nml_file,   &
                                  do_nml_term, to_upper, file_exist,            &
                                  nmlfileunit, find_namelist_in_file,           &
                                  check_namelist_read,                          &
                                  E_ERR, E_ALLMSG

!use location_io_mod,      only :  nc_write_location_atts, nc_write_location
!use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
!                                 nc_check, nc_add_global_creation_time,        &
!                                 nc_begin_define_mode, nc_end_define_mode,     &
!                                 nc_open_file_readonly, nc_open_file_readwrite,&
!                                 nc_close_file, nc_add_attribute_to_variable,  &
!                                 nc_get_dimension_size, nc_get_variable_size,  &
!                                 nc_get_variable,                              &
!                                 nc_put_variable,                              &
!                                 nc_get_variable_dimension_names,              &
!                                 nc_define_dimension, nc_variable_exists,      &
!                                 nc_define_integer_variable, &
!                                 nc_define_real_variable, &
!                                 nc_define_double_variable, &
!                                 NF90_MAX_NAME, NF90_MAX_VAR_DIMS
use netcdf_utilities_mod


use         obs_kind_mod,  only : QTY_STATE_VARIABLE,             &
                                  QTY_SOIL_TEMPERATURE,           &
                                  QTY_SOIL_MOISTURE,              &
!                                  QTY_PAR_DIRECT,                 &
!                                  QTY_SOLAR_INDUCED_FLUORESCENCE, &
                                  QTY_LEAF_AREA_INDEX,            &
                                  QTY_GROSS_PRIMARY_PROD_FLUX,    &
                                  get_index_for_quantity,         &
                                  get_name_for_quantity


use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain,  get_domain_size,       &
                                  get_index_start, get_index_end,     &
                                  get_num_domains, get_num_variables, &
                                  get_num_dims, get_dim_name,         &
                                  get_dim_length, get_domain_size,    &     
                                  get_model_variable_indices,         &
                                  get_variable_name, get_varid_from_kind 

use default_model_mod,     only : end_model, pert_model_copies, nc_write_model_vars  

use typesizes

use netcdf


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
          read_model_time,        &
          write_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = "dvmdostem/model_mod.f90"

! message strings
character(len=512) :: string1, string2, string3

! static_init_model
logical, save :: module_initialized = .false.

!integer            :: dom_id ! used to access the state structure
integer(i8)        :: model_size
type(time_type)    :: time_step

! Missing value in dvmdostem restart file
real(r8), parameter    ::  MISSING_TEM_R8 = -9999.0_R8
integer(i8), parameter ::  MISSING_TEM_I8 = -9999_I8


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
character(len=256) :: model_restart_file = 'model_bkgd.nc'
character(len=256) :: model_para_file  = 'cpara_bkgd.nc'
character(len=256) :: model_grid_file  = 'run-mask.nc'
integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
integer            :: assimilation_date = 201000101     ! DA time YYYYMMDD
!logical            :: debug = .true.        ! true = print out debugging message    
!logical            :: para_1d = .false.      ! true = 1d parameter estimation
character(len=metadatalength) :: tem_variables(max_state_variables * num_state_table_columns ) = ' '

namelist /model_nml/ model_restart_file,          &
                     model_para_file,             &
                     model_grid_file,             & 
                     assimilation_period_days,    &
                     assimilation_period_seconds, &
                     assimilation_date,           &
                     calendar,                    &
                     tem_variables                

                     
!--------------------------------------------------------------------------
! Dimensions of dvmdostem restart file
integer :: ngrid     = -1  ! number of grid cells 
integer :: nx        = -1  ! X
integer :: ny        = -1  ! Y
integer :: npfts     = -1  ! number of pfts
integer :: nlevsoi   = -1  ! Number of 'soil' layers

! below are dimensions in restart file but currently we don't need
integer :: pftpart   =  3 
integer :: snowlayer =  6  ! snow layer
integer :: rootlayer = 10  ! root layer
integer :: rocklayer =  5 
integer ::    fronts = 10  
integer ::   prevten = 10 
integer :: prevtwelve = 12


!----------------------------------------------------------------------
! dimensions defined from run-mask.nc
integer :: nrun     = -1  ! number of activated grid cells (run = 1)
integer :: nlon     = -1  ! number of longotude  
integer :: nlat     = -1  ! number of latitude

!------------------------------------------------------------------------
! The lon and lat is at the center of each grid cells
! lon(:) and lat(:) saves the lon and lat for those "ran" gridcells

real(r8), allocatable ::  lon(:,:)             ! used grid longitude
real(r8), allocatable ::  lat(:,:)             ! used grid latitude
real(r8), allocatable ::  level(:)             ! depth distributed to 1D state vector
real(r8), allocatable ::  depth(:,:,:)         ! depth, additionally calculate from DZsoil
integer,  allocatable ::  run(:,:) 
integer,  allocatable ::  gridnlev(:,:)        ! nlev for each gridcell (X, Y) 
integer,  allocatable ::  xindx(:), yindx(:)   ! indx for lon and lat

real(r8), allocatable ::  areafrac_pft(:)      ! fraction/weights for each indx value 
real(r8), allocatable ::  varwt_pft(:,:,:)     ! fraction/weights of different pfts/areas for each grid

!----------------------------------------------------------------------
! Domain ID
integer :: dom_restart = -1


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
integer     :: i, j, iunit, io, nvars
integer     :: idom, ivar, rank
integer(i8) :: indx
integer     :: ncid, close_nc

! Variable tables
integer                      :: qty_list(max_state_variables)
logical                      :: update_list(max_state_variables)
character(len=obstypelength) :: var_names(max_state_variables)
real(r8)                     :: var_ranges(max_state_variables,2)
logical                      :: var_update(max_state_variables)
integer                      :: var_qtys(  max_state_variables)

! dimension variables
integer                      :: dimlens(NF90_MAX_VAR_DIMS)
character(len=obstypelength) :: dimnames(NF90_MAX_VAR_DIMS)

character(len=obstypelength) :: varname
character(len=*), parameter  :: routine = 'static_init_model'
logical                      :: debug = .true.


! Initialization Flag
!--------------------------------------------------------
if ( module_initialized ) return

module_initialized = .true.


! Read the namelist
!--------------------------------------------------------
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")


! Output the namelist values if requested
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)


! Set calendar and timesteps
!-----------------------------------------------------------
call set_calendar_type( calendar ) 

time_step = set_time(assimilation_period_seconds, assimilation_period_days)
! TODO: CLM uses get_time --> set_time assumes [0 86400] seconds
!       check if dvm-dos-tem  

! Define dimensions
!----------------------------------------------------------
! All soil column related dimensions(depth, gridnlev) calculated here:
 call get_restart_dims(ncid, model_restart_file, nx, ny, npfts, nlevsoi)

! Allocate matrix --> read run-mask.nc, save lon & lat
!-----------------------------------------------------------
 allocate(lon(nx,ny))
 allocate(lat(nx,ny))
 allocate(run(nx,ny))

! Define dimension & location
!---------------------------------------------------
! here we read and save lon, lat, run
 call get_grid_dims(ncid, ngrid)
 ncid = 0


! Define & check the variable table
!---------------------------------------------------
call parse_variable_table(tem_variables, nfields, variable_table, qty_list, update_list) 

! Define domain index for model variables

dom_restart = add_domain(model_restart_file, nfields, &
                    var_names = variable_table(1:nfields, VAR_NAME_INDEX), &
                    kind_list = qty_list(1:nfields), &
                    update_list = update_list(1:nfields))

model_size = get_domain_size(dom_restart)
nvars = get_num_variables(dom_restart)
idom = get_num_domains()

allocate(level(model_size))        ! level for state vector
allocate(xindx(model_size))        ! X indx
allocate(yindx(model_size))        ! Y indx 
allocate(areafrac_pft(model_size)) ! weighting for the variable within gridcell (pft only)

level(:) = 0.0_r8 ! Initialize all levels to surface

write(*,*)" model_size    = ", model_size
write(*,*)" num of var    = ", nvars
write(*,*)" num of domain = ", idom
write(*,*)" dom_restart   = ", dom_restart

! Fill the array for different dimensions 
!--------------------------------------------------
ncid = 0
! open the restart file to read-in variables that need to calculate fraction(weighting) 
! for each grid cell
ncid = nc_open_file_readonly(trim(model_restart_file), routine)

do ivar = 1, get_num_variables(dom_restart)

   varname = get_variable_name(idom,ivar)
   indx    = get_index_start(idom,ivar)
   rank    = get_num_dims(idom,ivar)

   do j = 1,rank
      dimnames(j) = get_dim_name(idom,ivar,j)
      dimlens( j) = get_dim_length(idom,ivar,j)
   enddo

   
   ! Get gridcell weighting (varwt_pft) for each variable with pft dimension 
   if ( trim(dimnames(1)) == 'pft') then   
        !! Here we assumes the variable with pft dimensions are all rank = 3
        allocate(varwt_pft(dimlens(1), dimlens(2), dimlens(3)))   
        varwt_pft(:,:,:) = 1 ! default value
        ! calculate varwt_pft
        call get_varwt_pft(ncid, varname, dimlens)
   endif


   ! Currently we only have rank=3 data in restart file 
   call fill_rank3_metadata(varname, dimnames(1:3), dimlens(1:3), indx)

   if(debug)then
         write(*,*) "====== Testing for fill_ranks_metadata ======"  
         write(*,*) "         varname = ", varname 
         write(*,*) "            rank = ", rank
         write(*,*) "            indx = ", indx
         write(*,*) "    dimnames(1:3)= ", dimnames(1:3)  
         write(*,*) "     dimlens(1:3)= ", dimlens(1:3)
         write(*,*) "varwt_pft(:,1,2) = ", varwt_pft(1:dimlens(1),2,1)
   endif

enddo 

call nc_close_file(ncid)

end subroutine static_init_model



subroutine fill_rank3_metadata(varname, dim_name, dim_length, indx)
!---------------------------------------------------------------
! Fill the rank 3 metadata array
!
! Rank 3 VARIABLES in restart file:
! lai   (Y, X, pft)
! TSsoil(Y, X, soillayer)
!
!! lai(Y, X, pft) --> in ncdump
!! lai(pft, X, Y) --> in fortran/DART


character(len=*), intent(in)    :: varname
character(len=*), intent(in)    :: dim_name(3)
integer,          intent(in)    :: dim_length(3)
integer(i8),      intent(inout) :: indx

integer                         :: i, j, k
logical, parameter              :: debug = .true.


if(debug)then
   write(*,*)'===== debugging from [fill_rank3_metadata] ====='     
   write(*,*)'variable name: ',trim(varname)
   write(*,*)'dimension 1 (i) ',dim_name(1),dim_length(1)
   write(*,*)'dimension 2 (j) ',dim_name(2),dim_length(2)
   write(*,*)'dimension 3 (k) ',dim_name(3),dim_length(3)   
   write(*,*)'          indx =',indx
endif


! here the input indx is the starting point of the variable in the state vector
! so we don't need to initialize it from zero
SELECT CASE ( trim(dim_name(1)) )
   CASE ("pft")
      write(*,*) " Dimension with PFTs"
      ! For dimension of PFTs, calculate "areafrac_pft" 
      do i = 1, dim_length(3) ! loop Y
        do j = 1, dim_length(2) ! loop X
           do k = 1, dim_length(1) !loop pft
              level(indx) = 0.0  ! at surface
              xindx(indx) = j
              yindx(indx) = i
              areafrac_pft(indx) = varwt_pft(k, j, i)
              indx = indx + 1
           enddo
        enddo
      enddo
      
   CASE ("soillayer")
      write(*,*) " Dimension with vertical soil layer"
      do i = 1, dim_length(3) ! loop Y
         do j = 1, dim_length(2) ! loop X
            do k = 1, dim_length(1) ! loop soillayer
               level(indx) = depth(k, j, i) 
               xindx(indx) = j
               yindx(indx) = i
               areafrac_pft(indx) = 1
               indx = indx + 1
           enddo
        enddo
      enddo


   CASE DEFAULT
      write(*,*)'unsupported vertical dimension name "'//trim(dim_name(1))//'"'
           
END SELECT


end subroutine fill_rank3_metadata        




!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 
! This interface is required for all applications.

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.
! This interface is required for all applications.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model 

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

! Local variables
real(r8)                     :: loc_array(LocationDims)
real(r8)                     :: llon, llat, llev ! local lon, lat, lev(can be pft or depth)
character(len=*), parameter  :: routine = 'model_interpolate'
character(len=256)           :: qty_string
integer                      :: varid
logical, parameter           :: debug = .true. 


if ( .not. module_initialized ) call static_init_model

! initialize status var
istatus(:) = 0
expected_obs(:) = MISSING_R8

! check location
loc_array = get_location(location)
llon      = loc_array(1)
llat      = loc_array(2)
llev      = loc_array(3)

if(debug)then
    write(*,*) '=============== Debugging from [model_interpolate] ==============='    
    write(*,*) 'llon, llat, llev = ',llon, llat, llev    
endif

! check quantity existence
varid = get_varid_from_kind(dom_restart, iqty)

if (varid < 1) then
   istatus = 1
   return
endif


! select interpolation method for variables
select case( iqty )
     
    case( QTY_LEAF_AREA_INDEX, QTY_GROSS_PRIMARY_PROD_FLUX) ! 2D on grid cell 
        write(*,*) '2D Grid cell interpolation'
        call compute_gridcell_value(state_handle, ens_size, location, iqty, &
                                    expected_obs, istatus) 

    case( QTY_SOIL_MOISTURE, QTY_SOIL_TEMPERATURE) ! 3D variables            
        write(*,*) '3D interpolation'


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

! local variables
real(r8) :: glon, glat 
integer  :: var_id, dom_id
integer  :: ip, jp, kp, local_qty
logical, parameter :: debug = .true.


if ( .not. module_initialized ) call static_init_model


! from the dart index get the local variables indices
!===========================================================
if (present(qty_type)) then

   qty_type = MISSING_I

   call get_model_variable_indices(indx, ip, jp, kp, &
            var_id=var_id, dom_id=dom_id, kind_index=qty_type)

   if( qty_type == MISSING_I ) then
      write(string1,*) '[get_state_meta_data] Cannot find DART QTY  for indx ', indx
      write(string2,*) 'variable "'//trim(get_variable_name(dom_id, var_id))//'"'
      call error_handler(E_ERR,'get_state_meta_data',string1, source, text2=string2)
   endif
endif
! set_location subroutine varies according to model dimension: "oned" "twod" "threed"
!------------------------------------------------------------------------------------------
glon = lon(ip, jp)
glat = lat(ip, jp)
location = set_location(glon, glat, level(indx), VERTISLEVEL)

if(debug)then
        write(*,*) '[Debugging from get_state_meta_data]'
        write(*,*) 'ip, jp, indx = ',ip, jp, indx  
        write(*,*) 'lon, lat, level = ',glon, glat, level(indx)
endif


! For 1d site application (testing only)
!==============================
!location = state_loc(indx)
!if (present(qty_type)) qty_type = QTY_STATE_VARIABLE    ! default variable quantity (quick test)

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

!call nc_write_location_atts(ncid, msize)
call nc_end_define_mode(ncid)
!call nc_write_location(ncid, state_loc, msize)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts


subroutine get_restart_dims(ncid, filename, nx, ny, npfts, nlevsoi)
!==================================================================
! Read restart file dimensions and calculate depth for each level
!-=================================================================
 integer,intent(inout)          :: ncid
 character(len=*),intent(in)    :: filename ! restart file
 integer,intent(out)            :: nx, ny
 integer,intent(out)            :: npfts
 integer,intent(out)            :: nlevsoi
 
 integer                        :: i,j,k
 real(r8), allocatable          :: soil_thickness(:,:,:)
 real(r8)                       :: nsc



 character(len=*), parameter :: routine = 'get_restart_dims'
 logical, parameter :: debug = .true.

 ! read in dimension
 ncid = nc_open_file_readonly(filename, routine)

 ! get dimension names & lengths
 nx = nc_get_dimension_size(ncid,'X',routine)
 ny = nc_get_dimension_size(ncid,'Y',routine)
 npfts = nc_get_dimension_size(ncid,'pft',routine)
 nlevsoi = nc_get_dimension_size(ncid,'soillayer',routine)


 allocate(soil_thickness(nlevsoi, nx, ny))
 allocate(depth(nlevsoi, nx,ny))
 allocate(gridnlev(nx, ny))

! read DZsoil (thickness of soil layer)
 call nc_get_variable(ncid, 'DZsoil', soil_thickness, routine)

! read TYPEsoil (the type of soil)
! call nc_get_variable(ncid, 'TYPEsoil', var, routine)

gridnlev=0

! calculate the depth of soil
 do i = 1, nx
    do j =1, ny
       if(soil_thickness(1,i,j) == MISSING_R8)then
         ! unused column 
         depth(:, i, j) = MISSING_R8
         gridnlev(i, j) = 0
       else        
          nsc = 0.0       
          do k = 1, nlevsoi             
             if(k == 1 )then
                 depth(k, i, j)= 0.0
             else
                 nsc = nsc + soil_thickness(k-1, i, j)
                 depth(k, i, j) = nsc
              endif
              ! record number of used levels for each grid cell
              if(depth (k, i, j) >= 0)then
                 gridnlev(i, j) = gridnlev(i, j) + 1      
              endif

             if(debug)then              
                write(*,*) " depth (k,i,j) = ", depth(k, i, j)             
             endif
       enddo
    enddo
 enddo

 deallocate(soil_thickness)

 if(debug)then
     write(*,*)'[get_restart_dims] Dimension Outputs:'
     write(*,*)'nx       = ',nx
     write(*,*)'ny       = ',ny
     write(*,*)'nlevsoi  = ',nlevsoi
     write(*,*)'gridnlev = ',gridnlev(2,:)
     write(*,*)'npfts    = ',npfts
 endif

 call nc_close_file(ncid, routine)
 ncid = 0

end subroutine get_restart_dims


subroutine get_grid_dims(ncid, ngrid)
!=============================================================
! read run-mask.nc and define dimension
! CCC: only used in 2d/3d application
!=============================================================
 integer,intent(inout)       :: ncid
 integer, intent(out)        :: ngrid
 integer                     :: i, j, nxi, nyi
! integer(i8), allocatable    :: run_int64(:,:)
 integer, allocatable        :: ix(:), iy(:)

 character(len=*), parameter :: routine = 'get_grid_dims'
 logical, parameter :: debug = .true.


! open run-mask nc file
 ncid = nc_open_file_readonly(model_grid_file, routine)

! get dimension names & lengths
 nxi = nc_get_dimension_size(ncid,'X',routine)
 nyi = nc_get_dimension_size(ncid,'Y',routine)

 !allocate(run(nxi,nyi))
 allocate(ix(nxi*nyi), iy(nxi*nyi))

! read in data
 call nc_get_variable(ncid, 'lat', lat, routine)
 call nc_get_variable(ncid, 'lon', lon, routine)
 call nc_get_variable(ncid, 'run', run, routine)
! call nc_get_int_2d(ncid, 'run', run_int64)

! DART uses longitude [0 360]
 where (lon <   0.0_r8) lon = lon + 360.0_r8
 where (lon > 360.0_r8) lon = lon - 360.0_r8

 ngrid = 0
 do i = 1,nxi
    do j=1, nyi
      if(run(i,j) > 0)then
        ngrid = ngrid +1
        ix(ngrid) = i
        iy(ngrid) = j
      endif
    enddo
 enddo

 if (ngrid == 0)then
    write(*,*)"========= ERROR : No cell is used in run-mask.nc ======="
    return
 else
    write(*,*)"[get_grid_dims] grid cells used in run-mask : ",ngrid
    if(debug)then
       write(*,*)"lon = ",lon(ix(ngrid),iy(ngrid))
       write(*,*)"lat = ",lat(ix(ngrid),iy(ngrid))
    endif
 endif

 ! close nc file
 call nc_close_file(ncid)
 ncid = 0
! deallocate(run_int64)
 deallocate(ix)
 deallocate(iy)

end subroutine get_grid_dims


subroutine get_varwt_pft(ncid, varname, dimlens)
! Calculate the variable weight of different PFTs within gridcell
integer, intent(inout)          :: ncid
integer, intent(in)             :: dimlens(3)
character(len=*), intent(in)    :: varname

! local variables
real(r8), allocatable    :: var(:,:,:)
real(r8), allocatable    :: total_var(:,:)
integer                  :: i, j, k

character(len=*), parameter :: routine = 'get_varwt_pft'

! lai(Y, X, pft) ncdump = lai(pft, X, Y) in fortran


 allocate(var(dimlens(1),dimlens(2),dimlens(3)))
 allocate(total_var(dimlens(2),dimlens(3)))

 call nc_get_variable(ncid, trim(varname), var, routine)
 
 total_var(:,:) = 0

 ! calculate variable weighting
 do i = 1, dimlens(2)
    do j = 1, dimlens(3)
       do k =1, dimlens(1) !PFT loop
            if (var(k,i,j) > -999)then ! skip unused pfts
               total_var(i,j) = total_var(i,j) + var(k,i,j)
            endif
       enddo
    enddo
 enddo

! write(*,*) 'total_var = ', total_var

 ! default weighting
 varwt_pft = 1.0

 ! calculate variable weighting
 do i = 1, dimlens(2)
    do j = 1, dimlens(3)
       do k =1, dimlens(1) !PFT loop
            if (var(k,i,j) > -999)then ! skip unused pfts
               varwt_pft(k,i,j)= var(k,i,j)/total_var(i,j)
            endif  
       enddo
    enddo
 enddo

 !write(*,*) 'varwt_pft(k,i,j) = ', varwt_pft

 deallocate(var)
 deallocate(total_var)
  

end subroutine get_varwt_pft


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



function read_model_time(filename)
! Function to read dvm-dos-tem file time
! CCC: temperately read time from GPP monthly diagnosis files as trial
!      
! The ultimate goal is to update the "restart files"
! will need to modify this subroutine to adapt to restart files.


character(len=*), intent(in) :: filename
type(time_type) :: read_model_time

character(len=*), parameter :: routine = 'read_model_time'
integer, allocatable :: time(:)
integer :: ncid, nt, i, n
integer :: datetime, YYYYMMDD
integer :: year, month, day, hour, minute, second
integer :: yyyy, mm, dd


if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'read_model_time',string1,source)
endif

! Set hour, minute, second = 0 
hour = 0
minute = 0
second = 0




! For now, there's no time dimension in dvmdostem restart file
! In order to preliminarily test our code, we temporarily set
! time equal to 20100101
year = 2010
month = 01
day = 01

!! Below are the code for file that has "time" dimension in netcdf
! open model nc file
! ncid = nc_open_file_readonly(trim(filename), routine)
! ! if we have timeseries
! nt = nc_get_dimension_size(ncid,'time',routine) !get dimension

! ! allocate and read data
! allocate(time(nt))
! call nc_get_variable(ncid, 'time', time, routine)

!============ For multiple time steps =============
! The GPP monthly file has multiple epochtime time 
!n=0  
!do i =1, nt
!   call epochtime_to_date(time(i), 1, datetime, yyyy, mm, dd)
!   if(datetime == assimilation_date)then
!      YYYYMMDD = datetime
!      year = yyyy
!      month = mm
!      day = dd
!      n=1
!      exit
!   endif   
!enddo
!
! ! check if the code correctly find assimilation time in file
! if(n == 0)then
!   write(string1,*) '[read_model_time] NO avaiable time step in the given file'
!   return   
! endif
!
!
!=========  For single time step file ===========
! after read-in "time" in nc file
! call epochtime_to_date(time, 1, datetime, year, month, day)
!
! if(datetime /= assimilation_date)then
!    write(string1,*) '[read_model_time] NO avaiable time step in the given file'
!    return
! endif

!call nc_close_file(ncid)

! Return model time
read_model_time=set_date(year, month, day, hour, minute, second)

end function read_model_time


! dvmdostem (current version) doesn't have time variable in its restart file
! Besides, this write_model_time is not supported for NOLEAP calendar in DART
! We need to build it by our own

subroutine write_model_time(ncid, dart_time)
! This subroutine write out DA time to output file from scratch

integer,             intent(in) :: ncid
type(time_type),     intent(in) :: dart_time

character(len=*), parameter :: routine = 'write_model_time'

integer :: year, month, day, hour, minute, second
integer :: rst_time_yyyymmdd, rst_time_hhmmss
integer :: ymdVarID, hmsVarID
integer :: io, io1

! currently we only output YYYYMMDD
! we can add hourly information later


call get_date(dart_time, year, month, day, hour, minute, second)

rst_time_yyyymmdd = year*10000 + month*100 + day
rst_time_hhmmss  = hour*3600 + minute*60  + second

io1 = nf90_inq_varid(ncid, 'timemgr_rst_curr_ymd', ymdVarID)

! Define time variable if it does not already exist
if (.not. nc_variable_exists(ncid, 'rst_time_yyyymmdd')) then

   call nc_begin_define_mode(ncid)
   
   io = nf90_def_var(ncid,'rst_time_yyyymmdd',NF90_INT,ymdVarID)
   call nc_check(io, routine, 'defining rst_time_yyyymmdd')
   call nc_add_attribute_to_variable(ncid,'rst_time_yyyymmdd','long_name','start date',routine)
   call nc_add_attribute_to_variable(ncid,'rst_time_yyyymmdd','units','YYYYMMDD',routine)
   
   call nc_end_define_mode(ncid)

endif

io = nf90_put_var(ncid, ymdVarID, rst_time_yyyymmdd)
call nc_check(io, routine, 'put_var rst_time_yyyymmdd')

end subroutine write_model_time


subroutine epochtime_to_date(epochtime, run_type, YYYYMMDD, year, month, day)
!==============================================
! convert dvmdostem  epoch time to date time
!
! NOTE: Current dvm-dos-tem uses 365 days calendar for 1 year
!       NO LEAPS YEAR!! 
!
! run_type = 1  'tr' data  
!          = 2  'sc' data
!----------------------------------------------
 integer, intent(in)             :: epochtime
 integer, intent(in)             :: run_type   ! transient('tr') or scenario('sc')          
 integer, intent(out)            :: YYYYMMDD
 integer, intent(out)            :: year, month, day

 character(len=20)               :: date_string
 integer                         :: base_year, days_remaining
 integer, dimension(12)          :: days_in_month  
 logical :: debug = .false.

 ! Define days in months for a common year
  days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


 ! Define the based date 
 if( run_type == 1 )then
    base_year = 1901
 elseif ( run_type == 2 ) then
    base_year = 2016
 else
    write(*,*) "[epochtime_to_date] Can't regonize the run_type!"     
    return
 endif

 ! Convert epochtime to date
  ! year 
  if (epochtime > 365) then
     year = base_year + day / 365
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

  ! convert date to YYYYMMDD
  YYYYMMDD = year*10000 + month*100 + day

  if (debug)then
    WRITE (date_string, '(I4.4,"-",I2.2,"-",I2.2)') year, month, day
    WRITE (*,*) YYYYMMDD
  endif

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

 integer :: i
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


subroutine compute_vertical_value(state_handle, ens_size, location, qty_index, interp_val, istatus)
! Purpose: model_interpolation for vertical variables(no horizontal interpolation)

type(ensemble_type), intent(in)  :: state_handle ! state vector
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
integer,             intent(in)  :: qty_index    ! QTY in DART state needed for interpolation
real(r8),            intent(out) :: interp_val(ens_size)   ! area-weighted result
integer,             intent(out) :: istatus(ens_size)      ! error code (0 == good)

! Local variables
integer     :: varid, imem
integer(i8) :: index1, indexi, indexN
integer     :: xi, yi
real(r8)    :: loc_lat, loc_lon, loc_lev
real(r8)    :: state(ens_size)
real(r8)    :: depthbelow, depthabove
integer     :: levabove, levbelow
real(r8)    :: total(ens_size)
real(r8)    :: total_frac(ens_size)
real(r8), dimension(LocationDims) :: loc
character(len=obstypelength) :: varname

character(len=*), parameter :: routine = 'compute_vertical_value'
logical, parameter          :: debug = .true.



interp_val = MISSING_R8  ! the DART bad value flag
istatus    = 0          

! Get varname
varid = get_varid_from_kind(dom_restart, qty_index)
if (varid < 0) then
   istatus = 1
   return
else
   varname = get_variable_name(dom_restart,varid)
endif

! Get target location
loc        = get_location(location)  ! loc is in DEGREES
loc_lon    = loc(1)
loc_lat    = loc(2)
loc_lev    = loc(3)

call Find_Closet_Indx(loc_lon, loc_lat, xi, yi)


! get index slot
index1 = get_index_start(dom_restart, varname)
indexN = get_index_end(dom_restart, varname)


! Define where the target level located in the vertical column:
call Find_above_below_lev(loc_lev, xi, yi, levabove, levbelow, depthabove, depthbelow)









end subroutine compute_vertical_value


subroutine Find_above_below_lev(loc_lev, xi, yi, levabove, levbelow, depthabove, depthbelow)
real(r8), intent(in)     :: loc_lev
integer,  intent(in)     :: xi, yi                 ! gridcell indx
real(r8), intent(out)    :: depthbelow, depthabove ! depth in meter
integer, intent(out)     :: levabove, levbelow     ! level indx

! local variables
integer                  :: k, nllev
        
! Define where the target level located in the vertical column:
!
!  _______________________________________________________ surface = 0
!            ++++++++++++++                |   |
!                                          |   |
! model lev  ++++++++++++++ levbove (indx) |   |depthabove(m)
!                                          |
!                X --------> loc_lev       |
!                                          |
!            ++++++++++++++ levbelow       |depthbelow (m)
!            ++++++++++++++
! 
! Assumptions:
! 1) If loc_lev(obs lev) fail above surface, interpolate to the surface
! 2) If loc_lev(obs lev) fail below deepest layer, interpolate to the deepest layer
!
! I made these assumptions for simplicity, can modify it if neccessary


! nlev for target grid cell (xi, yi)
nllev = gridnlev(xi,yi)


if (loc_lev  <= 0.0 ) then  ! In TEM, we assumes surface = 0 (top layer)
   depthabove = 0.0
   depthbelow = 0.0
   levabove = 1
   levbelow = 1
elseif (loc_lev >= depth(nllev,xi,yi)) then  ! assume at deepest level
   depthabove    = depth(nllev,xi,yi)         
   depthbelow    = depth(nllev,xi,yi)     
   levabove    = nllev
   levbelow    = nllev
else
   LAYERS : do k = 2, nllev
      if (loc_lev   < depth(k,xi,yi)) then
         levabove = k - 1
         levbelow = k
         depthabove = depth(k-1,xi,yi)
         depthbelow = depth(k,xi,yi)
         exit LAYERS
      endif
   enddo LAYERS
endif


end subroutine 



!---------------------------------------------------------------------------------
! Stolen from DART-CLM... modified to dvm-dos-tem version
! Jan 2024 Note: will be used in PFT patches

subroutine compute_gridcell_value(state_handle, ens_size, location, qty_index, &
                                  interp_val, istatus)

! Purpose: model_interpolation for PFTs variables(no vertical interpolation)

! Passed variables
type(ensemble_type), intent(in)  :: state_handle ! ensemble handler
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
integer,             intent(in)  :: qty_index    ! QTY in DART state needed for interpolation
real(r8),            intent(out) :: interp_val(ens_size)   ! area-weighted result
integer,             intent(out) :: istatus(ens_size)      ! error code (0 == good)

character(len=*), parameter :: routine = 'compute_gridcell_value'

! Local storage
integer     :: varid, imem
integer(i8) :: index1, indexi, indexN
integer     :: xi, yi
real(r8)    :: loc_lat, loc_lon
real(r8)    :: state(ens_size)
real(r8)    :: total(ens_size)
real(r8)    :: total_frac(ens_size)
real(r8), dimension(LocationDims) :: loc
character(len=obstypelength) :: varname

logical, parameter           :: debug = .true.

!-------------------------------------------------------------------
! error message 
! istatus     = 0  Good
! istatus     = 1  Can't find target variable
! istatus     = 2  Errors in weighted-area calculation  

if ( .not. module_initialized ) call static_init_model

interp_val = MISSING_R8  ! the DART bad value flag
istatus    = 0

! This is the target location we want to interpolate to:
loc        = get_location(location)  ! loc is in DEGREES
loc_lon    = loc(1)
loc_lat    = loc(2)

! Next, find the closest grid cell (lon, lat) to the target location
! NOTE: For now, we don't have horizontal interpolation, so here we
!       simply find the closet location
!       while I don't think this is an ideal way to deal with it...

call Find_Closet_Indx(loc_lon, loc_lat, xi, yi)


!! Check if variables are in state
! initialization
varname = "missing_varname"
varid = -1

varid = get_varid_from_kind(dom_restart, qty_index)
if (varid < 0) then
   istatus = 1
   return
else
   varname = get_variable_name(dom_restart,varid)
endif

if(debug)then
   write(*,*) '======= Debugging from [ compute_gridcell_value ] ===== '
   write(*,*)'Working on "',trim(varname),'"'
   write(*,*)'targetlon, lon, loc_lon is ',loc_lon
   write(*,*)'targetlat, lat, loc_lat is ',loc_lat
endif


index1 = get_index_start(dom_restart, varid)
indexN = get_index_end(dom_restart, varid)

total      = 0.0_r8      ! temp storage for state vector
total_frac = 0.0_r8      ! temp storage for area fraction


ELEMENTS : do indexi = index1, indexN

   ! check if it's in the same gridcell
   if ( xindx(indexi) /= xi  ) cycle ELEMENTS
   if ( yindx(indexi) /= yi  ) cycle ELEMENTS
   if ( areafrac_pft(indexi) ==   0.0_r8  ) cycle ELEMENTS

   state = get_state(indexi, state_handle)
!   write(*,*)'Working on index           : ',indexi
!   write(*,*)'state(:) is             : ',state(:)


   MEMBERS: do imem = 1, ens_size

      if(state(imem) == MISSING_R8) cycle MEMBERS
     
      total(imem)      = total(imem)      + state(imem)*areafrac_pft(indexi)
      total_frac(imem) = total_frac(imem) + areafrac_pft(indexi)

      if(debug)then
         write(*,*) '======= Debugging from MEMBERS loops in [ compute_gridcell_value ] ===== '
         write(*,*)'Working on ensemble (imem) : ',imem
         write(*,*)'Working on index           : ',indexi
         write(*,*)'state(imem) is             : ',state(imem)
         write(*,*)'areafrac_pft (indexi)      : ',areafrac_pft(indexi)
      endif

   enddo MEMBERS

!   if(debug)then
!         write(*,*) '======= Debugging from ELEMENTS loops in [ compute_gridcell_value ] ===== '
!         write(*,*)'Working on index           : ',indexi
!         write(*,*)'total_frac(:) is           : ',total_frac(:)
          write(*,*)'total(:)                   : ',total(:)
!    endif


enddo ELEMENTS

do imem = 1,ens_size
   if (total_frac(imem) > 0.0_r8 .and. istatus(imem) == 0) then
      interp_val(imem) = total(imem)
   else
      interp_val(imem) = MISSING_R8
      istatus(imem)    = 2
   endif
enddo

end subroutine compute_gridcell_value


subroutine Find_Closet_Indx(target_lon, target_lat, closest_x, closest_y)
! subroutine to find the closet lon & lat indx 
 real(r8), intent(in)    :: target_lon, target_lat
 integer,  intent(out)   :: closest_x, closest_y ! index for closest point
 integer                 :: i, j
 real(r8)                :: dist, min_dist


min_dist = 1.0E30  ! Initialize with a large number 
! Loop through each point in the arrays
  do i = 1, nx
    do j = 1, ny

      ! Calculate the Euclidean distance to the target location
      dist = sqrt((lon(i, j) - target_lon)**2 + (lat(i, j) - target_lat)**2)

      ! Check if this is the smallest distance so far
      if (dist < min_dist) then
        min_dist = dist
        closest_x = i
        closest_y = j
      end if
    end do
  end do

end subroutine Find_Closet_Indx



!=================================================
! Not neccessary for dvmdostem
!=================================================
subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

end subroutine adv_1step


subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

x = MISSING_R8

end subroutine init_conditions


subroutine init_time(time)

type(time_type), intent(out) :: time

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

