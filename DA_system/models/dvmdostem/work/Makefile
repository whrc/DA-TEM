# Makefile created by mkmf $Id$ 

include /home/cchang/DART/DA-TEM/DA_system/build_templates/mkmf.template
.DEFAULT:
	-touch $@
all: fixsys model_mod_check
fixsys: 
	@ (cd /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities; ./fixsystem $(FC) )
cov_cutoff_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/cov_cutoff_mod.f90 types_mod.o utilities_mod.o location_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/cov_cutoff_mod.f90
time_manager_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/time_manager_mod.f90
reg_factor_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/reg_factor_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/reg_factor_mod.f90
forward_operator_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/forward_operator_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o mpi_utilities_mod.o obs_sequence_mod.o obs_def_mod.o obs_kind_mod.o ensemble_manager_mod.o distributed_state_mod.o quality_control_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/forward_operator_mod.f90
sampling_error_correction_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/sampling_error_correction_mod.f90 types_mod.o utilities_mod.o netcdf_utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/sampling_error_correction_mod.f90
quad_utils_mod.o: /home/cchang/DART/DA-TEM/DA_system/models/utilities//quad_utils_mod.f90 types_mod.o location_mod.o utilities_mod.o options_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/models/utilities//quad_utils_mod.f90
obs_kind_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_kind_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_kind_mod.f90
state_structure_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_structure_mod.f90 utilities_mod.o obs_kind_mod.o types_mod.o netcdf_utilities_mod.o sort_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_structure_mod.f90
default_location_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//default_location_mod.f90 types_mod.o ensemble_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//default_location_mod.f90
quality_control_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/quality_control_mod.f90 types_mod.o utilities_mod.o location_mod.o obs_sequence_mod.o obs_def_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/quality_control_mod.f90
types_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/types_mod.f90
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/types_mod.f90
sort_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/sort_mod.f90 types_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/sort_mod.f90
io_filenames_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/io_filenames_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o mpi_utilities_mod.o state_structure_mod.o ensemble_manager_mod.o netcdf_utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/io_filenames_mod.f90
obs_impact_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/obs_impact_mod.f90 types_mod.o utilities_mod.o obs_kind_mod.o parse_args_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/obs_impact_mod.f90
location_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/oned/location_mod.f90 types_mod.o utilities_mod.o random_seq_mod.o ensemble_manager_mod.o default_location_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/oned/location_mod.f90
assim_model_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_model_mod.f90 types_mod.o time_manager_mod.o ensemble_manager_mod.o model_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_model_mod.f90
model_check_utilities_mod.o: /home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/model_check_utilities_mod.f90 types_mod.o utilities_mod.o location_mod.o obs_kind_mod.o ensemble_manager_mod.o model_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/model_check_utilities_mod.f90
netcdf_utilities_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/netcdf_utilities_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/netcdf_utilities_mod.f90
filter_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/filter_mod.f90 types_mod.o options_mod.o obs_sequence_mod.o obs_def_mod.o obs_def_utilities_mod.o time_manager_mod.o utilities_mod.o assim_model_mod.o assim_tools_mod.o obs_model_mod.o ensemble_manager_mod.o adaptive_inflate_mod.o mpi_utilities_mod.o random_seq_mod.o state_vector_io_mod.o io_filenames_mod.o direct_netcdf_mod.o state_structure_mod.o forward_operator_mod.o quality_control_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/filter_mod.f90
direct_netcdf_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/direct_netcdf_mod.f90 types_mod.o options_mod.o ensemble_manager_mod.o time_manager_mod.o utilities_mod.o netcdf_utilities_mod.o mpi_utilities_mod.o state_structure_mod.o io_filenames_mod.o assim_model_mod.o model_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/direct_netcdf_mod.f90
tem_to_dart.o: /home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/tem_to_dart.f90 types_mod.o utilities_mod.o netcdf_utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/tem_to_dart.f90
dart_time_io_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/dart_time_io_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o netcdf_utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/dart_time_io_mod.f90
obs_model_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/obs_model_mod.f90 utilities_mod.o assim_model_mod.o state_vector_io_mod.o obs_sequence_mod.o obs_def_mod.o time_manager_mod.o ensemble_manager_mod.o mpi_utilities_mod.o io_filenames_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/obs_model_mod.f90
test_interpolate_oned.o: /home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/test_interpolate_oned.f90 types_mod.o utilities_mod.o netcdf_utilities_mod.o location_mod.o obs_kind_mod.o ensemble_manager_mod.o model_check_utilities_mod.o model_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/test_interpolate_oned.f90
options_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/options_mod.f90 types_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/options_mod.f90
obs_def_utilities_mod.o: /home/cchang/DART/DA-TEM/DA_system//observations/forward_operators/obs_def_utilities_mod.f90 types_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system//observations/forward_operators/obs_def_utilities_mod.f90
obs_sequence_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_sequence_mod.f90 types_mod.o location_mod.o obs_def_mod.o obs_kind_mod.o time_manager_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_sequence_mod.f90
no_cray_win_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/no_cray_win_mod.f90 mpi_utilities_mod.o types_mod.o ensemble_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/no_cray_win_mod.f90
model_mod.o: /home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/model_mod.f90 types_mod.o time_manager_mod.o location_mod.o utilities_mod.o location_io_mod.o netcdf_utilities_mod.o obs_kind_mod.o ensemble_manager_mod.o distributed_state_mod.o state_structure_mod.o default_model_mod.o dart_time_io_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/model_mod.f90
obs_def_mod.o: /home/cchang/DART/DA-TEM/DA_system/observations/forward_operators/obs_def_mod.f90 types_mod.o utilities_mod.o location_mod.o obs_kind_mod.o assim_model_mod.o cov_cutoff_mod.o ensemble_manager_mod.o obs_def_utilities_mod.o time_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/observations/forward_operators/obs_def_mod.f90
random_seq_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/random_seq_mod.f90 types_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/random_seq_mod.f90
location_io_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//location_io_mod.f90 types_mod.o netcdf_utilities_mod.o location_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//location_io_mod.f90
distributed_state_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/distributed_state_mod.f90 mpi_utilities_mod.o types_mod.o ensemble_manager_mod.o no_cray_win_mod.o utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/distributed_state_mod.f90
state_vector_io_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_vector_io_mod.f90 adaptive_inflate_mod.o direct_netcdf_mod.o types_mod.o mpi_utilities_mod.o ensemble_manager_mod.o utilities_mod.o time_manager_mod.o io_filenames_mod.o model_mod.o state_structure_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_vector_io_mod.f90
utilities_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/utilities_mod.f90 types_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/utilities_mod.f90
assert_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/assert_mod.f90
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/assert_mod.f90
adaptive_inflate_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/adaptive_inflate_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o random_seq_mod.o ensemble_manager_mod.o mpi_utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/adaptive_inflate_mod.f90
default_model_mod.o: /home/cchang/DART/DA-TEM/DA_system/models/utilities//default_model_mod.f90 types_mod.o time_manager_mod.o location_mod.o utilities_mod.o netcdf_utilities_mod.o ensemble_manager_mod.o dart_time_io_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/models/utilities//default_model_mod.f90
parse_args_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/parse_args_mod.f90 utilities_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/parse_args_mod.f90
assim_tools_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_tools_mod.f90 types_mod.o options_mod.o utilities_mod.o sort_mod.o random_seq_mod.o obs_sequence_mod.o obs_def_mod.o obs_kind_mod.o cov_cutoff_mod.o reg_factor_mod.o obs_impact_mod.o sampling_error_correction_mod.o location_mod.o ensemble_manager_mod.o mpi_utilities_mod.o adaptive_inflate_mod.o time_manager_mod.o assim_model_mod.o distributed_state_mod.o quality_control_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_tools_mod.f90
mpi_utilities_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/mpi_utilities_mod.f90
ensemble_manager_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/ensemble_manager_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o mpi_utilities_mod.o sort_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/ensemble_manager_mod.f90
model_mod_check.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/programs/model_mod_check/model_mod_check.f90 types_mod.o utilities_mod.o mpi_utilities_mod.o location_mod.o obs_kind_mod.o obs_sequence_mod.o assim_model_mod.o time_manager_mod.o ensemble_manager_mod.o state_vector_io_mod.o state_structure_mod.o io_filenames_mod.o distributed_state_mod.o model_mod.o test_interpolate_oned.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/programs/model_mod_check/model_mod_check.f90
schedule_mod.o: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/schedule_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(MPIFC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/schedule_mod.f90
./netcdf_utilities_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/netcdf_utilities_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/netcdf_utilities_mod.f90 .
./obs_sequence_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_sequence_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_sequence_mod.f90 .
./adaptive_inflate_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/adaptive_inflate_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/adaptive_inflate_mod.f90 .
./mpi_utilities_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/mpi_utilities_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/mpi_utilities_mod.f90 .
./default_model_mod.f90: /home/cchang/DART/DA-TEM/DA_system/models/utilities//default_model_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/models/utilities//default_model_mod.f90 .
./distributed_state_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/distributed_state_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/distributed_state_mod.f90 .
./obs_model_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/obs_model_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/obs_model_mod.f90 .
./filter_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/filter_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/filter_mod.f90 .
./quality_control_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/quality_control_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/quality_control_mod.f90 .
./sampling_error_correction_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/sampling_error_correction_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/sampling_error_correction_mod.f90 .
./no_cray_win_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/no_cray_win_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/no_cray_win_mod.f90 .
./state_vector_io_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_vector_io_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_vector_io_mod.f90 .
./tem_to_dart.f90: /home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/tem_to_dart.f90
	cp /home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/tem_to_dart.f90 .
./parse_args_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/parse_args_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/parse_args_mod.f90 .
./io_filenames_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/io_filenames_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/io_filenames_mod.f90 .
./cov_cutoff_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/cov_cutoff_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/cov_cutoff_mod.f90 .
./reg_factor_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/reg_factor_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/reg_factor_mod.f90 .
./obs_kind_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_kind_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_kind_mod.f90 .
./random_seq_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/random_seq_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/random_seq_mod.f90 .
./time_manager_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/time_manager_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/time_manager_mod.f90 .
./obs_def_utilities_mod.f90: /home/cchang/DART/DA-TEM/DA_system//observations/forward_operators/obs_def_utilities_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system//observations/forward_operators/obs_def_utilities_mod.f90 .
./dart_time_io_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/dart_time_io_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/dart_time_io_mod.f90 .
./utilities_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/utilities_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/utilities_mod.f90 .
./obs_def_mod.f90: /home/cchang/DART/DA-TEM/DA_system/observations/forward_operators/obs_def_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/observations/forward_operators/obs_def_mod.f90 .
./types_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/types_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/types_mod.f90 .
./sort_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/sort_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/sort_mod.f90 .
./location_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/oned/location_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/oned/location_mod.f90 .
./test_interpolate_oned.f90: /home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/test_interpolate_oned.f90
	cp /home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/test_interpolate_oned.f90 .
./model_mod_check.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/programs/model_mod_check/model_mod_check.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/programs/model_mod_check/model_mod_check.f90 .
./model_check_utilities_mod.f90: /home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/model_check_utilities_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/model_check_utilities_mod.f90 .
./ensemble_manager_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/ensemble_manager_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/ensemble_manager_mod.f90 .
./assim_model_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_model_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_model_mod.f90 .
./options_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/options_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/options_mod.f90 .
./forward_operator_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/forward_operator_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/forward_operator_mod.f90 .
./obs_impact_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/obs_impact_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/obs_impact_mod.f90 .
./direct_netcdf_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/direct_netcdf_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/direct_netcdf_mod.f90 .
./location_io_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//location_io_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//location_io_mod.f90 .
./schedule_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/schedule_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/schedule_mod.f90 .
./state_structure_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_structure_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_structure_mod.f90 .
./default_location_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//default_location_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//default_location_mod.f90 .
./assim_tools_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_tools_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_tools_mod.f90 .
./assert_mod.f90: /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/assert_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/assert_mod.f90 .
./quad_utils_mod.f90: /home/cchang/DART/DA-TEM/DA_system/models/utilities//quad_utils_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/models/utilities//quad_utils_mod.f90 .
./model_mod.f90: /home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/model_mod.f90
	cp /home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/model_mod.f90 .
SRC = /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/cov_cutoff_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/time_manager_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/reg_factor_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/forward_operator_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/sampling_error_correction_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/utilities//quad_utils_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_kind_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_structure_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//default_location_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/quality_control_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/types_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/sort_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/io_filenames_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/obs_impact_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/oned/location_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_model_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/model_check_utilities_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/netcdf_utilities_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/filter_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/direct_netcdf_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/tem_to_dart.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/dart_time_io_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/obs_model_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/test_interpolate_oned.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/options_mod.f90 /home/cchang/DART/DA-TEM/DA_system//observations/forward_operators/obs_def_utilities_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_sequence_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/no_cray_win_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/model_mod.f90 /home/cchang/DART/DA-TEM/DA_system/observations/forward_operators/obs_def_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/random_seq_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//location_io_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/distributed_state_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_vector_io_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/utilities_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/assert_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/adaptive_inflate_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/utilities//default_model_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/parse_args_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_tools_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/mpi_utilities_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/ensemble_manager_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/programs/model_mod_check/model_mod_check.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/schedule_mod.f90
OBJ = cov_cutoff_mod.o time_manager_mod.o reg_factor_mod.o forward_operator_mod.o sampling_error_correction_mod.o quad_utils_mod.o obs_kind_mod.o state_structure_mod.o default_location_mod.o quality_control_mod.o types_mod.o sort_mod.o io_filenames_mod.o obs_impact_mod.o location_mod.o assim_model_mod.o model_check_utilities_mod.o netcdf_utilities_mod.o filter_mod.o direct_netcdf_mod.o tem_to_dart.o dart_time_io_mod.o obs_model_mod.o test_interpolate_oned.o options_mod.o obs_def_utilities_mod.o obs_sequence_mod.o no_cray_win_mod.o model_mod.o obs_def_mod.o random_seq_mod.o location_io_mod.o distributed_state_mod.o state_vector_io_mod.o utilities_mod.o assert_mod.o adaptive_inflate_mod.o default_model_mod.o parse_args_mod.o assim_tools_mod.o mpi_utilities_mod.o ensemble_manager_mod.o model_mod_check.o schedule_mod.o
OFF = /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/netcdf_utilities_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_sequence_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/adaptive_inflate_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/mpi_utilities_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/utilities//default_model_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/distributed_state_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/obs_model_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/filter_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/quality_control_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/sampling_error_correction_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/no_cray_win_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_vector_io_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/tem_to_dart.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/parse_args_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/io_filenames_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/cov_cutoff_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/reg_factor_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/obs_kind_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/random_seq_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/time_manager_mod.f90 /home/cchang/DART/DA-TEM/DA_system//observations/forward_operators/obs_def_utilities_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/dart_time_io_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/utilities_mod.f90 /home/cchang/DART/DA-TEM/DA_system/observations/forward_operators/obs_def_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/types_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/sort_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/oned/location_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/test_interpolate_oned.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/programs/model_mod_check/model_mod_check.f90 /home/cchang/DART/DA-TEM/DA_system/models/model_mod_tools/model_check_utilities_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/ensemble_manager_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_model_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/options_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/observations/forward_operator_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/obs_impact_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/direct_netcdf_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//location_io_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/schedule_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/io/state_structure_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/location/utilities//default_location_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/assimilation/assim_tools_mod.f90 /home/cchang/DART/DA-TEM/DA_system/assimilation_code/modules/utilities/assert_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/utilities//quad_utils_mod.f90 /home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/model_mod.f90
clean: neat
	-rm -f .cppdefs *.mod $(OBJ) model_mod_check
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
model_mod_check: $(OBJ)
	$(MPILD) $(OBJ) -o model_mod_check  $(LDFLAGS)
