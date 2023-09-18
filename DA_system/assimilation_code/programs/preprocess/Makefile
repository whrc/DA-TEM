# Makefile created by mkmf $Id$ 

include /home/cchang/DART/build_templates/mkmf.template
.DEFAULT:
	-touch $@
all: fixsys /home/cchang/DART/assimilation_code/programs/preprocess/preprocess
fixsys: 
	@ (cd /home/cchang/DART/assimilation_code/modules/utilities; ./fixsystem $(FC) )
time_manager_mod.o: /home/cchang/DART/assimilation_code/modules/utilities/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/assimilation_code/modules/utilities/time_manager_mod.f90
parse_args_mod.o: /home/cchang/DART/assimilation_code/modules/utilities/parse_args_mod.f90 utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/assimilation_code/modules/utilities/parse_args_mod.f90
preprocess.o: preprocess.f90 types_mod.o utilities_mod.o parse_args_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	preprocess.f90
null_mpi_utilities_mod.o: /home/cchang/DART/assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
types_mod.o: /home/cchang/DART/assimilation_code/modules/utilities/types_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/assimilation_code/modules/utilities/types_mod.f90
utilities_mod.o: /home/cchang/DART/assimilation_code/modules/utilities/utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	/home/cchang/DART/assimilation_code/modules/utilities/utilities_mod.f90
./null_mpi_utilities_mod.f90: /home/cchang/DART/assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
	cp /home/cchang/DART/assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 .
./utilities_mod.f90: /home/cchang/DART/assimilation_code/modules/utilities/utilities_mod.f90
	cp /home/cchang/DART/assimilation_code/modules/utilities/utilities_mod.f90 .
./types_mod.f90: /home/cchang/DART/assimilation_code/modules/utilities/types_mod.f90
	cp /home/cchang/DART/assimilation_code/modules/utilities/types_mod.f90 .
./time_manager_mod.f90: /home/cchang/DART/assimilation_code/modules/utilities/time_manager_mod.f90
	cp /home/cchang/DART/assimilation_code/modules/utilities/time_manager_mod.f90 .
./parse_args_mod.f90: /home/cchang/DART/assimilation_code/modules/utilities/parse_args_mod.f90
	cp /home/cchang/DART/assimilation_code/modules/utilities/parse_args_mod.f90 .
SRC = /home/cchang/DART/assimilation_code/modules/utilities/time_manager_mod.f90 /home/cchang/DART/assimilation_code/modules/utilities/parse_args_mod.f90 preprocess.f90 /home/cchang/DART/assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 /home/cchang/DART/assimilation_code/modules/utilities/types_mod.f90 /home/cchang/DART/assimilation_code/modules/utilities/utilities_mod.f90
OBJ = time_manager_mod.o parse_args_mod.o preprocess.o null_mpi_utilities_mod.o types_mod.o utilities_mod.o
OFF = /home/cchang/DART/assimilation_code/modules/utilities/null_mpi_utilities_mod.f90 /home/cchang/DART/assimilation_code/modules/utilities/utilities_mod.f90 /home/cchang/DART/assimilation_code/modules/utilities/types_mod.f90 /home/cchang/DART/assimilation_code/modules/utilities/time_manager_mod.f90 /home/cchang/DART/assimilation_code/modules/utilities/parse_args_mod.f90
clean: neat
	-rm -f .cppdefs *.mod $(OBJ) /home/cchang/DART/assimilation_code/programs/preprocess/preprocess
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
/home/cchang/DART/assimilation_code/programs/preprocess/preprocess: $(OBJ)
	$(LD) $(OBJ) -o /home/cchang/DART/assimilation_code/programs/preprocess/preprocess  $(LDFLAGS)
