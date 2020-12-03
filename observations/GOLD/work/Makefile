# Makefile created by mkmf $Id$ 

include ../../../mkmf/mkmf.template


.DEFAULT:
	-touch $@
all: fixsys obs_sequence_tool
fixsys: 
	@ (cd ../../../mpi_utilities; ./fixsystem $(FC) )
utilities_mod.o: ../../../utilities/utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../utilities/utilities_mod.f90
model_mod.o: ../../../models/lorenz_63/model_mod.f90 types_mod.o time_manager_mod.o location_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/lorenz_63/model_mod.f90
cov_cutoff_mod.o: ../../../cov_cutoff/cov_cutoff_mod.f90 types_mod.o utilities_mod.o location_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../cov_cutoff/cov_cutoff_mod.f90
types_mod.o: ../../../common/types_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../common/types_mod.f90
obs_kind_mod.o: ../../../obs_kind/obs_kind_mod.f90 utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_kind/obs_kind_mod.f90
assim_model_mod.o: ../../../assim_model/assim_model_mod.f90 types_mod.o location_mod.o time_manager_mod.o utilities_mod.o model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assim_model/assim_model_mod.f90
null_mpi_utilities_mod.o: ../../../mpi_utilities/null_mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../mpi_utilities/null_mpi_utilities_mod.f90
random_seq_mod.o: ../../../random_seq/random_seq_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../random_seq/random_seq_mod.f90
time_manager_mod.o: ../../../time_manager/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../time_manager/time_manager_mod.f90
location_mod.o: ../../../location/oned/location_mod.f90 types_mod.o utilities_mod.o random_seq_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../location/oned/location_mod.f90
obs_def_mod.o: ../../../obs_def/obs_def_mod.f90 types_mod.o utilities_mod.o location_mod.o assim_model_mod.o obs_kind_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_def/obs_def_mod.f90
obs_sequence_mod.o: ../../../obs_sequence/obs_sequence_mod.f90 types_mod.o location_mod.o obs_def_mod.o obs_kind_mod.o time_manager_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_sequence/obs_sequence_mod.f90
obs_sequence_tool.o: ../../../obs_sequence/obs_sequence_tool.f90 types_mod.o utilities_mod.o location_mod.o obs_def_mod.o obs_kind_mod.o time_manager_mod.o obs_sequence_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_sequence/obs_sequence_tool.f90
./assim_model_mod.f90: ../../../assim_model/assim_model_mod.f90
	cp ../../../assim_model/assim_model_mod.f90 .
./null_mpi_utilities_mod.f90: ../../../mpi_utilities/null_mpi_utilities_mod.f90
	cp ../../../mpi_utilities/null_mpi_utilities_mod.f90 .
./utilities_mod.f90: ../../../utilities/utilities_mod.f90
	cp ../../../utilities/utilities_mod.f90 .
./time_manager_mod.f90: ../../../time_manager/time_manager_mod.f90
	cp ../../../time_manager/time_manager_mod.f90 .
./cov_cutoff_mod.f90: ../../../cov_cutoff/cov_cutoff_mod.f90
	cp ../../../cov_cutoff/cov_cutoff_mod.f90 .
./model_mod.f90: ../../../models/lorenz_63/model_mod.f90
	cp ../../../models/lorenz_63/model_mod.f90 .
./obs_kind_mod.f90: ../../../obs_kind/obs_kind_mod.f90
	cp ../../../obs_kind/obs_kind_mod.f90 .
./obs_sequence_tool.f90: ../../../obs_sequence/obs_sequence_tool.f90
	cp ../../../obs_sequence/obs_sequence_tool.f90 .
./random_seq_mod.f90: ../../../random_seq/random_seq_mod.f90
	cp ../../../random_seq/random_seq_mod.f90 .
./location_mod.f90: ../../../location/oned/location_mod.f90
	cp ../../../location/oned/location_mod.f90 .
./obs_def_mod.f90: ../../../obs_def/obs_def_mod.f90
	cp ../../../obs_def/obs_def_mod.f90 .
./types_mod.f90: ../../../common/types_mod.f90
	cp ../../../common/types_mod.f90 .
./obs_sequence_mod.f90: ../../../obs_sequence/obs_sequence_mod.f90
	cp ../../../obs_sequence/obs_sequence_mod.f90 .
SRC = ../../../utilities/utilities_mod.f90 ../../../models/lorenz_63/model_mod.f90 ../../../cov_cutoff/cov_cutoff_mod.f90 ../../../common/types_mod.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../assim_model/assim_model_mod.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../random_seq/random_seq_mod.f90 ../../../time_manager/time_manager_mod.f90 ../../../location/oned/location_mod.f90 ../../../obs_def/obs_def_mod.f90 ../../../obs_sequence/obs_sequence_mod.f90 ../../../obs_sequence/obs_sequence_tool.f90
OBJ = utilities_mod.o model_mod.o cov_cutoff_mod.o types_mod.o obs_kind_mod.o assim_model_mod.o null_mpi_utilities_mod.o random_seq_mod.o time_manager_mod.o location_mod.o obs_def_mod.o obs_sequence_mod.o obs_sequence_tool.o
OFF = ../../../assim_model/assim_model_mod.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../utilities/utilities_mod.f90 ../../../time_manager/time_manager_mod.f90 ../../../cov_cutoff/cov_cutoff_mod.f90 ../../../models/lorenz_63/model_mod.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../obs_sequence/obs_sequence_tool.f90 ../../../random_seq/random_seq_mod.f90 ../../../location/oned/location_mod.f90 ../../../obs_def/obs_def_mod.f90 ../../../common/types_mod.f90 ../../../obs_sequence/obs_sequence_mod.f90
clean: neat
	-rm -f .cppdefs $(OBJ) obs_sequence_tool
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
obs_sequence_tool: $(OBJ)
	$(LD) $(OBJ) -o obs_sequence_tool  $(LDFLAGS)
