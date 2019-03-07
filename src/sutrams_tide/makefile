# SUTRA-MS example makefile:

SOURCELOC = 
UTILITYLOC = 
NEWMOD = 
PROGRAM = SUTRA_MS.exe
INCDIR= .

# Define the Fortran compiler flags

#F90FLAGS= -O3 -tpp7 -axW -Vaxlib -cm

#F90= ifc

#F90FLAGS= -O2 -fno-second-underscore -r8
#F90FLAGS= -g -fbounds-check -fno-second-underscore -Wall -ftrace=full -r8
#F90FLAGS= -g -fbounds-check -fno-second-underscore -Wall -ftrace=full
#F90FLAGS= -O2 -fno-second-underscore 
#F90= g95a

#F90FLAGS= -O2 -fno-second-underscore -fdefault-real-8 -m64
#F90FLAGS= -O2 -fno-second-underscore 
#F90FLAGS= -g -fbounds-check -fno-second-underscore -Wall -ffpe-trap=invalid,zero 
#F90= gfortran

#release
#F90FLAGS= -O3 -heap-arrays:0 -fpe:0 -traceback
F90FLAGS= -O3 

#debug
#F90FLAGS= -debug:full -traceback -heap-arrays:0

F90= gfortran


# 
# Define the C compile flags
# -D_UF defines UNIX naming conventions for mixed language compilation.
# 
#CFLAGS= -D_UF -O3 -ansi -pedantic
#CFLAGS= -D_UF -O3 
#CC= gcc
CFLAGS= 
CC= cl

# Define the libraries

#SYSLIBS= -lmisalign -ldgc -lm 
#SYSLIBS= -lc
USRLIB  = 

# Define all object files which make up sutrams

SOURCES = \
        SLAPSolver.for \
        MSutraMSPrecision.f90 \
        TecplotIOModule.f90 \
        MFileUnits.f90 \
        MSutraMSErrorHandler.f90 \
        MTotalStorageCalculations.f90 \
        MAutomaticTimeStep.f90 \
        MCommonItems.f90 \
        MUserSpecifiedOutputTime.f90 \
        MSutraHydraulicZones.f90 \
        MSutraStorage.f90 \
        MColumnStorage.f90 \
        MSolverStorage.f90 \
        MTriad2Column.f90 \
        MTecplot.f90 \
        MSpecifiedObservationLocations.f90 \
        MTransientBC.f90 \
        SolverBandedGaussianElimination.f90 \
        NodalProdSorpData.f90 \
        GlobalTriadFormat.f90 \
        Rotate.f90 \
        BoundaryConditionMultipleSpecies.f90 \
        OutputNodeDataToLst_2d.f90 \
        ErrorIO.f90 \
        BandwidthCalculation.f90 \
        SourceMultipleSpecies.f90 \
        InputData2.f90 \
        CloseFiles.f90 \
        ZeroValues.f90 \
        GlobalAn.f90 \
        SutraMSMainProgram.f90 \
        ParseWords.f90 \
        BoundaryCondition1Species.f90 \
        SolverPreparation.f90 \
        Observations.f90 \
        Element_2d.f90 \
        Source.f90 \
        InputData0.f90 \
        BoundaryConditionsInput.f90 \
        Tensor.f90 \
        RotateMatrix.f90 \
        FileOpen.f90 \
        Adsorption.f90 \
        StoreRestartData.f90 \
        OutputNodeDataToLst_3d.f90 \
        BasisCalculations_2d.f90 \
        SLAPSolverWrapper.f90 \
        Nodal.f90 \
        Connectivity.f90 \
        GlobalColumnFormat.f90 \
        BoundaryConditionAssembly.f90 \
        SutraMSSubroutine.f90 \
        PointerSet.f90 \
        Element_3d.f90 \
        Source1Species.f90 \
        Chengji_OutputNodeData.f90 \
		Chengji_AverageVelocity.f90 \
		Chengji_InputET.f90 \
        OutputNodeData.f90 \
        InputData1.f90 \
        AnisotropicDispersion.f90 \
        TriadFormatSetup.f90 \
        SkipCommment.f90 \
        BudgetCalculation.f90 \
        UserDefinedSubroutines.f90 \
        OutputObservationData.f90 \
        Global27NodeMolecule.f90 \
        BasisCalculations_3d.f90 \
        SutraEnd.f90 \
        OutputVelocityData.f90 \
        DimensionWork.f90 \

OBJECTS_1=$(SOURCES:.f90=.o)
OBJECTS=$(OBJECTS_1:.for=.o)

#OBJECTS=$(SOURCES:.f=.o)
## Libraries
#LIBS = -L/usr/lib -llapack -lblas
#LIBS = -L/usr/lib -llapack.so.3 -lblas.so.3
#LIBS = /usr/lib64/atlas/liblapack.so.3 /usr/lib64/libblas.so.3
LIBS  = 

## make all
all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(F90) $(F90FLAGS) $^ -o $@ $(LIBS)

%.o: %.for
	$(F90) -c $(F90FLAGS) $^ 
%.o: %.f90
	$(F90) -c $(F90FLAGS) $^ 


## make clean
clean:
	@rm -f *.o *.mod *~ $(PROGRAM)
	@echo " "
	@echo "cleaning OK"
	@echo "-----------"
