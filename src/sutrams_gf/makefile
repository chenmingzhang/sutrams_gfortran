# SUTRA-MS example makefile:

SOURCELOC = 
UTILITYLOC = 
NEWMOD = 
PROGRAM = SUTRA-MS_1.1_x64.exe
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
#F90FLAGS= -O2 -heap-arrays:0 -fpe:0 -traceback
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

OBJECTS = \
        SLAPSolver.o \
        MSutraMSPrecision.o \
        MFileUnits.o \
        TecplotIOModule.o \
        MSutraMSErrorHandler.o \
        MTotalStorageCalculations.o \
        MAutomaticTimeStep.o \
        MCommonItems.o \
        MUserSpecifiedOutputTime.o \
        MSutraHydraulicZones.o \
        MSutraStorage.o \
        MColumnStorage.o \
        MSolverStorage.o \
        MTriad2Column.o \
        MTecplot.o \
        MSpecifiedObservationLocations.o \
        MTransientBC.o \
        SolverBandedGaussianElimination.o \
        NodalProdSorpData.o \
        GlobalTriadFormat.o \
        Rotate.o \
        BoundaryConditionMultipleSpecies.o \
        OutputNodeDataToLst_2d.o \
        ErrorIO.o \
        BandwidthCalculation.o \
        SourceMultipleSpecies.o \
        InputData2.o \
        CloseFiles.o \
        ZeroValues.o \
        GlobalAn.o \
        SutraMSMainProgram.o \
        ParseWords.o \
        BoundaryCondition1Species.o \
        SolverPreparation.o \
        Observations.o \
        Element_2d.o \
        Source.o \
        InputData0.o \
        BoundaryConditionsInput.o \
        Tensor.o \
        RotateMatrix.o \
        FileOpen.o \
        Adsorption.o \
        StoreRestartData.o \
        OutputNodeDataToLst_3d.o \
        BasisCalculations_2d.o \
        SLAPSolverWrapper.o \
        Nodal.o \
        Connectivity.o \
        GlobalColumnFormat.o \
        BoundaryConditionAssembly.o \
        SutraMSSubroutine.o \
        PointerSet.o \
        Element_3d.o \
        Source1Species.o \
        OutputNodeData.o \
        InputData1.o \
        AnisotropicDispersion.o \
        TriadFormatSetup.o \
        SkipCommment.o \
        BudgetCalculation.o \
        UserDefinedSubroutines.o \
        OutputObservationData.o \
        Global27NodeMolecule.o \
        BasisCalculations_3d.o \
        SutraEnd.o \
        OutputVelocityData.o \
        DimensionWork.o \

install: sutrams

# Define Task Function Program Modtools

all: sutrams

# Define what sutrams is

sutrams: $(OBJECTS)
	-$(F90) $(F90FLAGS) -o $(PROGRAM) $(OBJECTS) $(USRLIB) $(SYSLIBS)

# Pth_Object codes of sutrams

.f.o:
	$(F90) $(F90FLAGS) -c $<

.for.o:
	$(F90) $(F90FLAGS) -c $<

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.F90.o:
	$(F90) $(F90FLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	- del *.o del *.mod
#
#  end
