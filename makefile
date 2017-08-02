# SUTRA-MS example makefile:

SOURCELOC = 
UTILITYLOC = 
NEWMOD = 
PROGRAM = ..\bin\SUTRA-MS_1.1_x64.exe
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
F90FLAGS= -O2 -heap-arrays:0 -fpe:0 -traceback
#debug
#F90FLAGS= -debug:full -traceback -heap-arrays:0

F90= ifort


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
        SLAPSolver.obj \
        MSutraMSPrecision.obj \
        TecplotIOModule.obj \
        MFileUnits.obj \
        MSutraMSErrorHandler.obj \
        MTotalStorageCalculations.obj \
        MAutomaticTimeStep.obj \
        MCommonItems.obj \
        MUserSpecifiedOutputTime.obj \
        MSutraHydraulicZones.obj \
        MSutraStorage.obj \
        MColumnStorage.obj \
        MSolverStorage.obj \
        MTriad2Column.obj \
        MTecplot.obj \
        MSpecifiedObservationLocations.obj \
        MTransientBC.obj \
        SolverBandedGaussianElimination.obj \
        NodalProdSorpData.obj \
        GlobalTriadFormat.obj \
        Rotate.obj \
        BoundaryConditionMultipleSpecies.obj \
        OutputNodeDataToLst_2d.obj \
        ErrorIO.obj \
        BandwidthCalculation.obj \
        SourceMultipleSpecies.obj \
        InputData2.obj \
        CloseFiles.obj \
        ZeroValues.obj \
        GlobalAn.obj \
        SutraMSMainProgram.obj \
        ParseWords.obj \
        BoundaryCondition1Species.obj \
        SolverPreparation.obj \
        Observations.obj \
        Element_2d.obj \
        Source.obj \
        InputData0.obj \
        BoundaryConditionsInput.obj \
        Tensor.obj \
        RotateMatrix.obj \
        FileOpen.obj \
        Adsorption.obj \
        StoreRestartData.obj \
        OutputNodeDataToLst_3d.obj \
        BasisCalculations_2d.obj \
        SLAPSolverWrapper.obj \
        Nodal.obj \
        Connectivity.obj \
        GlobalColumnFormat.obj \
        BoundaryConditionAssembly.obj \
        SutraMSSubroutine.obj \
        PointerSet.obj \
        Element_3d.obj \
        Source1Species.obj \
        OutputNodeData.obj \
        InputData1.obj \
        AnisotropicDispersion.obj \
        TriadFormatSetup.obj \
        SkipCommment.obj \
        BudgetCalculation.obj \
        UserDefinedSubroutines.obj \
        OutputObservationData.obj \
        Global27NodeMolecule.obj \
        BasisCalculations_3d.obj \
        SutraEnd.obj \
        OutputVelocityData.obj \
        DimensionWork.obj \

install: sutrams

# Define Task Function Program Modtools

all: sutrams

# Define what sutrams is

sutrams: $(OBJECTS)
	-$(F90) $(F90FLAGS) -o $(PROGRAM) $(OBJECTS) $(USRLIB) $(SYSLIBS)

# Pth_Object codes of sutrams

.f.obj:
	$(F90) $(F90FLAGS) -c $<

.for.obj:
	$(F90) $(F90FLAGS) -c $<

.f90.obj:
	$(F90) $(F90FLAGS) -c $<

.c.obj:
	$(CC) $(CFLAGS) -c $<

clean:
	- del *.obj del *.mod
#
#  end
