FC = ifort
FCFLAGS = -m64 -traceback -O3 -qopenmp -implicitnone  -Wl,-stack_size,0x100000000 -L/Volumes/FILES/Projects/Fortran/SuiteSparse/lib -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas 
LDFLAFS = -m64 -traceback -O3 -qopenmp -implicitnone  -Wl,-stack_size,0x100000000 -L/Volumes/FILES/Projects/Fortran/SuiteSparse/lib -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas 
# -O3

PROG = $(OUT)

MOD = Parameters.o Globals.o umfpack.o Procedures.o 

SUBR = 	AllocateArrays.o SetParameters.o Grids.o IterateBellman.o HJBUpdate.o cumnor.o rtsec.o StationaryDistribution.o SaveSteadyStateOutput.o DistributionStatistics.o rtbis.o rtflsp.o InitialSteadyState.o FinalSteadyState.o SolveSteadyStateEqum.o Calibration.o MomentConditions.o dfovec.o newuoa-h.o newuob-h.o update.o trsapp-h.o biglag.o bigden.o mnbrak.o golden.o sort2.o  CumulativeConsumption.o  FnDiscountRate.o  OptimalConsumption.o FnHoursBC.o  ImpulseResponses.o IRFSequence.o Transition.o  SaveIRFOutput.o IterateTransitionStickyRb.o IterateTransOneAssetStickyRb.o FnCapitalEquity.o CumulativeConsTransition.o DiscountedMPC.o DiscountedMPCTransition.o


OBJ = $(MOD) $(SUBR)

$(PROG).out: $(OBJ) Main.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
Main.o: $(MOD)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
