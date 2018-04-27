FC = ifort
FCFLAGS = -O3 -qopenmp
LDFLAGS = -qopenmp -Wl,-stack_size,0x100000000 -lumfpack

MOD = Parameters.o Globals.o umfpack.o Procedures.o 

SUBR = 	AllocateArrays.o SetParameters.o Grids.o IterateBellman.o HJBUpdate.o cumnor.o rtsec.o StationaryDistribution.o SaveSteadyStateOutput.o DistributionStatistics.o rtbis.o rtflsp.o InitialSteadyState.o FinalSteadyState.o SolveSteadyStateEqum.o Calibration.o MomentConditions.o dfovec.o newuoa-h.o newuob-h.o update.o trsapp-h.o biglag.o bigden.o mnbrak.o golden.o sort2.o  CumulativeConsumption.o  FnDiscountRate.o  OptimalConsumption.o FnHoursBC.o  ImpulseResponses.o IRFSequence.o Transition.o  SaveIRFOutput.o IterateTransitionStickyRb.o IterateTransOneAssetStickyRb.o FnCapitalEquity.o CumulativeConsTransition.o DiscountedMPC.o DiscountedMPCTransition.o


OBJ = $(MOD) $(SUBR)

Main: $(OBJ) Main.o
	$(FC) $(LDFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
