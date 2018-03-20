FC = ifort
FCFLAGS = -m64 -traceback -O3 -qopenmp -implicitnone  -Wl,-stack_size,0x100000000
LDFLAFS = -m64 -traceback -O3 -qopenmp -implicitnone  -Wl,-stack_size,0x100000000

PROG = $(OUT)

MOD = Parameters.o Globals.o Procedures.o random.o

SUBR = 	MakeGuess.o Estimation.o TransitionMatrix.o Simulate.o ComputeMoments.o SetParameters.o dfovec.o newuoa-h.o newuob-h.o update.o trsapp-h.o biglag.o bigden.o InvertMatrix.o SaveOutput.o CombinedProcess.o


OBJ = $(MOD) $(SUBR)

$(PROG).out: $(OBJ) Main.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
Main.o: $(MOD)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
