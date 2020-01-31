FC = ifort
FCFLAGS = -O3 -qopenmp -fpconstant
LDFLAGS = -qopenmp

PROG = $(OUT)

MOD = Parameters.o Globals.o Procedures.o random.o

SUBR = 	MakeGuess.o Estimation.o Simulate.o ComputeMoments.o SetParameters.o dfovec.o newuoa-h.o newuob-h.o update.o trsapp-h.o biglag.o bigden.o InvertMatrix.o SaveOutput.o


OBJ = $(MOD) $(SUBR)

Main: $(OBJ) Main.o
	$(FC) $(LDFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
