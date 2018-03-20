FC = ifort
# FCFLAGS = -m64 -g -debug all -implicitnone -Wl,-stack_size,0x100000000 -save-temps -warn all -fp-stack-check -ftrapuv -traceback  -L/usr/local/lib -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas -check all
# LDFLAFS = -m64 -g -Wl,-stack_size,0x100000000 -L/usr/local/lib -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas -check all

FCFLAGS = -m64 -g -debug all -implicitnone -Wl,-stack_size,0x100000000 -save-temps -warn all -fp-stack-check -ftrapuv -traceback -check all
LDFLAFS = -m64 -g -Wl,-stack_size,0x100000000 -check all

PROG = $(OUT)

MOD = Parameters.o Globals.o Procedures.o random.o

SUBR = 	MakeGuess.o Estimation.o Simulate.o ComputeMoments.o SetParameters.o dfovec.o newuoa-h.o newuob-h.o update.o trsapp-h.o biglag.o bigden.o InvertMatrix.o SaveOutput.o

OBJ = $(MOD) $(SUBR)

$(PROG).out: $(OBJ) Main.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
Main.o: $(MOD)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

