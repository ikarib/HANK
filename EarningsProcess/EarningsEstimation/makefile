release:
	# ENV = /opt/intel/bin/compilervars.sh intel64
	make -f release.makefile

debug:
	# ENV = /opt/intel/bin/compilervars.sh intel64
	make -f debug.makefile
	dsymutil $(OUT).out

hpcrelease:
	make -f hpcrelease.makefile

pbs:
	make -f pbs.makefile

compilesubmit:
	make -f hpcrelease.makefile
	make -f pbs.makefile
	qsub $(OUT).pbs

.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD *genmod* *~ 
	rm -fr *.dSYM