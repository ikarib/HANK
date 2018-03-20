ENV = /opt/intel/bin/compilervars.sh intel64

release:
	make -f release.makefile

debug:
	make -f debug.makefile
	dsymutil $(OUT).out

.PHONY: clean 

clean:
	rm -f *.o *.mod *.MOD *genmod*
