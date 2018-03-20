NAME = $(OUT)
SRCDIR="/home/gwk210/hank-main"
RUNDIR="/scratch/gwk210/hank_runs/$(NAME)"
WORK="/work/gwk210/hank_output"

.PHONY: $(NAME).pbs

$(NAME).pbs:

	#create run directory on SCRATCH	
	rm -f $(NAME).pbs
	mkdir -p $(RUNDIR)
	cp $(SRCDIR)/$(NAME).out $(RUNDIR)/
	cp $(SRCDIR)/SetParameters.f90 $(RUNDIR)/
	cp $(SRCDIR)/Parameters.f90 $(RUNDIR)/
	
	#make pbs file
	echo "#! /bin/bash" >> $(NAME).pbs
	
	# set up PBS environement
	echo "#PBS -l nodes=1:ppn=12,walltime=12:00:00" >> $(NAME).pbs
	echo "#PBS -l mem=20GB" >> $(NAME).pbs
	echo "#PBS -M gkaplan@princeton.edu"  >> $(NAME).pbs
	echo "#PBS -m be" >> $(NAME).pbs
	echo "#PBS -j oe"  >> $(NAME).pbs
	
	#run model and zip results
	echo "cd $(RUNDIR)" >> $(NAME).pbs
	echo "./$(NAME).out > output_$(NAME).txt" >> $(NAME).pbs
	echo "cp -R Output/* ./" >> $(NAME).pbs
	echo "rm -r Output" >> $(NAME).pbs
	
	#copy output to directory on WORK
	echo "cd .." >> $(NAME).pbs
	echo "zip -r results_$(NAME).zip $(NAME)" >> $(NAME).pbs
	echo "mkdir -p $(WORK)" >> $(NAME).pbs
	echo "cp results_$(NAME).zip $(WORK)/results_$(NAME).zip" >> $(NAME).pbs
	