------------------------------
ANALYSIS ON SCF DATA
------------------------------
The directory STATA contains two subfolders. The folder prog contains the DO file SCFgenstats.do and the subfolder data contains the DTA file SCF_04_small.dta. 

The latter is an extract from the Survey of Consumer Finances (SCF) 2004 (with individual ages restricted to ages 22-79) based on the data constructed in 
Kaplan-Violante (ECMA 2014) and Kaplan-Violante-Weidner (BPEA, 2014). Please, contact the authors for additional details. The raw SCF files can be downloaded here:

https://www.federalreserve.gov/econres/scfindex.htm

The file SCFgenstats.do generates the statistics in Table 5 and Table C.1. To run it, you need to specify the path for "rawdir" on line 16. 
While running, the file installs STATA packages "inequal7" and "pshare" needed to compute some of the statistics.


------------------------------
MODEL SOLUTION AND SIMULATION
------------------------------
STEP 1: Fortran Code to Solve Model

The fortran code is contained in Fortran/. It consists of a number of subroutines and functions. The main program is Main.f90 and most parameters are set in SetParameters.f90. The current configuration solves for the steady-state of the model and the transition following the baseline monetary shock in which lump sum transfers adjust. 

The relevant compile line to compile the code depends on the Fortran compiler being used and the system it is being run on. An example makefile is inlcuded in the folder. The file release.makefile is an example compile command that was used with Intel Fortran on a desktop running OSX. The file debug.makefile is an example compile command without optimizations and that creates debug objects. The file hpcrelase.makefile is an example compile command that was used with Intel Fortran on a Unix server.

The code requires BLAS, LAPACK and UMFPACK libraries. These are available in SuiteSparse, which can be downloaded here: http://faculty.cse.tamu.edu/davis/suitesparse.html

The code creates a large number of text files in various subdirectories. You need to specify the path to the desired location for the enclosing output directory in line 13 of SetParameters.f90.

The code uses as input the discretized earnings process. The text files with this information are contained in the folder Fortran/earnings_input/. These are loaded automatically by the Fortran code.

STEP 2: Create Matlab Workspaces from text file output

The matlab file Matlab/MakeWorkspaces.m is then used to create matlab workspaces from the text files that are outputted by the Fortran program. 
In Line 6 you need to specify the path to the directory that contains the directory where the Fortran output is stored
In Line 9 you need to specify the name of the directory that contains the Fortran output

In this way the code can create matlab workspaces corresponding to multiple different sets of output that are contained in separate directories in the same parent directory.

MakeWorkspaces.m will call the functions SteadyStateFigures_fun.m and TransitionFigures_fun.m and will create the followingfiles:
Steadystate_workspace.mat
IRF_Monetary_NOFS_workspace.mat
IRF_Monetary_PE1_workspace.mat
IRF_Monetary_PE2_workspace.mat
etc...
IRF_Monetary_PE15_workspace.mat

The files will be placed in the same directory where the Fortran output was stored. Note that the code will only produce these workspaces if the relevant output files exist in the Fortran output directory.

STEP 3: Run remaining Matlab files to create figures and tables in the paper as follows. All Matlab files are contained in the Matlab/ subdirectory

STEP 3.1: Run MakePlotsPaper.m to make the following figures 1a, 1b, 2a, 2b, 3a, 3b, 4a, 4b, D3a, D3b. The scalars for Table 5 will be output to the screen. You need to change line 6 to the path to the directory where the Fortran output is stored, and line 7 to the path to the directory where you want to store the figures.

STEP 3.2: Run Decomposition.m to get the values for Tables 6 and 7. Each column of these tables refers to a different version of the model. To produce the values for the baseline (Table 6, Column 1), put the path to the Fortran output in line 6. The values for the table are contained in the variable 'finaltable', and will display automatically.

To produce the output for each of the remaining columns of Table 6 and Table 7, you need to re-run the Fortran code with a different set of parameter values, store the output in a folder, and then re-run Decomposition.m with the variable InputDir on line 6 pointing to the relevant folder. For each version of the model, all the changes that need to be made to the Fortran code are contained in the file SetParameters.f90. They are as follows:
 Table 6, Column 2: Set profdistfrac = 1.0 on line 66
 Table 6, Column 3: Set profdistfrac = 0.1 on line 66
 Table 6, Column 4: Set theta = 50.0 on line 179
 Table 6, Column 5: Set phitaylor = 2.0 on line 180
 Table 6, Column 6: Set frisch = 0.5 on line 255 and set rho = 0.0133348071 on line 143
 Table 7, Column 1: Same as baseline
 Table 7, Column 2: Set AdjGovBudgetConstraint = 1 on line 71
 Table 7, Column 3: Set AdjGovBudgetConstraint = 4 on line 71
 Table 7, Column 4: Set AdjGovBudgetConstraint = 3 on line 71
 
STEP 3.3: Run ConDecompDist.m to make the following figures 5a, 5b, 6a, 6b.  You need to change line 7 to the path to the directory where the Fortran output is stored and line 8 to the path to the directory where you want to store the figures.

STEP 3.4: Solve the one asset models. In SetParameters.f90 you need to make the following changes:
 set OneAssetNoCapital = 1 on line 53
 set the discount rate rho to the desired level on line 143
When you have the output you can run the Matlab code PlotOneAssModels.m to make figures 7a and 7b. You need to change line 6 to the path to the enclosing directory where the Fortran output is stored, and line 7 to the path to the directory where you want to store the figures. You need to change lines 10 to 19 to reflect the names of the subdirectories for each value of the discount rate rho.

STEP 3.5 Run PersistencePlots.m to produce figures 8a and 8b. You first need to re-run the Fortran code to produce output under alternative values for the persistence of the monetary shock, in both the B adjust and T adjust case. The output directories should be put on lines 8 to 33 of PersistencePlots.m. When running the Fortran code you set the persistence of the monetary shock by changing the value of MonetaryShockPers on line 80 of SetParameters.f90. For example, to run the case with persistance equal to 0.10 set MonetaryShockPers = exp(-0.10). To run the B adjust case Set AdjGovBudgetConstraint = 3 on line 71 of SetParameters.f90.

STEP 3.6 Run PlotPhillipsCurve.m to produce figures 9a, 9b and 9c. You first need to re-run the Fortran code to produce output under alternative values for the size of the shock to the Taylor rule, in both the B adjust and T adjust case. The output directories should be named *_i for i \in 1...14 where each of the 14 runs corresponds to a difference size shock.  When running the Fortran code you set the size of the monetary shock by changing the value of MonetaryShockSize on line 79 of SetParameters.f90. You need to put the stub for the directory names in lines 9 and 10 of PlotPhillipsCurve.m

------------------------------
EARNINGS PROCESS ESTIMATION
------------------------------

STEP 1: Use the Fortran code in EarningsProcess/EarningsEstimation to estimate the parameters of the income process. The folder contains an example makefile. The output is placed in the directory EarningsProcess/EarningsEstimation/earnings_estimation_output. The actual output used in the paper is included in EarningsProcess/EarningsEstimation/earnings_estimation_paper. These might differ by a small amount due to simulation error. The file parameters.txt contains the parameter estimates. The targeted moments are in SetParameters.f90. The moments of the fitted process are contained in moments.txt. 

STEP 2: Use the Fortran code in EarningsProcess/EarningsDiscretization to discretize the income process. The folder contains an example makefile. The output is placed in the directory EarningsProcess/EarningsDiscretization/earnings_discretization_output. The actual discretization used in the paper is included in EarningsProcess/EarningsEstimation/earnings_discretization_paper. These might differ by a small amount due to simulation error. The file ygrid_combined.txt contains the combined 33 point grid that is used in the main model. The associated continuous time Markov matrix is contained in ymarkov_combined.txt. The associated ergodic distribution is contained in ydist_combined.txt. For convenience, these files are also contained in /Fortran/earnings_input which is where they are read by the full model.

The relevant tables and figures are constructed as follows:
	Table D1: "Model Estimated" is found in EarningsProcess/EarningsEstimation/earnings_estimation_paper/moments.txt
	Table D1: "Model Discretized" is found in EarningsProcess/EarningsDiscretization/earnings_discretization_paper/moments.txt
	Figure D1: Use the Stata file EarningsProcess/EarningsEstimation/PlotDistributions.do to create the figures.
	Figure D2: This figure is created from the file EarningsProcess/EarningsEstimation/earnings_estimation_paper/yannsim_lorenz.txt

