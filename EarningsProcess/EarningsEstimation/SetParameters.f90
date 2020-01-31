SUBROUTINE SetParameters

USE Parameters
USE Globals
USE random
!USE Procedures

IMPLICIT NONE


INTEGER	   		:: iseed(2),isize
INTEGER			:: in,it

OutputDir = "earnings_estimation_output/"

Display					= 0
SaveSimulations 		= 1
RevertToMedianWithSkew 	= 0
UseNormalDist 			= 1
UseDoubleParetoDist 	= 0
Include2ndProcess 		= 1
AdditiveDrift 			= 0


EstimateLambda1 = 1
EstimateZeta1P 	= 0
EstimateZeta1N 	= 0
EstimateLambda2 = 1
EstimateZeta2P 	= 0
EstimateZeta2N 	= 0
EstimateRho1 	= 0
EstimateRho2 	= 0
EstimateSigma1 	= 1
EstimateSigma2 	= 1
EstimateDelta1 	= 1
EstimateDelta2 	= 1

MatchVarLogY 	= 1
MatchVarD1LogY 	= 1
MatchSkewD1LogY = 0
MatchKurtD1LogY = 1
MatchVarD5LogY  = 1
MatchSkewD5LogY = 0
MatchKurtD5LogY = 1
MatchFracD1Less5 = 0
MatchFracD1Less10 = 1
MatchFracD1Less20 = 1
MatchFracD1Less50 = 1

TargetVarLogY 	 = 0.760
TargetVarD1LogY	 = 0.217
TargetSkewD1LogY = -0.587
TargetKurtD1LogY = 13.377
TargetVarD5LogY  = 0.437
TargetSkewD5LogY = -0.378
TargetKurtD5LogY = 8.782
TargetFracD1Less5 = 0.34
TargetFracD1Less10 = 0.51
TargetFracD1Less20 = 0.68
TargetFracD1Less50 = 0.85

ndfls = 3
rhobeg = 5.0
rhoend = 1.0D-4
iprint = 3
lambdamax = 2.0	!jump process cant arrive more frequently than quarterly on average
sigmamax = 2.0
zetamax = 1000.0
zetamin = 1.0
deltamax = 1.0

CALL system ("mkdir -p " // trim(OutputDir))	

allocate(y1jumprand(nsim,Tsim))
allocate(y2jumprand(nsim,Tsim))
allocate(y1rand(nsim,Tsim))
allocate(y2rand(nsim,Tsim))
allocate(ysim(nsim,Tsim))
allocate(ylevsim(nsim,Tsim))
allocate(y1sim(nsim,Tsim))
allocate(y2sim(nsim,Tsim))
allocate(y1jumpI(nsim,Tsim))
allocate(y2jumpI(nsim,Tsim))
allocate(yannsim(nsim,5))
allocate(yannlevsim(nsim,5))

!set random seed
isize = 2
iseed(1) = 7755
iseed(2) = 7744
CALL RANDOM_SEED(size = isize)
CALL RANDOM_SEED(put = iseed)
CALL RANDOM_NUMBER(y1jumprand)
CALL RANDOM_NUMBER(y2jumprand)
IF (UseNormalDist==1) THEN
	DO in = 1,nsim
	DO it = 1,Tsim
			y1rand(in,it) = random_normal()
			y2rand(in,it) = random_normal()	
	END DO
	END DO
ELSE IF (UseDoubleParetoDist==1) THEN
	CALL RANDOM_NUMBER(y1rand)
	CALL RANDOM_NUMBER(y2rand)
END IF

!OPEN(3, FILE = trim(OutputDir) // 'yjumprand.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nsim,Tsim*2,(/y1jumprand,y2jumprand/));
!OPEN(3, FILE = trim(OutputDir) // 'yrand.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nsim,Tsim*2,(/y1rand,y2rand/));

END SUBROUTINE SetParameters
