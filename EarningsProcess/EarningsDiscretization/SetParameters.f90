SUBROUTINE SetParameters

USE Parameters
USE Globals

IMPLICIT NONE

OutputDir = "earnings_discretization_output/"

Display				= 0
RestrictJumpDomain 	= 0
SaveSimulations 	= 1!0
DriftPointsVersion 	= 2	!2 for the 2 point version, 3 for the 3 point version 
AddPointsCloseZeroY2 = 0

EstimateY1width 	= 1
EstimateY2width 	= 1
EstimateY1GridPar 	= 1
EstimateY2GridPar 	= 1
DefaultY1GridPar1	= 0
DefaultY2GridPar1	= 0

lambda1 = 0.07987
beta1 	= 0.7616
sigma1 	= 1.7352
rho1 	= 0.0

lambda2 = 0.00656
beta2 	= 0.00939
sigma2 	= 1.5259 !1.52585
rho2 	= 0.0

deltaforapprox 	= 1.0

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

TargetVarLogY 	 = 0.7
TargetVarD1LogY	 = 0.23	
TargetSkewD1LogY = -1.35
TargetKurtD1LogY = 17.8
TargetVarD5LogY  = 0.46
TargetSkewD5LogY = -1.01
TargetKurtD5LogY = 11.55
TargetFracD1Less5 = 0.35
TargetFracD1Less10 = 0.54
TargetFracD1Less20 = 0.71
TargetFracD1Less50 = 0.86

ndfls = 4 !6
rhobeg = 5.0 !2.0
rhoend = 1.0D-4	!1.0D-8
iprint = 4


CALL system ("mkdir -p " // trim(OutputDir))	

END SUBROUTINE SetParameters
