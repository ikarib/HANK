SUBROUTINE dfovec(n,m,x,f)
!objective function for DFLS minimization
!output f vector of least squares objective

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER, INTENT(IN)		:: m,n
REAL(8), INTENT(IN)		:: x(n)
REAL(8), INTENT(OUT)	:: f(m)
REAL(8)					:: lp(12),lweight(m),lf
INTEGER					:: im,ip

objeval = objeval+1
write(4,*) '*********************************'
write(4,*) 'EVALUATION NUMBER: ',objeval

!untransform paramaters and create full parameter vector
ip = 0
IF(EstimateLambda1==1) THEN
	ip = ip+1
	lp(1) = lambdamax*logistic(x(ip))
	write(4,*) ' lambda1: ',lp(1)
ELSE
	lp(1) = lambda1guess
END IF
IF(EstimateZeta1P==1) THEN
	ip = ip+1
	lp(2) = zetamin + (zetamax-zetamin)*logistic(x(ip))			
	write(4,*) ' zeta1P: ',lp(2)
ELSE
	lp(2) = zeta1guess
END IF
IF(EstimateZeta1N==1) THEN
	ip = ip+1
	lp(3) = -(zetamin + (zetamax-zetamin)*logistic(x(ip)))		
	write(4,*) ' zeta1N: ',lp(3)		
ELSE
	lp(3) = -lp(2)
END IF
IF(EstimateLambda2==1) THEN
	ip = ip +1
	lp(4) = lambdamax*logistic(x(ip))	
	write(4,*) ' lambda2: ',lp(4)		
ELSE
! 	lp(4) = 0.0
	lp(4) = lambda2guess	
END IF
IF(EstimateZeta2P==1) THEN
	ip = ip+1
	lp(5) = zetamin + (zetamax-zetamin)*logistic(x(ip))		
	write(4,*) ' zeta2P: ',lp(5)		
ELSE
	lp(5) = zeta2guess
END IF
IF(EstimateZeta2N==1) THEN
	ip = ip+1
	lp(6) = -(zetamin + (zetamax-zetamin)*logistic(x(ip)))		
	write(4,*) ' zeta1N: ',lp(6)
ELSE
	lp(6) = -lp(5)
END IF
IF(EstimateRho1==1) THEN
	ip = ip+1
	lp(7) = logistic(x(ip))
	write(4,*) ' rho1: ',lp(7)		
ELSE
	lp(7) = rho1guess
END IF
IF(EstimateRho2==1) THEN
	ip = ip+1
	lp(8) = logistic(x(ip))
	write(4,*) ' rho2: ',lp(8)		
ELSE
	lp(8) = rho2guess
END IF
IF(EstimateSigma1==1) THEN
	ip = ip+1
	lp(9) = sigmamax*logistic(x(ip))		
	write(4,*) ' sigma1: ',lp(9)		
ELSE
	lp(9) = sigma1guess
END IF
IF(EstimateSigma2==1) THEN
	ip = ip+1
	lp(10) = sigmamax*logistic(x(ip))		
	write(4,*) ' sigma2: ',lp(10)		
ELSE
	lp(10) = sigma2guess
END IF
IF(EstimateDelta1==1) THEN
	ip = ip+1
	lp(11) = deltamax*logistic(x(ip))
	write(4,*) ' delta1: ',lp(11)		
ELSE
	lp(11) = delta1guess
END IF
IF(EstimateDelta2==1) THEN
	ip = ip+1
	lp(12) = deltamax*logistic(x(ip))
	write(4,*) ' delta2: ',lp(12)		
ELSE
	lp(12) = delta2guess
END IF

!extract parameters
lambda1 = lp(1)
zeta1P = lp(2)
zeta1N = lp(3)
lambda2 = lp(4)
zeta2P = lp(5)
zeta2N = lp(6)
rho1 = lp(7)
rho2 = lp(8)
sigma1 = lp(9)
sigma2 = lp(10)
delta1 = lp(11)
delta2 = lp(12)

CALL Simulate
CALL ComputeMoments

lweight = 1.0
im = 0
IF (MatchVarLogY == 1) THEN
	im = im +1
	f(im) = sqrt(lweight(im))*(mu2y/TargetVarLogY -1.0)
	write(4,*) ' VarLogY, target: ',TargetVarLogY, ' model: ',mu2y	
END IF	
IF (MatchVarD1LogY == 1) THEN
	im = im +1
	f(im) = sqrt(lweight(im))*(mu2dy1/TargetVarD1LogY -1.0)
	write(4,*) ' VarD1LogY, target: ',TargetVarD1LogY, ' model: ',mu2dy1
END IF	
IF (MatchSkewD1LogY == 1) THEN
	im = im +1
	f(im) = sqrt(lweight(im))*(gam3dy1/TargetSkewD1LogY -1.0)
	write(4,*) ' SkewD1LogY, target: ',TargetSkewD1LogY, ' model: ',gam3dy1
END IF	
IF (MatchKurtD1LogY == 1) THEN
	im = im +1
	lweight(im) = 0.5
	f(im) = sqrt(lweight(im))*(gam4dy1/TargetKurtD1LogY -1.0)
	write(4,*) ' KurtD1LogY, target: ',TargetKurtD1LogY, ' model: ',gam4dy1
END IF	
IF (MatchVarD5LogY == 1) THEN
	im = im +1
	f(im) = sqrt(lweight(im))*(mu2dy5/TargetVarD5LogY -1.0)
	write(4,*) ' VarD5LogY, target: ',TargetVarD5LogY, ' model: ',mu2dy5
END IF	
IF (MatchSkewD5LogY == 1) THEN
	im = im +1
	lweight(im) = 0.5	
	f(im) = sqrt(lweight(im))*(gam3dy5/TargetSkewD5LogY -1.0)
	write(4,*) ' VarSkewD5LogY, target: ',TargetSkewD5LogY, ' model: ',gam3dy5
END IF	
IF (MatchKurtD5LogY == 1) THEN
	im = im +1
	lweight(im) = 0.5
	f(im) = sqrt(lweight(im))*(gam4dy5/TargetKurtD5LogY -1.0)
	write(4,*) ' KurtD5LogY, target: ',TargetKurtD5LogY, ' model: ',gam4dy5
END IF	
IF (MatchFracD1Less5 == 1) THEN
	im = im +1
	f(im) = sqrt(lweight(im))*(fracdy1less5/TargetFracD1Less5 -1.0)
	write(4,*) ' FracD1Less5, target: ',TargetFracD1Less5, ' model: ',fracdy1less5
END IF	
IF (MatchFracD1Less10 == 1) THEN
	im = im +1
	f(im) = sqrt(lweight(im))*(fracdy1less10/TargetFracD1Less10 -1.0)
	write(4,*) ' FracD1Less10, target: ',TargetFracD1Less10, ' model: ',fracdy1less10
END IF	
IF (MatchFracD1Less20 == 1) THEN
	im = im +1
	f(im) = sqrt(lweight(im))*(fracdy1less20/TargetFracD1Less20 -1.0)
	write(4,*) ' FracD1Less20, target: ',TargetFracD1Less20, ' model: ',fracdy1less20
END IF	
IF (MatchFracD1Less50 == 1) THEN
	im = im +1
	f(im) = sqrt(lweight(im))*(fracdy1less50/TargetFracD1Less50 -1.0)
	write(4,*) ' FracD1Less50, target: ',TargetFracD1Less50, ' model: ',fracdy1less50
END IF	

lweight = lweight/sum(lweight)

lf = sqrt(sum((f**2.0)*lweight)/sum(lweight))
write(4,*) ' objective fun: ',lf
write(4,*) ' '


END SUBROUTINE dfovec