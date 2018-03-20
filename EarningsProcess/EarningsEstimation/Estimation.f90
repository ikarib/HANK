SUBROUTINE Estimation

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER	   				:: ip,j,npt
REAL(8), ALLOCATABLE	:: x(:),w(:)


!set up parameters
nparam = 0
IF (EstimateLambda1 == 1) 	nparam = nparam+1
IF (EstimateZeta1P == 1) 	nparam = nparam+1
IF (EstimateZeta1N == 1) 	nparam = nparam+1
IF (EstimateLambda2 == 1) 	nparam = nparam+1
IF (EstimateZeta2P == 1) 	nparam = nparam+1
IF (EstimateZeta2N == 1) 	nparam = nparam+1
IF (EstimateRho1 == 1) 	nparam = nparam+1
IF (EstimateRho2 == 1) 	nparam = nparam+1
IF (EstimateSigma1 == 1) 	nparam = nparam+1
IF (EstimateSigma2 == 1) 	nparam = nparam+1
IF (EstimateDelta1 == 1) 	nparam = nparam+1
IF (EstimateDelta2 == 1) 	nparam = nparam+1

ALLOCATE(paramguess(nparam),paramout(nparam),paramscale(nparam))

!assign guesses
ip = 0
IF(EstimateLambda1==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic(lambda1guess/lambdamax)	
END IF
IF(EstimateZeta1P==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic((zeta1guess-zetamin)/(zetamax-zetamin))	
END IF
IF(EstimateZeta1N==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic((zeta1guess-zetamin)/(zetamax-zetamin))	
END IF
IF(EstimateLambda2==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic(lambda2guess/lambdamax)	
END IF
IF(EstimateZeta2P==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic((zeta2guess-zetamin)/(zetamax-zetamin))		
END IF
IF(EstimateZeta2N==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic((zeta2guess-zetamin)/(zetamax-zetamin))		
END IF
IF(EstimateRho1==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic(rho1guess)
END IF
IF(EstimateRho2==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic(rho2guess)
END IF
IF(EstimateSigma1==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic(sigma1guess/sigmamax)	
END IF
IF(EstimateSigma2==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic(sigma2guess/sigmamax)	
END IF
IF(EstimateDelta1==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic(delta1guess/deltamax)
END IF
IF(EstimateDelta2==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic(delta2guess/deltamax)
END IF

paramscale = 1.0

!set up moments
nmoments = 0
IF (MatchVarLogY == 1)		nmoments = nmoments+1
IF (MatchVarD1LogY == 1)	nmoments = nmoments+1
IF (MatchSkewD1LogY == 1)	nmoments = nmoments+1
IF (MatchKurtD1LogY == 1)	nmoments = nmoments+1
IF (MatchVarD5LogY == 1)	nmoments = nmoments+1
IF (MatchSkewD5LogY == 1)	nmoments = nmoments+1
IF (MatchKurtD5LogY == 1)	nmoments = nmoments+1
IF (MatchFracD1Less5 == 1)	nmoments = nmoments+1
IF (MatchFracD1Less10 == 1)	nmoments = nmoments+1
IF (MatchFracD1Less20 == 1)	nmoments = nmoments+1
IF (MatchFracD1Less50 == 1)	nmoments = nmoments+1


!dfls estimation
ALLOCATE(x(nparam))
x = paramguess*paramscale
npt = 2*nparam+1
IF(ALLOCATED(w)) DEALLOCATE(w)
ALLOCATE(w((npt+11)*(npt+nparam) +nparam*(3*nparam+11)/2) )
maxfun = 500*(nparam+1) 
objeval = 0
OPEN(4, FILE = trim(OutputDir) // 'iterations' //   '.txt', STATUS = 'replace')
DO j = 1,ndfls
	write(4,*) '********************************** '
	write(4,*) 'DFLS/DFBOLS MINIMIZATION ATTEMPT ', j
 	CALL NEWUOA_H(nparam,npt,x,rhobeg,rhoend,iprint,maxfun,w,nmoments)
	write(4,*) '********************************** '
END DO
CLOSE(4)
paramout = x/paramscale

END SUBROUTINE Estimation