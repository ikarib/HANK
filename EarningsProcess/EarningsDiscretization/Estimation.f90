SUBROUTINE Estimation

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER	   				:: iseed(2),isize,ip,j,npt
REAL(8), ALLOCATABLE	:: x(:),w(:),f(:)

!draw random numbers
isize = 2
iseed(1) = 7755
iseed(2) = 7744
CALL RANDOM_SEED(size = isize)
CALL RANDOM_SEED(put = iseed)
CALL RANDOM_NUMBER(y1rand)
CALL RANDOM_NUMBER(y2rand)

!set up parameters
nparam = 0
IF (EstimateY1width == 1) 	nparam = nparam+1
IF (EstimateY2width == 1) 	nparam = nparam+1
IF (EstimateY1GridPar == 1) 	nparam = nparam+1
IF (EstimateY2GridPar == 1) 	nparam = nparam+1

ALLOCATE(paramguess(nparam),paramout(nparam),paramscale(nparam))

!assign guesses
ip = 0
IF(EstimateY1width==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic((y1widthguess-y1widthmin)/(y1widthmax-y1widthmin))
END IF
IF(EstimateY2width==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic((y2widthguess-y2widthmin)/(y2widthmax-y2widthmin))
END IF
IF(EstimateY1GridPar==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic((y1gridparguess-y1gridmin)/(1.0-y1gridmin))
END IF
IF(EstimateY2GridPar==1) THEN
	ip = ip+1
	paramguess(ip) = invlogistic((y2gridparguess-y2gridmin)/(1.0-y2gridmin))	
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

!
! ALLOCATE(f(nmoments))
! CALL dfovec(nparam,nmoments,x,f)

END SUBROUTINE Estimation