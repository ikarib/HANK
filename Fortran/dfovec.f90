SUBROUTINE dfovec(n,m,x,f)
!objective function for DFLS minimization
!output f vector of least squares objective

USE Parameters
USE Globals

IMPLICIT NONE

INTEGER, INTENT(IN)		:: m,n
REAL(8), INTENT(IN)		:: x(n)
REAL(8), INTENT(OUT)	:: f(m)
REAL(8)					:: lf,lfvec(m)


CALL MomentConditions(n,x,lf,lfvec)
f = lfvec*sqrt(diagweight/real(m))	

END SUBROUTINE dfovec