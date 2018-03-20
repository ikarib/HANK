REAL(8) FUNCTION GOLDEN(AX,BX,CX,F,TOL,XMIN)
!Given a function F, and given a bracketing triplet of abscissas 
!AX, BX, CX (such that BX is between AX and CX, and F(BX) is less 
!than both F(AX) and F(CX)), this routine performs a golden section 
!search for the minimum, isolating it to a fractional precision of 
!about TOL. The abscissa of the minimum is returned as XMIN, and the minimum
!function value is returned as GOLDEN, the returned function value.

IMPLICIT NONE

REAL(8), INTENT(IN)	:: AX,BX,CX,TOL
REAL(8), INTENT(OUT)	:: XMIN

REAL(8)	:: F1,F2,X0,X1,X2,X3,R,C,F
PARAMETER(R=.61803399,C=1.-R)  !Golden ratios

EXTERNAL F


X0=AX  !At any given time we will keep trace of 4 points: 
X3=CX  !X0,X1,X2,X3. 
IF(ABS(CX-BX).GT.ABS(BX-AX)) THEN
	X1=BX
	X2=BX+C*(CX-BX)
ELSE
	X2=BX
	X1=BX-C*(BX-AX)
ENDIF

F1=F(X1)
F2=F(X2)  !Initial function evaluations

1 IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2))) THEN
	IF(F2.LT.F1) THEN
	  X0=X1
	  X1=X2
	  X2=R*X1+C*X3
	  F1=F2
	  F2=F(X2)
	ELSE
	  X3=X2
	  X2=X1
	  X1=R*X2+C*X0
	  F2=F1
	  F1=F(X1)
	ENDIF
GOTO 1
ENDIF

IF(F1.LT.F2) THEN
	GOLDEN=F1
	XMIN=X1
ELSE
	GOLDEN=F2
	XMIN=X2
ENDIF

RETURN
END
