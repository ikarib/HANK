SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC)
!Given a function FUNC(X), and given distinct initial points AX and
!BX, this routine searches in the downhill direction (defined by the
!function as evaluated at the initial points) and returns new points
!AX, BX, CX which bracket a minimum of the function. Also returned
!are the function values at the three points, FA, FB and FC.


IMPLICIT NONE

REAL(8),INTENT(INOUT)	:: AX,BX
REAL(8),INTENT(OUT)		:: CX,FA,FB,FC

REAL(8)	::DUM,FU,Q,R,U,ULIM,GOLD,GLIMIT,TINY,FUNC

INTEGER	:: ITER	!added by greg to avoid getting stuck in an infinite loop
PARAMETER(GOLD=1.618034,GLIMIT=100.,TINY=1.D-20)
!The first parameter is the default ratio by which successive intervals
!are magnified; the second is the maximum magnification allowed for
!a parabolic-fit step.

EXTERNAL FUNC

FA=FUNC(AX)
FB=FUNC(BX)
IF(FB.GT.FA) THEN
  DUM=AX
  AX=BX
  BX=DUM
  DUM=FB
  FB=FA
  FA=DUM
ENDIF
CX=BX+GOLD*(BX-AX)
FC=FUNC(CX)
ITER = 1
1 IF((FB.GE.FC) .and. ITER<=25) THEN
	ITER = ITER +1
  R=(BX-AX)*(FB-FC)
  Q=(BX-CX)*(FB-FA)
  U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
  ULIM=BX+GLIMIT*(CX-BX)
  IF((BX-U)*(U-CX).GT.0) THEN
    FU=FUNC(U)
    IF(FU.LT.FC) THEN
	  AX=BX
	  FA=FB
	  BX=U
	  FB=FU
	  GOTO 1
    ELSE IF(FU.GT.FB) THEN
	  CX=U
	  FC=FU
	  GOTO 1
    ENDIF
	U=CX+GOLD*(CX-BX)
	FU=FUNC(U)
  ELSE IF((CX-U)*(U-ULIM).GT.0) THEN
    FU=FUNC(U)
	IF(FU.LT.FC) THEN
	  BX=CX
	  CX=U
	  U=CX+GOLD*(CX-BX)
	  FB=FC
	  FC=FU
	  FU=FUNC(U)
    ENDIF
  ELSE IF((U-ULIM)*(ULIM-CX).GE.0) THEN
    U=ULIM
	FU=FUNC(U)
  ELSE
    U=CX+GOLD*(CX-BX)
	FU=FUNC(U)
  ENDIF
  AX=BX
  BX=CX
  CX=U
  FA=FB
  FB=FC
  FC=FU
  GOTO 1
ENDIF

IF(ITER>=100) THEN
 CX=AX
 FC=FA 
END IF

RETURN
END
