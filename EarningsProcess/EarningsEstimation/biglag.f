
C
C   Important Notice:
C   This GIGLAG are provided in the software NEWUOA, authored by M. J. D. Powell.
C
      SUBROUTINE BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1  KNEW,DELTA,D,ALPHA,GW,HCOL,W)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),
     1  GW(*),HCOL(*),W(*)
C
C     N is the number of variables.
C     NPT is the number of interpolation equations.
C     XOPT is the best interpolation point so far.
C     XPT contains the coordinates of the current interpolation points.
C     BMAT provides the last N columns of H.
C     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     KOPT is the index of the optimal interpolation point.
C     KNEW is the index of the interpolation point that is going to be moved.
C     DELTA is the current trust region bound.
C     D will be set to the step from XOPT to the new point.
C     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
C     GW, HCOL and W will be used for working space.
C
C     The step D is calculated in a way that attempts to maximize the modulus
C     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
C     the KNEW-th Lagrange function.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(ONE)
      DSQ=DELTA*DELTA
      NPTM=NPT-N-1
C
C     Set the first D and GW, where GW is the gradient of LFUNC at XOPT. The
C     first NPT components of HCOL and W will be set to the leading elements
C     of the KNEW-th column of H and the scalar products (XPT(K,.),XOPT),
C     K=1,2,...,NPT, respectively. DGD will be set to the curvature of LFUNC
C     in the direction D.
C
      ITERC=0
      DO 10 I=1,N
      D(I)=XPT(KNEW,I)-XOPT(I)
   10 GW(I)=BMAT(KNEW,I)
      DO 20 K=1,NPT
   20 HCOL(K)=ZERO
      DO 30 J=1,NPTM
      TEMP=ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 30 K=1,NPT
   30 HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
      DGD=ZERO
      DO 50 K=1,NPT
      W(K)=ZERO
      SUM=ZERO
      DO 40 J=1,N
      W(K)=W(K)+XOPT(J)*XPT(K,J)
   40 SUM=SUM+D(J)*XPT(K,J)
      TEMP=HCOL(K)*W(K)
      DGD=DGD+HCOL(K)*SUM*SUM
      DO 50 I=1,N
   50 GW(I)=GW(I)+TEMP*XPT(K,I)
      ALPHA=HCOL(KNEW)
C
C     Step along the direction D or -D if the usefulness of GW is doubtful,
C     where the tests depend on the angle between GW and D and on ||GW||.
C
      DD=ZERO
      GWSQ=ZERO
      DGW=ZERO
      DO 60 I=1,N
      DD=DD+D(I)**2
      DGW=DGW+D(I)*GW(I)
   60 GWSQ=GWSQ+GW(I)**2
      SCALE=DELTA/DSQRT(DD)
      IF (DGW*DGD .LT. ZERO) SCALE=-SCALE
      DO 70 I=1,N
   70 D(I)=SCALE*D(I)
      DGW=SCALE*DGW
      DENOM=DSQ*GWSQ-DGW*DGW
      IF (DENOM .LE. 0.01D0*DSQ*GWSQ) GOTO 150
      VLNEW=DGW+HALF*SCALE*SCALE*DGD
      IF (DSQ*GWSQ .LT. 0.01D0*VLNEW*VLNEW) GOTO 150
C
C     Begin the iteration by making GW orthogonal to D and of length DELTA.
C
   80 ITERC=ITERC+1
      DENOM=DSQRT(DENOM)
      DO 90 I=1,N
   90 GW(I)=(DSQ*GW(I)-DGW*D(I))/DENOM
C
C     Find the elements of W_check, and accumulate their contributions to
C     the coefficients of TAU, which is the restriction of LFUNC to a two
C     dimensional part of the boundary of the trust region.
C
      CF1=ZERO
      CF2=ZERO
      CF3=ZERO
      CF4=ZERO
      CF5=ZERO
      DO 110 K=1,NPT
      TEMPA=ZERO
      TEMPB=ZERO
      DO 100 I=1,N
      TEMPA=TEMPA+XPT(K,I)*D(I)
  100 TEMPB=TEMPB+XPT(K,I)*GW(I)
      TMPA=TEMPA*HCOL(K)
      TMPB=TEMPB*HCOL(K)
      CF1=CF1+HALF*TMPB*TEMPB
      CF2=CF2+TMPA*W(K)
      CF3=CF3+TMPB*W(K)
      CF4=CF4+HALF*(TMPA*TEMPA-TMPB*TEMPB)
  110 CF5=CF5+TMPA*TEMPB
      DO 120 I=1,N
      TEMP=BMAT(KNEW,I)
      CF2=CF2+TEMP*D(I)
  120 CF3=CF3+TEMP*GW(I)
C
C     Seek the value of the angle that maximizes the modulus of TAU.
C
      TAUBEG=CF1+CF2+CF4
      TAUMAX=TAUBEG
      TAUOLD=TAUBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO 130 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
      IF (DABS(TAU) .GT. DABS(TAUMAX)) THEN
          TAUMAX=TAU
          ISAVE=I
          TEMPA=TAUOLD
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=TAU
      END IF
  130 TAUOLD=TAU
      IF (ISAVE .EQ. 0) TEMPA=TAU
      IF (ISAVE .EQ. IU) TEMPB=TAUBEG
      STEP=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-TAUMAX
          TEMPB=TEMPB-TAUMAX
          STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+STEP)
C
C     Calculate the new D and test for convergence.
C
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
      DO 140 I=1,N
  140 D(I)=CTH*D(I)+STH*GW(I)
      IF (ITERC .GE. N) GOTO 200
      IF (ITERC .EQ. 1) TAUBEG=ZERO
      IF (DABS(TAU) .LE. 1.1D0*DABS(TAUBEG)) GOTO 200
C
C     Set GW to the gradient of LFUNC at the new displacement D from XOPT.
C     Then branch for the next iteration unless GW and D are nearly parallel.
C
  150 DO 160 I=1,N
  160 GW(I)=BMAT(KNEW,I)
      DO 180 K=1,NPT
      SUM=W(K)
      DO 170 J=1,N
  170 SUM=SUM+D(J)*XPT(K,J)
      TEMP=HCOL(K)*SUM
      DO 180 I=1,N
  180 GW(I)=GW(I)+TEMP*XPT(K,I)
      GWSQ=ZERO
      DGW=ZERO
      DO 190 I=1,N
      GWSQ=GWSQ+GW(I)**2
  190 DGW=DGW+D(I)*GW(I)
      DENOM=DSQ*GWSQ-DGW*DGW
      IF (DENOM .GE. 1.0D-8*DSQ*GWSQ) GOTO 80
  200 RETURN
      END

