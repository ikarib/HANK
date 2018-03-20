C
C   Important Notice:
C   This BIGDEN are provided in the software NEWUOA, authored by M. J. D. Powell.
C
      SUBROUTINE BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1  KNEW,DELTA,D,VLAG,BETA,GW,W,WVEC,PROD)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),
     1  VLAG(*),GW(*),W(*),WVEC(NDIM,*),PROD(NDIM,*)
      DIMENSION DEN(9),DENEX(9),WW(9)
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
C     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
C     BETA will be set to the value that will occur in the updating formula
C       when the KNEW-th interpolation point is moved to its new position.
C     GW, W, WVEC, PROD and the private arrays DEN, DENEX and WW will be
C       used for working space, but on return W will be set to W_check.
C       The space for W may be part of the space for PROD.
C
C     D is calculated in a way that should provide a denominator with a large
C     modulus in the updating formula when the KNEW-th interpolation point is
C     shifted to the new position XOPT+D.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      QUART=0.25D0
      TWO=2.0D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(ONE)
      NPTM=NPT-N-1
      DSQ=DELTA*DELTA
C
C     Set the initial D and G, where G is the gradient of the KNEW-th
C     Lagrange function at the trust region centre. The gradient of this
C     function at XPT(KNEW,.) is set in W. The array VLAG will hold the
C     second derivative coefficients of the Lagrange function.
C
      DO 10 I=1,N
      D(I)=XPT(KNEW,I)-XOPT(I)
      GW(I)=BMAT(KNEW,I)
   10 W(I)=BMAT(KNEW,I)
      DO 20 K=1,NPT
   20 VLAG(K)=ZERO
      DO 30 J=1,NPTM
      TEMP=ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 30 K=1,NPT
   30 VLAG(K)=VLAG(K)+TEMP*ZMAT(K,J)
      DO 50 K=1,NPT
      SUM=ZERO
      TEMPB=ZERO
      DO 40 J=1,N
      SUM=SUM+XPT(K,J)*XOPT(J)
   40 TEMPB=TEMPB+XPT(K,J)*XPT(KNEW,J)
      IF (K .EQ. KOPT) XOPTSQ=SUM
      TEMPA=VLAG(K)*SUM
      TEMPB=VLAG(K)*TEMPB
      DO 50 I=1,N
      GW(I)=GW(I)+TEMPA*XPT(K,I)
   50 W(I)=W(I)+TEMPB*XPT(K,I)
      ALPHA=VLAG(KNEW)
C
C     Revise G if its modulus seems to be unusually small.
C
      TEMP=ZERO
      TEMPA=ZERO
      TEMPB=ZERO
      DO 60 I=1,N
      TEMP=TEMP+D(I)**2
      TEMPA=TEMPA+GW(I)**2
   60 TEMPB=TEMPB+W(I)**2
      IF (TEMPB*DSQ .GT. 1.0D4*TEMP*TEMPA) THEN
          DO 70 I=1,N
   70     GW(I)=W(I)
      END IF
C
C     Begin the iteration by making D and G orthogonal and of length DELTA.
C
      ITERC=0
   80 ITERC=ITERC+1
      DD=ZERO
      DG=ZERO
      DO 90 I=1,N
      DD=DD+D(I)**2
   90 DG=DG+D(I)*GW(I)
      GG=ZERO
      DO 100 I=1,N
      GW(I)=DD*GW(I)-DG*D(I)
  100 GG=GG+GW(I)**2
      TEMPD=DELTA/DSQRT(DD)
      TEMPG=DELTA/DSQRT(GG)
      XOPTD=ZERO
      XOPTG=ZERO
      DO 110 I=1,N
      D(I)=TEMPD*D(I)
      GW(I)=TEMPG*GW(I)
      XOPTD=XOPTD+XOPT(I)*D(I)
  110 XOPTG=XOPTG+XOPT(I)*GW(I)
C
C     Set the coefficients of the first two terms of BETA.
C
      TEMPA=HALF*XOPTD*XOPTD
      TEMPB=HALF*XOPTG*XOPTG
      DEN(1)=DSQ*(XOPTSQ+HALF*DSQ)+TEMPA+TEMPB
      DEN(2)=TWO*XOPTD*DSQ
      DEN(3)=TWO*XOPTG*DSQ
      DEN(4)=TEMPA-TEMPB
      DEN(5)=XOPTD*XOPTG
      DO 120 I=6,9
  120 DEN(I)=ZERO
C
C     Put the coefficients of Wcheck in WVEC.
C
      DO 140 K=1,NPT
      TEMPA=ZERO
      TEMPB=ZERO
      TEMPC=ZERO
      DO 130 I=1,N
      TEMPA=TEMPA+XPT(K,I)*D(I)
      TEMPB=TEMPB+XPT(K,I)*GW(I)
  130 TEMPC=TEMPC+XPT(K,I)*XOPT(I)
      WVEC(K,1)=QUART*(TEMPA*TEMPA+TEMPB*TEMPB)
      WVEC(K,2)=TEMPA*TEMPC
      WVEC(K,3)=TEMPB*TEMPC
      WVEC(K,4)=QUART*(TEMPA*TEMPA-TEMPB*TEMPB)
  140 WVEC(K,5)=HALF*TEMPA*TEMPB
      DO 150 I=1,N
      IP=I+NPT
      WVEC(IP,1)=ZERO
      WVEC(IP,2)=D(I)
      WVEC(IP,3)=GW(I)
      WVEC(IP,4)=ZERO
  150 WVEC(IP,5)=ZERO
C
C     Put the coefficents of THETA*Wcheck in PROD.
C
      DO 220 JC=1,5
      KU=NPT
      IF (JC .EQ. 2 .OR. JC .EQ. 3) KU=NDIM
      DO 160 I=1,NPT
  160 PROD(I,JC)=ZERO
      DO 180 J=1,NPTM
      SUM=ZERO
      DO 170 I=1,NPT
  170 SUM=SUM+ZMAT(I,J)*WVEC(I,JC)
      IF (J .LT. IDZ) SUM=-SUM
      DO 180 I=1,NPT
  180 PROD(I,JC)=PROD(I,JC)+SUM*ZMAT(I,J)
      IF (KU .EQ. NDIM) THEN
          DO 200 I=1,NPT
          SUM=ZERO
          DO 190 J=1,N
  190     SUM=SUM+BMAT(I,J)*WVEC(NPT+J,JC)
  200     PROD(I,JC)=PROD(I,JC)+SUM
      END IF
      DO 220 I=1,N
      SUM=ZERO
      DO 210 K=1,KU
  210 SUM=SUM+BMAT(K,I)*WVEC(K,JC)
  220 PROD(NPT+I,JC)=SUM
C
C     Include in DEN the part of BETA that depends on THETA.
C
      DO 240 K=1,NDIM
      SUM=ZERO
      DO 230 I=1,5
      WW(I)=HALF*PROD(K,I)*WVEC(K,I)
  230 SUM=SUM+WW(I)
      DEN(1)=DEN(1)-WW(1)-SUM
      TEMPA=PROD(K,1)*WVEC(K,2)+PROD(K,2)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,4)+PROD(K,4)*WVEC(K,2)
      TEMPC=PROD(K,3)*WVEC(K,5)+PROD(K,5)*WVEC(K,3)
      DEN(2)=DEN(2)-TEMPA-HALF*(TEMPB+TEMPC)
      DEN(6)=DEN(6)-HALF*(TEMPB-TEMPC)
      TEMPA=PROD(K,1)*WVEC(K,3)+PROD(K,3)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,5)+PROD(K,5)*WVEC(K,2)
      TEMPC=PROD(K,3)*WVEC(K,4)+PROD(K,4)*WVEC(K,3)
      DEN(3)=DEN(3)-TEMPA-HALF*(TEMPB-TEMPC)
      DEN(7)=DEN(7)-HALF*(TEMPB+TEMPC)
      TEMPA=PROD(K,1)*WVEC(K,4)+PROD(K,4)*WVEC(K,1)
      DEN(4)=DEN(4)-TEMPA-WW(2)+WW(3)
      TEMPA=PROD(K,1)*WVEC(K,5)+PROD(K,5)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,3)+PROD(K,3)*WVEC(K,2)
      DEN(5)=DEN(5)-TEMPA-HALF*TEMPB
      DEN(8)=DEN(8)-WW(4)+WW(5)
      TEMPA=PROD(K,4)*WVEC(K,5)+PROD(K,5)*WVEC(K,4)
  240 DEN(9)=DEN(9)-HALF*TEMPA
C
C     Extend DEN so that it holds all the coefficients of DENOM.
C
      SUM=ZERO
      DO 250 I=1,5
      WW(I)=HALF*PROD(KNEW,I)**2
  250 SUM=SUM+WW(I)
      DENEX(1)=ALPHA*DEN(1)+WW(1)+SUM
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,2)
      TEMPB=PROD(KNEW,2)*PROD(KNEW,4)
      TEMPC=PROD(KNEW,3)*PROD(KNEW,5)
      DENEX(2)=ALPHA*DEN(2)+TEMPA+TEMPB+TEMPC
      DENEX(6)=ALPHA*DEN(6)+TEMPB-TEMPC
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,3)
      TEMPB=PROD(KNEW,2)*PROD(KNEW,5)
      TEMPC=PROD(KNEW,3)*PROD(KNEW,4)
      DENEX(3)=ALPHA*DEN(3)+TEMPA+TEMPB-TEMPC
      DENEX(7)=ALPHA*DEN(7)+TEMPB+TEMPC
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,4)
      DENEX(4)=ALPHA*DEN(4)+TEMPA+WW(2)-WW(3)
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,5)
      DENEX(5)=ALPHA*DEN(5)+TEMPA+PROD(KNEW,2)*PROD(KNEW,3)
      DENEX(8)=ALPHA*DEN(8)+WW(4)-WW(5)
      DENEX(9)=ALPHA*DEN(9)+PROD(KNEW,4)*PROD(KNEW,5)
C
C     Seek the value of the angle that maximizes the modulus of DENOM.
C
      SUM=DENEX(1)+DENEX(2)+DENEX(4)+DENEX(6)+DENEX(8)
      DENOLD=SUM
      DENMAX=SUM
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      WW(1)=ONE
      DO 280 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      WW(2)=DCOS(ANGLE)
      WW(3)=DSIN(ANGLE)
      DO 260 J=4,8,2
      WW(J)=WW(2)*WW(J-2)-WW(3)*WW(J-1)
  260 WW(J+1)=WW(2)*WW(J-1)+WW(3)*WW(J-2)
      SUMOLD=SUM
      SUM=ZERO
      DO 270 J=1,9
  270 SUM=SUM+DENEX(J)*WW(J)
      IF (DABS(SUM) .GT. DABS(DENMAX)) THEN
          DENMAX=SUM
          ISAVE=I
          TEMPA=SUMOLD
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=SUM
      END IF
  280 CONTINUE
      IF (ISAVE .EQ. 0) TEMPA=SUM
      IF (ISAVE .EQ. IU) TEMPB=DENOLD
      STEP=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-DENMAX
          TEMPB=TEMPB-DENMAX
          STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+STEP)
C
C     Calculate the new D and test for convergence.
C
      WW(2)=DCOS(ANGLE)
      WW(3)=DSIN(ANGLE)
      DO 290 I=1,N
  290 D(I)=WW(2)*D(I)+WW(3)*GW(I)
      DO 300 J=4,8,2
      WW(J)=WW(2)*WW(J-2)-WW(3)*WW(J-1)
  300 WW(J+1)=WW(2)*WW(J-1)+WW(3)*WW(J-2)
      BETA=ZERO
      DENMAX=ZERO
      TAU=ZERO
      DO 310 J=1,9
      BETA=BETA+DEN(J)*WW(J)
      DENMAX=DENMAX+DENEX(J)*WW(J)
  310 IF (J .LE. 5) TAU=TAU+PROD(KNEW,J)*WW(J)
      IF (ITERC .GE. N) GOTO 390
      IF (ITERC .EQ. 1) DENOLD=ZERO
      IF (DABS(DENMAX) .LE. 1.2D0*DABS(DENOLD)) GOTO 390
C
C     Set G to half the gradient of DENOM with respect to D. Then branch
C     for the next iteration.
C
      DO 320 I=1,N
  320 GW(I)=TAU*BMAT(KNEW,I)
      DO 370 K=1,NDIM
      PRVAL=ZERO
      DO 330 J=1,5
  330 PRVAL=PRVAL+PROD(K,J)*WW(J)
      IF (K .LE. NPT) THEN
          SUM=ZERO
          DO 340 I=1,N
  340     SUM=SUM+XPT(K,I)*(XOPT(I)+D(I))
          IF (K .EQ. KOPT) THEN
              TEMPA=ALPHA*(SUM-XOPTSQ+DSQ)
              TEMPB=TEMPA+ALPHA*SUM
              DO 350 I=1,N
  350         GW(I)=GW(I)+TEMPA*XOPT(I)+TEMPB*D(I)
          END IF
          TEMP=(TAU*VLAG(K)-ALPHA*PRVAL)*SUM
          DO 360 I=1,N
  360     GW(I)=GW(I)+TEMP*XPT(K,I)
      ELSE
          GW(K-NPT)=GW(K-NPT)-ALPHA*PRVAL
      END IF
  370 CONTINUE
      GG=ZERO
      DG=ZERO
      DO 380 I=1,N
      GG=GG+GW(I)**2
  380 DG=DG+D(I)*GW(I)
      TEMP=DG*DG/(DSQ*GG)
      IF (TEMP .LE. ONE-1.0D-8) GOTO 80
C
C     Set the vector VLAG before the RETURN from the subroutine.
C
  390 DO 410 K=1,NDIM
      VLAG(K)=ZERO
      SUM=ZERO
      DO 400 J=1,5
      VLAG(K)=VLAG(K)+PROD(K,J)*WW(J)
  400 SUM=SUM+WVEC(K,J)*WW(J)
  410 W(K)=SUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
      RETURN
      END

