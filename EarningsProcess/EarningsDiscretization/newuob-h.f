C
C   Important Notice:
C   This NEWUOB_H are modifications and based on the subroutine NEWUOB in the software NEWUOA, authored by M. J. D. Powell.
C 
      SUBROUTINE NEWUOB_H(N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,XBASE,
     1  XOPT,XNEW,XPT,GQ,HQ,PQ,BMAT,ZMAT,NDIM,D,VLAG,W,mv)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION X(*),XBASE(*),XOPT(*),XNEW(*),XPT(NPT,*),GQ(*),
     1  HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),VLAG(*),W(*)
 
      integer nmax, mmax, nptmax
      parameter (nmax = 100, mmax=400, nptmax=2*nmax+1)
      dimension GQV(mmax,nmax),HQV(mmax,(nmax+1)*nmax/2),
     &          PQV(mmax,nptmax),
     &          WV(mmax,nmax), 
     &          v_err(mmax), v_beg(mmax), v_temp(mmax),
     &          DIFFV(mmax), v_opt(mmax), v_vquad(mmax),
     &          v_base(mmax), HD1(nmax), GQV_opt(mmax,nmax)
      logical model_update, opt_update, debug

      if (n.gt.nmax) then
        print *, "in newuob_h.f increase the dimension 
     &            nmax to be at least", n
        stop
      endif
      if (mv.gt.mmax) then
        print *, "in newuob_h.f increase the dimension 
     &            mmax to be at least", mv
        stop
      endif   

      debug = .false.
C
C     Set some constants.
C
      model_update = .true.
      opt_update = .true.
      HALF=0.5D0
      ONE=1.0D0
      TENTH=0.1D0
      ZERO=0.0D0
      NP=N+1
      NPTM=NPT-NP
      NFTEST=MAX0(MAXFUN,1)
 
C
C     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
C
      DO 20 J=1,N
      XBASE(J)=X(J)
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 J=1,(N*NP)/2
      do 25 m1=1, mv
   25 HQV(m1,j)=zero
   30 HQ(J)=ZERO
      DO 40 K=1,NPT
      do 35 m1=1, mv
   35 PQV(m1,k)=zero
      PQ(K)=ZERO
      DO 40 J=1,NPTM
   40 ZMAT(K,J)=ZERO
C
C     Begin the initialization procedure. NF becomes one more than the number
C     of function values so far. The coordinates of the displacement of the
C     next initial interpolation point from XBASE are set in XPT(NF,.).
C
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      RECIQ=DSQRT(HALF)/RHOSQ
      NF=0
   50 NFM=NF
      NFMM=NF-N
      NF=NF+1
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              XPT(NF,NFM)=RHOBEG
          ELSE IF (NFM .GT. N) THEN
              XPT(NF,NFMM)=-RHOBEG
          END IF
      END IF
C
C     Calculate the next value of F, label 70 being reached immediately
C     after this calculation. The least function value so far and its index
C     are required.
C
      DO 60 J=1,N
   60 X(J)=XPT(NF,J)+XBASE(J)
      GOTO 310
   70 W(NF)=F
      IF (NF .EQ. 1) THEN
          FBEG=F
          FOPT=F
          do 75 m1=1,mv
          t = v_err(m1)
          v_base(m1) = t
          v_beg(m1) = t
   75     v_opt(m1) = t
          KOPT=1
      ELSE IF (F .LT. FOPT) THEN
          FOPT=F
          do 76 m1=1,mv
   76     v_opt(m1) = v_err(m1)
          KOPT=NF
      END IF
C
C     Set the nonzero initial elements of BMAT and the quadratic model in
C     the cases when NF is at most 2*N+1.
C
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              do 78 m1=1,mv
   78         GQV(m1,NFM)=(v_err(m1)-v_beg(m1))/RHOBEG
              IF (NPT .LT. NF+N) THEN
                  BMAT(1,NFM)=-ONE/RHOBEG
                  BMAT(NF,NFM)=ONE/RHOBEG
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              END IF
          ELSE IF (NFM .GT. N) THEN
              BMAT(NF-N,NFMM)=HALF/RHOBEG
              BMAT(NF,NFMM)=-HALF/RHOBEG
              ZMAT(1,NFMM)=-RECIQ-RECIQ
              ZMAT(NF-N,NFMM)=RECIQ
              ZMAT(NF,NFMM)=RECIQ
              IH=(NFMM*(NFMM+1))/2
              do 79 m1=1, mv
                TEMP=(v_beg(m1)-v_err(m1))/RHOBEG
                HQV(m1,IH)=(GQV(m1,NFMM)-TEMP)/RHOBEG
                GQV(m1,NFMM)=HALF*(GQV(m1,NFMM)+TEMP)
   79         continue              
          END IF

      END IF
      IF (NF .LT. NPT) GOTO 50
C
C     Begin the iterative procedure, because the initial model is complete.
C
      RHO=RHOBEG
      DELTA=RHO
c 
      IDZ=1
      DIFFA=ZERO
      DIFFB=ZERO
      XOPTSQ=ZERO
      DO 80 I=1,N
      XOPT(I)=XPT(KOPT,I)
   80 XOPTSQ=XOPTSQ+XOPT(I)**2
   90 NFSAV=NF
C
C     Generate the next trust region step and test its length. Set KNEW
C     to -1 if the purpose of the next F will be to improve the model.
C
  100 KNEW=0
      if (debug) then
        print *, " Before TRSAPP: delta=",delta, " rho=",rho
      endif
      CALL TRSAPP_H(N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,D,W,W(NP),
     1  W(NP+N),W(NP+2*N),CRVMIN,GQV,HQV,PQV,XBASE,vquad1,
     1  GQV_opt,v_opt,v_base,XOPTSQ,mv,model_update,opt_update)
      DSQ=ZERO
      DO 110 I=1,N
  110 DSQ=DSQ+D(I)**2
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      if (debug) then
       print *, " After TRSAPP: ||d||=",dsqrt(dsq)," vquad1=",vquad1
      endif
  111 IF (DNORM .LT. HALF*RHO) THEN
          KNEW=-1
          DELTA=TENTH*DELTA
          RATIO=-1.0D0
          IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO   
          IF (NF .LE. NFSAV+2) GOTO 460
          TEMP=0.125D0*CRVMIN*RHO*RHO
          IF (TEMP .LE. DMAX1(DIFFA,DIFFB,DIFFC)) GOTO 460
          GOTO 490
      END IF
C
C     Shift XBASE if XOPT may be too far from XBASE. First make the changes
C     to BMAT that do not depend on ZMAT.
C
  120 if (debug) print *, " DSQ=",DSQ, " XOPTSQ=",XOPTSQ
      IF (DSQ .LE. 1.0D-1*XOPTSQ) THEN
          if (debug) print *, " Xbase move"
          model_update = .true.
          TEMPQ=0.25D0*XOPTSQ
          DO 140 K=1,NPT
          SUM=ZERO
          DO 130 I=1,N
  130     SUM=SUM+XPT(K,I)*XOPT(I)
          do 132 m1=1,mv
            v_temp(m1)=PQV(m1,K)*SUM
  132     continue 
          SUM=SUM-HALF*XOPTSQ
          W(NPT+K)=SUM
          DO 140 I=1,N
          do 142 m1=1,mv
            GQV(m1,I)=GQV(m1,I)+v_temp(m1)*XPT(K,I)
  142     continue       
          XPT(K,I)=XPT(K,I)-HALF*XOPT(I)
          VLAG(I)=BMAT(K,I)
          W(I)=SUM*XPT(K,I)+TEMPQ*XOPT(I)
          IP=NPT+I
          DO 140 J=1,I
  140     BMAT(IP,J)=BMAT(IP,J)+VLAG(I)*W(J)+W(I)*VLAG(J)
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO 180 K=1,NPTM
          SUMZ=ZERO
          DO 150 I=1,NPT
          SUMZ=SUMZ+ZMAT(I,K)
  150     W(I)=W(NPT+I)*ZMAT(I,K)
          DO 170 J=1,N
          SUM=TEMPQ*SUMZ*XOPT(J)
          DO 160 I=1,NPT
  160     SUM=SUM+W(I)*XPT(I,J)
          VLAG(J)=SUM
          IF (K .LT. IDZ) SUM=-SUM
          DO 170 I=1,NPT
  170     BMAT(I,J)=BMAT(I,J)+SUM*ZMAT(I,K)
          DO 180 I=1,N
          IP=I+NPT
          TEMP=VLAG(I)
          IF (K .LT. IDZ) TEMP=-TEMP
          DO 180 J=1,I
  180     BMAT(IP,J)=BMAT(IP,J)+TEMP*VLAG(J)
C
C     The following instructions complete the shift of XBASE, including
C     the changes to the parameters of the quadratic model.
C
          IH=0
          DO 200 J=1,N
          do 182 m1=1,mv
            WV(m1,J)=zero
  182     continue
          DO 190 K=1,NPT
          do 192 m1=1,mv
            WV(m1,J)=WV(m1,J)+PQV(m1,K)*XPT(K,J)
  192     continue
  190     XPT(K,J)=XPT(K,J)-HALF*XOPT(J)
          DO 200 I=1,J
          IH=IH+1
          IF (I .LT. J) then
             do 196 m1=1,mv
               GQV(m1,J)=GQV(m1,J)+HQV(m1,IH)*XOPT(I)
  196        continue
          endif
          do 198 m1=1,mv
            GQV(m1,I)=GQV(m1,I)+HQV(m1,IH)*XOPT(J)
            HQV(m1,IH)=HQV(m1,IH)+WV(m1,I)*XOPT(J)
     &                 +XOPT(I)*WV(m1,J)
  198     continue        
  200     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 210 J=1,N
          XBASE(J)=XBASE(J)+XOPT(J)
  210     XOPT(J)=ZERO
          XOPTSQ=ZERO
          do 212 m1=1,mv
            v_base(m1)=v_opt(m1)
  212     continue
      END IF
C
C     Pick the model step if KNEW is positive. A different choice of D
C     may be made later, if the choice of D by BIGLAG causes substantial
C     cancellation in DENOM.
C
      IF (KNEW .GT. 0) THEN          
          CALL BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,KNEW,
     1      DSTEP,D,ALPHA,W,W(NP),W(NDIM+1))
      END IF
C
C     Calculate VLAG and BETA for the current choice of D. The first NPT
C     components of W_check will be held in W.
C
      DO 230 K=1,NPT
      SUMA=ZERO
      SUMB=ZERO
      SUM=ZERO
      DO 220 J=1,N
      SUMA=SUMA+XPT(K,J)*D(J)
      SUMB=SUMB+XPT(K,J)*XOPT(J)
  220 SUM=SUM+BMAT(K,J)*D(J)
      W(K)=SUMA*(HALF*SUMA+SUMB)
  230 VLAG(K)=SUM
      BETA=ZERO
      DO 250 K=1,NPTM
      SUM=ZERO
      DO 240 I=1,NPT
  240 SUM=SUM+ZMAT(I,K)*W(I)
      IF (K .LT. IDZ) THEN
          BETA=BETA+SUM*SUM
          SUM=-SUM
      ELSE
          BETA=BETA-SUM*SUM
      END IF
      DO 250 I=1,NPT
  250 VLAG(I)=VLAG(I)+SUM*ZMAT(I,K)
      BSUM=ZERO
      DX=ZERO
      DO 280 J=1,N
      SUM=ZERO
      DO 260 I=1,NPT
  260 SUM=SUM+W(I)*BMAT(I,J)
      BSUM=BSUM+SUM*D(J)
      JP=NPT+J
      DO 270 K=1,N
  270 SUM=SUM+BMAT(JP,K)*D(K)
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D(J)
  280 DX=DX+D(J)*XOPT(J)
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     If KNEW is positive and if the cancellation in DENOM is unacceptable,
C     then BIGDEN calculates an alternative model step, XNEW being used for
C     working space.
C
      IF (KNEW .GT. 0) THEN
          TEMP=ONE+ALPHA*BETA/VLAG(KNEW)**2
          IF (DABS(TEMP) .LE. 0.8D0) THEN
              CALL BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1          KNEW,DSTEP,D,VLAG,BETA,XNEW,W,W(5*NDIM+1),W)
          END IF
      END IF
C
C     Calculate the next value of the objective function.
C
  290 DO 300 I=1,N
      XNEW(I)=XOPT(I)+D(I)
  300 X(I)=XBASE(I)+XNEW(I)
      NF=NF+1
  310 IF (NF .GT. NFTEST) THEN
          NF=NF-1
          IF (IPRINT .GT. 0) PRINT 320
  320     FORMAT (/4X,'Return from NEWUOA because CALFUN has been',
     1      ' called MAXFUN times.')
          GOTO 530
      END IF
C
C     dfovec(n, mv, x, v_err) provides the values of the vector function v_err(x): R^n \to R^{mv}.
C     Here: n, mv, x \in R^n are input, v_err \in R^{mv} are output.
      call dfovec(n, mv, x, v_err)
C
C     f_value(mv,v_err,F) provides the value of the sum of the squres of the components of v_err(x)
C     i.e. F = sum_{i=1}^{mv} v_err_i (x)^2
C
      call f_value(mv,v_err,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 330, NF,F,(X(I),I=1,N)
  330      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      if (NF.eq.1) then
        iteropt=1
		FBEG=F !added by greg
      else
        if (F.le.FOPT) iteropt = NF
      endif

      if(F.le.dmax1(1.d-12,1.d-20*FBEG)) then
         print *, " F.le.dmax1(1.d-12,1.d-20*FBEG)"
         go to 530
      endif

      IF (NF .LE. NPT) GOTO 70
      IF (KNEW .EQ. -1) GOTO 530
C
C     Use the quadratic model to predict the change in F due to the step D,
C     and set DIFF to the error of this prediction.
C
      do 332 m1=1,mv
        v_vquad(m1)=zero
  332 continue   
      IH=0
      DO 340 J=1,N
      do 342 m1=1,mv
        v_vquad(m1)=v_vquad(m1)+D(J)*GQV(m1,J)
  342 continue    
      DO 340 I=1,J
      IH=IH+1
      TEMP=D(I)*XNEW(J)+D(J)*XOPT(I)
      IF (I .EQ. J) TEMP=HALF*TEMP
      do 340 m1=1,mv
        v_vquad(m1)=v_vquad(m1)+TEMP*HQV(m1,IH)
  340 continue
      DO 345 K=1,NPT
      do 345 m1=1,mv
        v_vquad(m1)=v_vquad(m1)+PQV(m1,K)*W(K)   
  345 continue
      do 350 m1=1,mv
        DIFFV(m1)=v_err(m1)-v_opt(m1)-v_vquad(m1)
  350 continue   
      if (debug) then      
        print *, " Knew=", Knew," vquad1 old=",vquad1
      endif
      if (knew.gt.0) then
        DO 351 I=1,N
  351   HD1(I) = zero
        IH=0
        DO 352 J=1,N
        DO 352 I=1,J
        IH=IH+1
        IF (I .LT. J) HD1(J)=HD1(J)+HQ(IH)*D(I)
  352   HD1(I)=HD1(I)+HQ(IH)*D(J)

        vquad1 = zero
        do 353 i=1,n
  353   vquad1 = vquad1 + D(i)*(GQ(i)+HALF*HD1(i)) 
      endif
      if (debug) then
        print *, " vquad1 new=",vquad1
      endif
      DIFF = F-FOPT-VQUAD1     
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM .GT. RHO) NFSAV=NF
C
C     Update FOPT and XOPT if the new F is the least value of the objective
C     function so far. The branch when KNEW is positive occurs if D is not
C     a trust region step.
C
      FSAVE=FOPT
      IF (F .LT. FOPT) THEN
          opt_update = .true.
          FOPT=F
          do 355 m1=1,mv
  355     v_opt(m1)=v_err(m1)
          XOPTSQ=ZERO
          DO 360 I=1,N
          XOPT(I)=XNEW(I)
  360     XOPTSQ=XOPTSQ+XOPT(I)**2
      END IF
      KSAVE=KNEW
      IF (KNEW .GT. 0) GOTO 410
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (VQUAD1 .GE. ZERO) THEN
          IF (IPRINT .GT. 0) PRINT 370
  370     FORMAT (/4X,'Return from NEWUOA because a trust',
     1      ' region step has failed to reduce Q.')
          GOTO 530
      END IF
      RATIO=(F-FSAVE)/VQUAD1
      if (debug) then
         print *, " Ratio=", ratio
      endif
      IF (RATIO .LE. TENTH) THEN
          DELTA=HALF*DNORM
      ELSE IF (RATIO. LE. 0.7D0) THEN
          DELTA=DMAX1(HALF*DELTA,DNORM)
      ELSE
          DELTA=DMIN1(DMAX1(2.d0*DELTA,4.d0*DNORM),1.d10)
      END IF
      IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
C
C     Set KNEW to the index of the next interpolation point to be deleted.
C
      RHOSQ=DMAX1(TENTH*DELTA,RHO)**2
      KTEMP=0
      DETRAT=ZERO
      IF (F .GE. FSAVE) THEN
          KTEMP=KOPT
          DETRAT=ONE
      END IF
      DO 400 K=1,NPT
      HDIAG=ZERO
      DO 380 J=1,NPTM
      TEMP=ONE
      IF (J .LT. IDZ) TEMP=-ONE
  380 HDIAG=HDIAG+TEMP*ZMAT(K,J)**2
      TEMP=DABS(BETA*HDIAG+VLAG(K)**2)
      DISTSQ=ZERO
      DO 390 J=1,N
  390 DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
      IF (DISTSQ .GT. RHOSQ) TEMP=TEMP*(DISTSQ/RHOSQ)**3
      IF (TEMP .GT. DETRAT .AND. K .NE. KTEMP) THEN
          DETRAT=TEMP
          KNEW=K
      END IF
  400 CONTINUE
      IF (KNEW .EQ. 0) GOTO 460
C
C     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
C     can be moved. Then make this move, and update the quadratic model,
C
  410 CALL UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
      model_update = .true.
      IH=0
      DO 420 I=1,N
      do 422 m1=1,mv
        v_temp(m1)=PQV(m1,KNEW)*XPT(KNEW,I)
  422 continue    
      DO 420 J=1,I
      IH=IH+1
      do 420 m1=1,mv
        HQV(m1,IH)=HQV(m1,IH)+v_temp(m1)*XPT(KNEW,J)
  420 continue   
      do 425 m1=1,mv
        PQV(m1,KNEW)=ZERO
  425 continue  
      DO 440 K=1,NPTM
      IF (ZMAT(KNEW,K) .NE. ZERO) THEN
          do 428 m1=1,mv
            v_temp(m1)=DIFFV(m1)*ZMAT(KNEW,K)
  428     continue     
          IF (K .LT. IDZ) then
             TEMP=-TEMP
             do 429 m1=1,mv
               v_temp(m1)=-v_temp(m1)
  429        continue
          ENDIF
          DO 430 J=1,NPT
          do 430 m1=1,mv
             PQV(m1,J)=PQV(m1,J)+v_temp(m1)*ZMAT(J,K)
  430     continue
      END IF
  440 CONTINUE
      DO 450 I=1,N
      XPT(KNEW,I)=XNEW(I)
      do 450 m1=1,mv
         GQV(m1,I)=GQV(m1,I)+DIFFV(m1)*BMAT(KNEW,I)
  450 continue
      IF (F .LT. FSAVE) KOPT=KNEW
C
C     If a trust region step has provided a sufficient decrease in F, then
C     branch for another trust region calculation. The case KSAVE>0 occurs
C     when the new function value was calculated by a model step.
C
      IF (F .LE. FSAVE+TENTH*VQUAD1) GOTO 100
      IF (KSAVE .GT. 0) GOTO 100
C
C     Alternatively, find out if the interpolation points are close enough
C     to the best point so far.
C
      KNEW=0
  460 DISTSQ=4.0D0*DELTA*DELTA
      DO 480 K=1,NPT
      SUM=ZERO
      DO 470 J=1,N
  470 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM .GT. DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  480 CONTINUE
C
C     If KNEW is positive, then set DSTEP, and branch back for the next
C     iteration, which will generate a "model step".
C
      IF (KNEW .GT. 0) THEN
          DSTEP=DMAX1(DMIN1(TENTH*DSQRT(DISTSQ),HALF*DELTA),RHO)
          DSQ=DSTEP*DSTEP
          GOTO 120
      END IF
      IF (RATIO .GT. ZERO) GOTO 100
c
c Knew =-1, indicating \|d\| \le 1/2 rho
c
      IF (DMAX1(DELTA,DNORM) .GT. RHO) THEN
           IF (KNEW.eq.-1.d0.and.delta.gt.dnorm) then
              KNEW = 0
              go to 111
           Endif 
          GOTO 100
      ENDIF
C
C     The calculations with the current value of RHO are complete. Pick the
C     next values of RHO and DELTA.
C
  490 IF (RHO .GT. RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO .LE. 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO .LE. 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)       
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) PRINT 500
  500         FORMAT (5X)
              PRINT 510, RHO,NF
  510         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              PRINT 520, FOPT,(XBASE(I)+XOPT(I),I=1,N)
  520         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          IF (KNEW.eq.-1.d0.and.delta.gt.dnorm) then
              NFSAV=NF
              KNEW = 0
              go to 111
          Endif 
          GOTO 90
      END IF
C
C     Return from the calculation, after another Newton-Raphson step, if
C     it is too short to have been tried before.
C
      IF (KNEW .EQ. -1) GOTO 290
  530 IF (FOPT .LE. F) THEN
          DO 540 I=1,N
  540     X(I)=XBASE(I)+XOPT(I)
          F=FOPT
      END IF
      IF (IPRINT .GE. 1) THEN
          PRINT 550, NF
  550     FORMAT (/4X,'At the return from NEWUOA',5X,
     1      'Number of function values =',I6)
          PRINT 520, F,(X(I),I=1,N)
          print *, "IterOpt=",iteropt
      END IF
      RETURN
      END

      subroutine f_value(mv,v_err,F)
      integer mv
      double precision v_err(*), F
      integer m1

      F=0.d0
      do 5 m1=1,mv
    5 F = F + v_err(m1)**2

      return
      end
