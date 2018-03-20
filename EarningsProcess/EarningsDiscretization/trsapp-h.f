C
C   Important Notice:
C   This TRSAPP_H are modifications and based on the subroutine TRSAPP in the software NEWUOA, authored by M. J. D. Powell.
C
      SUBROUTINE TRSAPP_H(N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,STEP,
     1  D,G,HD,HS,CRVMIN,GQV,HQV,PQV,XBASE,vquad,
     1  GQV_opt,v_opt,v_base,XOPTSQ,mv,model_update,opt_update)
      IMPLICIT double precision (A-H,O-Z)
      parameter (nmax = 100, mmax=400)
      DIMENSION XOPT(*),XPT(NPT,*),GQ(*),HQ(*),PQ(*),STEP(*),
     1  D(*),G(*),HD(*),HS(*),GQV(mmax,*),HQV(mmax,*),PQV(mmax,*),
     1  XBASE(*),v_opt(*),GQV_opt(mmax,*),v_base(*)
      logical model_update, opt_update, zero_res, debug
      double precision v_gtemp(mmax), gbeg(nmax), gtemp(nmax)  

      debug = .false.
      HALF=0.5D0
      ZERO=0.0D0

      if (n.gt.nmax) then
        print *, "in trsapp_h.f increase the dimension 
     &            nmax to be at least", n
        stop
      endif
      if (mv.gt.mmax) then
        print *, "in trsapp_h.f increase the dimension 
     &            mmax to be at least", mv
        stop
      endif 

      if ((.not.model_update).and.(.not.opt_update)) go to 8
       model_update = .false.
       opt_update = .false.


      if (dsqrt(XOPTSQ).gt.0.25d0*delta) then
c
c Use the gradient at xopt to formulate J^t J
c
       do m1=1,mv
         do i=1,n
           GQV_opt(m1,i) = GQV(m1,i)
         enddo
         do k=1, npt
           temp = zero
           do j=1,n
             temp = temp+XPT(k,j)*XOPT(j)
           enddo
           temp = temp*PQV(m1,k)
           do i=1,n
             GQV_opt(m1,i)=GQV_opt(m1,i)+temp*XPT(k,i)
           enddo
         enddo
         IH=0
         do j=1,n
           do i=1,j
             IH=IH+1
             if (i .lt. j) GQV_opt(m1,j)=GQV_opt(m1,j)
     &                                   +HQV(m1,IH)*XOPT(i)
                GQV_opt(m1,i)=GQV_opt(m1,i)+HQV(m1,IH)*XOPT(j)
           enddo
         enddo  
       enddo

       call f_grad(mv,v_opt,v_gtemp)
       gnorm2 = zero
       do i=1,n
         GQ(i) = zero
         do m1=1,mv
            GQ(i) = GQ(i) + v_gtemp(m1)*GQV_opt(m1,i)
         enddo
         gnorm2 = gnorm2 + GQ(i)**2
       enddo
c
c Calculate the explicite Hessian.
c
       f_opt = zero
       call f_value(mv,v_opt,f_opt)
       if (gnorm2.ge.1.d0.or.f_opt.le.dsqrt(gnorm2)) then
         zero_res = .true.
       else
         zero_res = .false.
       endif

       IH=0
       do j=1,n
         do i=1,j
           IH=IH+1
           if (zero_res) then
             t1 = zero
             do m1=1,mv  
               t1 = t1+GQV_opt(m1,i)*GQV_opt(m1,j)
             enddo
             HQ(IH) = 2.d0*t1
           else
             t1 = zero
             do m1=1,mv
               t2 = zero
               do k=1,npt
                 t2 = t2 + XPT(k,i)*PQV(m1,k)*XPT(k,j)
               enddo
               t2 = t2 + HQV(m1,IH)
               t1 = t1+(GQV_opt(m1,i)*GQV_opt(m1,j)+v_opt(m1)*t2) 
             enddo
             HQ(IH) = 2.d0*t1
           endif
         enddo
       enddo

      else
c
c Use the gradient at xbase to formulate J^t J
c
       call f_grad(mv,v_base,v_gtemp)
       gnorm2 = zero
       do i=1,n
         GQ(i) = zero
         do m1=1,mv
            GQ(i) = GQ(i) + v_gtemp(m1)*GQV(m1,i)
         enddo
         gnorm2 = gnorm2 + GQ(i)**2
       enddo
c
c Calculate the explicite Hessian.
c
       f_base = zero
       call f_value(mv,v_base,f_base)
       if (gnorm2.ge.1.d0.or.f_base.le.dsqrt(gnorm2)) then
         zero_res = .true.
       else
         zero_res = .false.
       endif

       IH=0
       do j=1,n
         do i=1,j
           IH=IH+1
           if (zero_res) then
             t1 = zero
             do m1=1,mv           
               t1 = t1+GQV(m1,i)*GQV(m1,j)
             enddo
             HQ(IH) = 2.d0*t1
           else           
             t1 = zero
             do m1=1,mv
               t2 = zero
               do k=1,npt
                 t2 = t2 + XPT(k,i)*PQV(m1,k)*XPT(k,j)
               enddo
               t2 = t2 + HQV(m1,IH)
               t1 = t1+(GQV(m1,i)*GQV(m1,j)+v_base(m1)*t2) 
             enddo
             HQ(IH) = 2.d0*t1
           endif
         enddo
       enddo
c calculte the gradient at xopt
       IH=0
       do j=1,n
         do i=1,j
           IH=IH+1
           if (i .lt. j) GQ(j)=GQ(j)+HQ(IH)*XOPT(i)
              GQ(i)=GQ(i)+HQ(IH)*XOPT(j)
         enddo
       enddo  

      endif

    8 HALF=0.5D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(1.0D0)
      DELSQ=DELTA*DELTA
      ITERC=0
      ITERMAX=N
      ITERSW=ITERMAX
      if (debug) then
        t = zero
        do i=1,n
          t = t + xopt(i)**2
        enddo
        print *, " ||xopt||=",dsqrt(t)
      endif
      gnorm2 = zero
      DO 10 I=1,N
      gnorm2 = gnorm2 + GQ(i)**2
   10 D(I)=zero 
      gnorm2 = dsqrt(gnorm2)  
      if (debug) print *, " gnorm2=", gnorm2 
      GOTO 170
C
C     Prepare for the first line search.
C
   20 continue
      QRED=ZERO
      DD=ZERO
      DO 30 I=1,N
      STEP(I)=ZERO
      HS(I)=ZERO
      G(I) = GQ(I)
      D(I)=-G(I)
      gbeg(i) = G(i)
   30 DD=DD+D(I)**2
      CRVMIN=ZERO
      IF (DD .EQ. ZERO) GOTO 160
      DS=ZERO
      SS=ZERO
      GG=DD
      GGBEG=GG
      if (debug) print *, " GGBEG=", GGBEG
C
C     Calculate the step to the trust region boundary and the product HD.
C
   40 ITERC=ITERC+1
      TEMP=DELSQ-SS
      BSTEP=TEMP/(DS+DSQRT(DS*DS+DD*TEMP))
c      BSTEP=(-DS+DSQRT(DS*DS+DD*TEMP))/DD
      if (debug) print *, " BSTEP=", BSTEP
      GOTO 170
   50 DHD=ZERO
      DO 60 J=1,N
   60 DHD=DHD+D(J)*HD(J)
C
C     Update CRVMIN and set the step-length ALPHA.
C
      ALPHA=BSTEP
      if (debug) then
         print *, " ITERC=",ITERC
         print *, " DHD/DD=", DHD/DD
      endif
      IF (DHD .GT. ZERO) THEN
          TEMP=DHD/DD
          IF (ITERC .EQ. 1) CRVMIN=TEMP
          CRVMIN=DMIN1(CRVMIN,TEMP)
          ALPHA=DMIN1(ALPHA,GG/DHD)
      END IF
      QADD=ALPHA*(GG-HALF*ALPHA*DHD)
      QRED=QRED+QADD
C
C     Update STEP and HS.
C
      GGSAV=GG
      GG=ZERO
      DO 70 I=1,N
      STEP(I)=STEP(I)+ALPHA*D(I)
      HS(I)=HS(I)+ALPHA*HD(I)
   70 GG=GG+(G(I)+HS(I))**2
      if (debug) print *, " GG=",GG

      IF (GG .LE. dmin1(1.0D-4*GGBEG,1.d-16)) GOTO 160 
      IF (GG .LE. 1.d-14*gnorm2) GOTO 160   

      IF (ITERC .EQ. ITERMAX) GOTO 160
C
C     Begin another conjugate direction iteration if required.
C
      IF (ALPHA .LT. BSTEP) THEN
          IF (QADD .LE. 1.D-6*QRED) GOTO 160
          TEMP=GG/GGSAV
          DD=ZERO
          DS=ZERO
          SS=ZERO
          DO 80 I=1,N
          D(I)=TEMP*D(I)-G(I)-HS(I)
          DD=DD+D(I)**2
          DS=DS+D(I)*STEP(I)
   80     SS=SS+STEP(I)**2
          IF (SS .LT. DELSQ) GOTO 40
      END IF
      CRVMIN=ZERO
      ITERSW=ITERC
C
C     Test whether an alternative iteration is required.
C
   90 IF (GG .LE. 1.0D-4*GGBEG) GOTO 160
      if (debug) print *, "curve search performed"
      SG=ZERO
      SHS=ZERO
      DO 100 I=1,N
      SG=SG+STEP(I)*G(I)
  100 SHS=SHS+STEP(I)*HS(I)
      SGK=SG+SHS
      ANGTEST=SGK/DSQRT(GG*DELSQ)
      IF (ANGTEST .LE. -0.99D0) GOTO 160
C
C     Begin the alternative iteration by calculating D and HD and some
C     scalar products.
C
      ITERC=ITERC+1
      TEMP=DSQRT(DELSQ*GG-SGK*SGK)
      TEMPA=DELSQ/TEMP
      TEMPB=SGK/TEMP
      DO 110 I=1,N
  110 D(I)=TEMPA*(G(I)+HS(I))-TEMPB*STEP(I)
      GOTO 170
  120 DG=ZERO
      DHD=ZERO
      DHS=ZERO
      DO 130 I=1,N
      DG=DG+D(I)*G(I)
      DHD=DHD+HD(I)*D(I)
  130 DHS=DHS+HD(I)*STEP(I)
C
C     Seek the value of the angle that minimizes Q.
C
      CF=HALF*(SHS-DHD)
      QBEG=SG+CF
      QSAV=QBEG
      QMIN=QBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO 140 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      QNEW=(SG+CF*CTH)*CTH+(DG+DHS*CTH)*STH
      IF (QNEW .LT. QMIN) THEN
          QMIN=QNEW
          ISAVE=I
          TEMPA=QSAV
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=QNEW
      END IF
  140 QSAV=QNEW
      IF (ISAVE .EQ. ZERO) TEMPA=QNEW
      IF (ISAVE .EQ. IU) TEMPB=QBEG
      ANGLE=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-QMIN
          TEMPB=TEMPB-QMIN
          ANGLE=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+ANGLE)
C
C     Calculate the new STEP and HS. Then test for convergence.
C
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      REDUC=QBEG-(SG+CF*CTH)*CTH-(DG+DHS*CTH)*STH
      GG=ZERO
      DO 150 I=1,N
      STEP(I)=CTH*STEP(I)+STH*D(I)
      HS(I)=CTH*HS(I)+STH*HD(I)
  150 GG=GG+(G(I)+HS(I))**2
      QRED=QRED+REDUC
      RATIO=REDUC/QRED
      IF (ITERC .LT. ITERMAX .AND. RATIO .GT. 0.01D0) GOTO 90
  160 continue
      do 161 i=1,n
         HD(i) = zero
  161 continue
      IH=0
      DO 162 J=1,N
      DO 162 I=1,J
      IH=IH+1
      IF (I .LT. J) HD(J)=HD(J)+HQ(IH)*step(I)
  162 HD(I)=HD(I)+HQ(IH)*step(J)
c      vquad = zero
c      do 163 i=1,n
c  163 vquad = vquad + step(i)*(GQ(i)+HALF*HD(i))
c     &        + XOPT(i)*HD(i)
      vquad= zero
      do 163 i=1,n
  163 vquad = vquad + step(i)*(gbeg(i)+HALF*HD(i))
      if (vquad.gt.zero) then
         print *," Warning: the TR subproblem was not well solved!"
         t = zero
         do i=1,n
           t = t + step(i)**2
         enddo
         print *, " vquad=", vquad, " Stepsize=",dsqrt(t)
         if (dsqrt(t).ge.half*DELTA) stop
      endif
      RETURN
C
C     The following instructions act as a subroutine for setting the vector
C     HD to the vector D multiplied by the second derivative matrix of Q.
C     They are called from three different places, which are distinguished
C     by the value of ITERC.
C
  170 continue

      DO 315 I=1,N
  315 HD(I) = ZERO
      IH=0
      DO 320 J=1,N
      DO 320 I=1,J
      IH=IH+1
      IF (I .LT. J) HD(J)=HD(J)+HQ(IH)*D(I)
  320 HD(I)=HD(I)+HQ(IH)*D(J) 

      IF (ITERC .EQ. 0) GOTO 20
      IF (ITERC .LE. ITERSW) GOTO 50
      GOTO 120
      END

      subroutine f_grad(mv,v_base,v_gtemp)  
      integer mv
      double precision v_base(*), v_gtemp(*)
      integer m1

      do 10 m1=1,mv
         v_gtemp(m1)=2.d0*v_base(m1)
   10 continue  

      return
      end

