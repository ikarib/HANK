SUBROUTINE HJBUpdate

USE Parameters
USE Globals
USE Procedures
USE mUMFPACK

IMPLICIT NONE

INTEGER			:: ia,ib,iy,iaby,iab,ii,iu,iw(naby)
REAL(8) 		:: lVaF,lVaB,lVbF,lVbB,lcF,lcB,lc0,lsF,lsB,ls0,lHcF,lHcB,lHc0,ldFB,ldBF,ldBB,lHdFB,lHdBF,lHdBB,lhF,lhB,lh0,llabdisutil
INTEGER 		:: validF,validB,validFB,validBF,validBB
REAL(8) 		:: ladriftB,ladriftF,lbdriftB,lbdriftF,lval,laudriftB,laudriftF,lbudriftB,lbudriftF
REAL(8), DIMENSION(nab) :: 	luvec,lbvec,lVvec,ldiag
TYPE(tCSR_di), DIMENSION(ngpy) 	 :: AUMF     !umfpack type

 c = 0.0
h = 0.0	
d = 0.0
u = 0.0
s = 0.0

!$OMP PARALLEL DO PRIVATE(ia,ib,iy,ii,iw,lVaF,lVaB,lVbF,lVbB,lcF,lcB,lc0,lsF,lsB,ls0,lHcF,lHcB,lHc0,ldFB,ldBF,ldBB,lHdFB,lHdBF,lHdBB,lhF,lhB,lh0,validF,validB,validFB,validBF,validBB,llabdisutil)
DO iaby = 1,naby
	ia = afromaby(iaby)
	ib = bfromaby(iaby)
	iy = yfromaby(iaby)

	!derivatives wrt a: forward
	IF (ia<ngpa) lVaF = (V(ia+1,ib,iy) - V(ia,ib,iy))/dagrid(ia)
	
	!derivatives wrt a: backward
	IF (ia>1) lVaB = (V(ia,ib,iy) - V(ia-1,ib,iy))/dagrid(ia-1)
	
	!derivatives wrt b: forward
	IF (ib<ngpb) lVbF = (V(ia,ib+1,iy) - V(ia,ib,iy))/dbgrid(ib)	

	!derivatives wrt b: backward
	IF (ib>1) lVbB = (V(ia,ib,iy) - V(ia,ib-1,iy))/dbgrid(ib-1)
	
	IF (ia<ngpa) lVaF = max(lVaF,dVamin)
	IF (ia>1) lVaB = max(lVaB,dVamin)
	IF (ib<ngpb) lVbF = max(lVbF,dVbmin)
	IF (ib>1) lVbB = max(lVbB,dVbmin)
	
	
	gidioprod = ygrid(iy)
	gnetwage = netwage*gidioprod
	gbdrift = bdrift(ib) + lumptransfer
	IF(DistributeProfitsInProportion==1 .and. TaxHHProfitIncome==0) gbdrift = gbdrift + (1.0-profdistfrac)*(1.0-corptax)*profit*ygrid(iy)/meanlabeff
	IF(DistributeProfitsInProportion==1 .and. TaxHHProfitIncome==1) gbdrift = gbdrift + (1.0-labtax)*(1.0-profdistfrac)*(1.0-corptax)*profit*ygrid(iy)/meanlabeff
	gill = agrid(ia)
	
	!consumption decision
	IF(ib<ngpb) THEN
		CALL OptimalConsumption(lVbF,lcF,lhF,lsF,lHcF)
	ELSE
		lsF = 0.0
		lHcF = -1.0e12
	END IF
	validF = 0
	IF(lsF>0.0) validF = 1
	
	IF (ib>1) THEN		
		CALL OptimalConsumption(lVbB,lcB,lhB,lsB,lHcB)
	ELSE
		CALL OptimalConsumption( -999.9_8 ,lcB,lhB,lsB,lHcB)
	END IF				 
	validB = 0
	IF(lsB<0.0) validB = 1
	
	CALL OptimalConsumption(-999.9_8,lc0,lh0,ls0,lHc0)
		
	IF (validF==1 .and. (validB==0 .or. lHcF>=lHcB) .and. lHcF>=lHc0) THEN	!forward
		c(ia,ib,iy) = lcF
		h(ia,ib,iy) = lhF
		s(ia,ib,iy) = lsF
	ELSE IF (validB==1 .and. (validF==0 .or. lHcB>=lHcF) .and. lHcB>=lHc0) THEN	!backward
		c(ia,ib,iy) = lcB
		h(ia,ib,iy) = lhB
		s(ia,ib,iy) = lsB	
	ELSE 
		c(ia,ib,iy) = lc0
		h(ia,ib,iy) = lh0
		s(ia,ib,iy) = ls0
	END IF
	
	!deposit decision
	IF(ia<ngpa .and. ib>1) THEN !a forward, b backward
		ldFB = adjcostfn1inv(lVaF/lVbB-1.0,agrid(ia))
		lHdFB = lVaF*ldFB - lVbB*(ldFB + adjcostfn(ldFB,agrid(ia)))
		validFB = 0
		IF(ldFB>0.0 .and. lHdFB>0.0) validFB = 1
	ELSE
		validFB = 0
		lHdFB = -1.0e12
	END IF

	IF (ia>1 .and. ib<ngpb) THEN !a backward, b forward
		ldBF = adjcostfn1inv(lVaB/lVbF-1.0,agrid(ia))
		lHdBF = lVaB*ldBF - lVbF*(ldBF + adjcostfn(ldBF,agrid(ia)))
		validBF = 0
		IF(ldBF <= -adjcostfn(ldBF,agrid(ia)) .and. lHdBF>0.0) validBF = 1
	ELSE
		validBF = 0
		lHdBF = -1.0e12
	END IF

	IF(ia>1 ) THEN !a backward, b backward
 		IF (ib==1) THEN
			IF(ScaleDisutilityIdio==0) llabdisutil = lhB**(1.0+1.0/frisch)/(1.0+1.0/frisch)
			IF(ScaleDisutilityIdio==1) llabdisutil = ygrid(iy)*lhB**(1.0+1.0/frisch)/(1.0+1.0/frisch)
			IF(NoLaborSupply==1 .or. LaborSupplySep==1) lVbB = utilfn1(lcB)
			IF(LaborSupplyGHH==1) lVbB = utilfn1(lcB-chi*llabdisutil)
		END IF
		ldBB = adjcostfn1inv(lVaB/lVbB-1.0,agrid(ia))
		lHdBB = lVaB*ldBB - lVbB*(ldBB + adjcostfn(ldBB,agrid(ia)))
		validBB = 0
		IF(ldBB > -adjcostfn(ldBB,agrid(ia)) .and. ldBB<=0.0 .and. lHdBB>0.0) validBB = 1
	ELSE 
		validBB = 0
		lHdBB = -1.0e12		
	END IF

	IF (validFB==1 .and. (validBF==0 .or. lHdFB>=lHdBF) .and. (validBB==0 .or. lHdFB>=lHdBB)) THEN !forward, backward
		d(ia,ib,iy) = ldFB					
	ELSE IF ((validFB==0 .or. lHdBF>=lHdFB) .and. validBF==1 .and. (validBB==0 .or. lHdBF>=lHdBB)) THEN !backward,forward
		d(ia,ib,iy) = ldBF
	ELSE IF ((validFB==0 .or. lHdBB>=lHdFB) .and. (validBF==0 .or. lHdBB>=lHdBF) .and. validBB==1) THEN !backward,backward
		d(ia,ib,iy) = ldBB	
	ELSE IF	(validFB==0 .and. validBF==0 .and. validBB==0) THEN !none
		d(ia,ib,iy) = 0.0
	ELSE !more than 1
		write(*,*) "should never be here!"
	END IF
	


	IF(NoLaborSupply==1) u(ia,ib,iy) = utilfn(c(ia,ib,iy))
	IF (ScaleDisutilityIdio==0) llabdisutil = h(ia,ib,iy)**(1.0+1.0/frisch)/(1.0+1.0/frisch)
	IF (ScaleDisutilityIdio==1) llabdisutil = ygrid(iy)*h(ia,ib,iy)**(1.0+1.0/frisch)/(1.0+1.0/frisch)
	IF(LaborSupplySep==1) u(ia,ib,iy) = utilfn(c(ia,ib,iy)) - chi*llabdisutil
	IF(LaborSupplyGHH==1) u(ia,ib,iy) = utilfn(c(ia,ib,iy) - chi*llabdisutil)
	bdot(ia,ib,iy) = s(ia,ib,iy)-d(ia,ib,iy)- adjcostfn(d(ia,ib,iy),agrid(ia))
	
END DO
!$OMP END PARALLEL DO


! construct ACOO(iy) matrix in sparse COO form and vectors by filling in coefficients on each V(ia,ib,iy)

!$OMP PARALLEL DO PRIVATE(ii,iu,iab,ia,ib,iw,luvec,lbvec,lVvec,ldiag,ladriftB,ladriftF,lbdriftB,lbdriftF,lval,laudriftB,laudriftF,lbudriftB,lbudriftF)
DO iy = 1,ngpy
	ACOO(iy)%row = 0
	ACOO(iy)%col = 0
	ACOO(iy)%val = 0.0_8
	AUCOO(iy)%row = 0
	AUCOO(iy)%col = 0
	AUCOO(iy)%val = 0.0_8

	ii = 0
	iu = 0
	DO iab = 1,nab !cannot be done in parallel
		ia = afromab(iab)
		ib = bfromab(iab)
	
		!vector of constants
		luvec(iab) = u(ia,ib,iy)
		lbvec(iab) = delta*luvec(iab) + V(ia,ib,iy) + delta*DOT_PRODUCT(ymarkovoff(iy,:),V(ia,ib,:))

		!a drifts
		ladriftB = min(d(ia,ib,iy),0.0_8) + min(adrift(ia),0.0_8)
		ladriftF = max(d(ia,ib,iy),0.0_8) + max(adrift(ia),0.0_8)
	
		!b drift
		lbdriftB = min(-d(ia,ib,iy) - adjcostfn(d(ia,ib,iy),agrid(ia)),0.0) + min(s(ia,ib,iy),0.0_8)
		lbdriftF = max(-d(ia,ib,iy) - adjcostfn(d(ia,ib,iy),agrid(ia)),0.0) + max(s(ia,ib,iy),0.0_8)

		!a drift, upwind
		IF(ib<ngpb) laudriftB = min(d(ia,ib,iy) + adrift(ia),0.0_8)
		IF(ib<ngpb) laudriftF = max(d(ia,ib,iy) + adrift(ia),0.0_8)
		IF(ib==ngpb) laudriftB = min(d(ia,ib-1,iy) + adrift(ia),0.0_8)
		IF(ib==ngpb) laudriftF = max(d(ia,ib-1,iy) + adrift(ia),0.0_8)

		!b drift,upwind
		IF(ib<ngpb) lbudriftB = min(s(ia,ib,iy) -d(ia,ib,iy) - adjcostfn(d(ia,ib,iy),agrid(ia)),0.0_8)
		IF(ib<ngpb) lbudriftF = max(s(ia,ib,iy) -d(ia,ib,iy) - adjcostfn(d(ia,ib,iy),agrid(ia)),0.0_8)
		IF(ib==ngpb) lbudriftB = min(s(ia,ib,iy) -d(ia,ib-1,iy) - adjcostfn(d(ia,ib-1,iy),agrid(ia)),0.0_8)
		IF(ib==ngpb) lbudriftF = max(s(ia,ib,iy) -d(ia,ib-1,iy) - adjcostfn(d(ia,ib-1,iy),agrid(ia)),0.0_8)
		
		!a-1
		lval = 0.0_8
		IF(ia>1) lval = lval - ladriftB/dagrid(ia-1)
		IF(lval .ne. 0.0_8) THEN
			ii = ii+1
			ACOO(iy)%row(ii) = iab
			ACOO(iy)%col(ii) = abfromab(ia-1,ib)
			ACOO(iy)%val(ii) = lval
		END IF

		!a-1, upwind
		lval = 0.0_8
		IF(ia>1) lval = lval - laudriftB/dagrid(ia-1)
		IF(lval .ne. 0.0_8) THEN
			iu = iu+1
			AUCOO(iy)%row(iu) = iab
			AUCOO(iy)%col(iu) = abfromab(ia-1,ib)
			AUCOO(iy)%val(iu) = lval
		END IF

		!b-1
		lval = 0.0_8
		IF(ib>1) lval = lval - lbdriftB/dbgrid(ib-1)
		IF(lval .ne. 0.0_8) THEN
			ii = ii+1
			ACOO(iy)%row(ii) = iab
			ACOO(iy)%col(ii) = abfromab(ia,ib-1)
			ACOO(iy)%val(ii) = lval
		END IF

		!b-1, upwind
		lval = 0.0_8
		IF(ib>1) lval = lval - lbudriftB/dbgrid(ib-1)
		IF(lval .ne. 0.0_8) THEN
			iu = iu+1
			AUCOO(iy)%row(iu) = iab
			AUCOO(iy)%col(iu) = abfromab(ia,ib-1)
			AUCOO(iy)%val(iu) = lval
		END IF
	
		!diagonal value
		lval = 0.0_8
		IF(ia>1 .and. ia<ngpa) lval = lval + ladriftB/dagrid(ia-1) - ladriftF/dagrid(ia)
		IF(ia==1) lval = lval - ladriftF/dagrid(ia)
		IF(ia==ngpa) lval = lval + ladriftB/dagrid(ia-1)
		IF(ib>1 .and. ib<ngpb) lval = lval + lbdriftB/dbgrid(ib-1) - lbdriftF/dbgrid(ib)
		IF(ib==1) lval = lval - lbdriftF/dbgrid(ib)
		IF(ib==ngpb) lval = lval + lbdriftB/dbgrid(ib-1)
		
		
		! do not impose that diagonal elements are non-zero, keep them even if zero since B matrix below will modify them
		ii = ii+1
		ACOO(iy)%row(ii) = iab
		ACOO(iy)%col(ii) = iab
		ACOO(iy)%val(ii) = lval	

		!diagonal value, upwind
		lval = 0.0_8
		IF(ia>1 .and. ia<ngpa) lval = lval + laudriftB/dagrid(ia-1) - laudriftF/dagrid(ia)
		IF(ia==1) lval = lval - laudriftF/dagrid(ia)
		IF(ia==ngpa) lval = lval + laudriftB/dagrid(ia-1)
		IF(ib>1 .and. ib<ngpb) lval = lval + lbudriftB/dbgrid(ib-1) - lbudriftF/dbgrid(ib)
		IF(ib==1) lval = lval - lbudriftF/dbgrid(ib)
		IF(ib==ngpb) lval = lval + lbudriftB/dbgrid(ib-1)
		iu = iu+1
		AUCOO(iy)%row(iu) = iab
		AUCOO(iy)%col(iu) = iab
		AUCOO(iy)%val(iu) = lval

		!b+1
		lval = 0.0_8
		IF(ib<ngpb) lval = lval + lbdriftF/dbgrid(ib)
		IF(lval .ne. 0.0_8) THEN
			ii = ii+1
			ACOO(iy)%row(ii) = iab
			ACOO(iy)%col(ii) = abfromab(ia,ib+1)
			ACOO(iy)%val(ii) = lval
		END IF	

		!b+1, upwind
		lval = 0.0_8
		IF(ib<ngpb) lval = lval + lbudriftF/dbgrid(ib)
		IF(lval .ne. 0.0_8) THEN
			iu = iu+1
			AUCOO(iy)%row(iu) = iab
			AUCOO(iy)%col(iu) = abfromab(ia,ib+1)
			AUCOO(iy)%val(iu) = lval
		END IF	

		!a+1
		lval = 0.0_8
		IF(ia<ngpa) lval = lval + ladriftF/dagrid(ia)
		IF(lval .ne. 0.0_8) THEN
			ii = ii+1
			ACOO(iy)%row(ii) = iab
			ACOO(iy)%col(ii) = abfromab(ia+1,ib)
			ACOO(iy)%val(ii) = lval
		END IF

		!a+1, upwind
		lval = 0.0_8
		IF(ia<ngpa) lval = lval + laudriftF/dagrid(ia)
		IF(lval .ne. 0.0_8) THEN
			iu = iu+1
			AUCOO(iy)%row(iu) = iab
			AUCOO(iy)%col(iu) = abfromab(ia+1,ib)
			AUCOO(iy)%val(iu) = lval
		END IF
	
	END DO
	ACOO(iy)%nz = ii
	AUCOO(iy)%nz = iu


	!allocate CSR matrix and convert COO to CSR
	IF(ALLOCATED(ACSR(iy)%val)) DEALLOCATE(ACSR(iy)%val)
	IF(ALLOCATED(ACSR(iy)%col)) DEALLOCATE(ACSR(iy)%col)
	IF(ALLOCATED(ACSR(iy)%row)) DEALLOCATE(ACSR(iy)%row)
	ACSR(iy)%n = nab
	ACSR(iy)%nz = ACOO(iy)%nz
	ALLOCATE(ACSR(iy)%val(ACSR(iy)%nz), ACSR(iy)%col(ACSR(iy)%nz), ACSR(iy)%row(ACSR(iy)%n+1))
	CALL coocsr (nab, ACOO(iy)%nz, ACOO(iy)%val(1:ACOO(iy)%nz), ACOO(iy)%row(1:ACOO(iy)%nz), ACOO(iy)%col(1:ACOO(iy)%nz), ACSR(iy)%val, ACSR(iy)%col, ACSR(iy)%row)

	IF(ALLOCATED(AUCSR(iy)%val)) DEALLOCATE(AUCSR(iy)%val)
	IF(ALLOCATED(AUCSR(iy)%col)) DEALLOCATE(AUCSR(iy)%col)
	IF(ALLOCATED(AUCSR(iy)%row)) DEALLOCATE(AUCSR(iy)%row)
	AUCSR(iy)%n = nab
	AUCSR(iy)%nz = AUCOO(iy)%nz
	ALLOCATE(AUCSR(iy)%val(AUCSR(iy)%nz), AUCSR(iy)%col(AUCSR(iy)%nz), AUCSR(iy)%row(AUCSR(iy)%n+1))
	CALL coocsr (nab, AUCOO(iy)%nz, AUCOO(iy)%val(1:AUCOO(iy)%nz), AUCOO(iy)%row(1:AUCOO(iy)%nz), AUCOO(iy)%col(1:AUCOO(iy)%nz), AUCSR(iy)%val, AUCSR(iy)%col, AUCSR(iy)%row)

	!construct B matrix = I +delta*(rho*I -A): assumes that all diagonal terms in A matrix are non-zero
	IF(ALLOCATED(BCSR(iy)%val)) DEALLOCATE(BCSR(iy)%val)
	IF(ALLOCATED(BCSR(iy)%col)) DEALLOCATE(BCSR(iy)%col)
	IF(ALLOCATED(BCSR(iy)%row)) DEALLOCATE(BCSR(iy)%row)
	BCSR(iy)%nz = ACSR(iy)%nz
	BCSR(iy)%n = nab
	ALLOCATE(BCSR(iy)%val(BCSR(iy)%nz), BCSR(iy)%col(BCSR(iy)%nz), BCSR(iy)%row(BCSR(iy)%n+1))
	BCSR(iy)%val = -delta*ACSR(iy)%val
	BCSR(iy)%row = ACSR(iy)%row
	BCSR(iy)%col = ACSR(iy)%col
	ldiag = 1.0 + delta*(rho+deathrate) - delta*ymarkovdiag(iy,iy)
	CALL apldia (BCSR(iy)%n, 0, BCSR(iy)%val, BCSR(iy)%col, BCSR(iy)%row, ldiag, BCSR(iy)%val, BCSR(iy)%col, BCSR(iy)%row, iw )
	AUMF(iy) = tCSR_di(BCSR(iy)%row-1,BCSR(iy)%col-1,BCSR(iy)%val) 
	lVvec = AUMF(iy) .umfpack. lbvec

	DO iab = 1,nab
		Vnew(afromab(iab),bfromab(iab),iy) = lVvec(iab)
	END DO

END DO
!$OMP END PARALLEL DO


!check for non-monotonicity: can't do in parallel
IF (ReportNonMonotonicity==1) THEN
	DO iaby = 1,naby
		IF(bfromaby(iaby)>1) THEN
			IF (Vnew(afromaby(iaby),bfromaby(iaby),yfromaby(iaby))<Vnew(afromaby(iaby),bfromaby(iaby)-1,yfromaby(iaby)) -1.0e-8) &
				write(*,*) 'non-monotonicity, b dimension, ib = ',bfromaby(iaby)-1,' ia = ',afromaby(iaby)
		END IF
		IF(afromaby(iaby)>1) THEN
			IF (Vnew(afromaby(iaby),bfromaby(iaby),yfromaby(iaby))<Vnew(afromaby(iaby)-1,bfromaby(iaby),yfromaby(iaby)) -1.0e-8) &
				write(*,*) 'non-monotonicity, a dimension, ia = ',afromaby(iaby)-1,' ib = ',bfromaby(iaby)
		END IF
	END DO
END IF
	

END SUBROUTINE HJBUpdate