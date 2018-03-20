SUBROUTINE DiscountedMPC

USE Parameters
USE Globals
USE Procedures
USE mUMFPACK

IMPLICIT NONE

INTEGER			:: iaby,it,iw(naby),ia,ib,iy,maxiter_discMPC
REAL(8) 		:: leta,ldiff,discMPCtol
REAL(8), DIMENSION(ngpa,ngpb,ngpy) 	:: lmpd_a,lmpd_b
REAL(8), DIMENSION(nab) 	:: ldiag,lvec
REAL(8), DIMENSION(nab,ngpy) 	:: lsubeff1ass,lsubeff1ass1,lsubeff2ass,lsubeff2ass1,lwealtheff1ass,lwealtheff1ass1,lwealtheff2ass,lwealtheff2ass1
REAL(8), DIMENSION(nab,ngpy) 	:: lmpcvec,lmpdavec,lmpdbvec,lmargadjcostvec,lbadjvec,leffdiscvec,lmucvec
TYPE(tCSR_di), DIMENSION(ngpy) 	:: AUMF     !umfpack type
TYPE(CSR)						:: lA

IF (Display>=1) write(*,*) "Computing discounted MPC and wealth in steady state"
maxiter_discMPC = 50
discMPCtol = 1.0e-10

leta = -log(MonetaryShockPers)


!compute MPCs
DO ia = 1,ngpa
	DO iy = 1,ngpy
		DO ib = 1,ngpb-1
			mpc(ia,ib,iy) = (c(ia,ib+1,iy)-c(ia,ib,iy)) / dbgrid(ib)
		END DO
		mpc(ia,ngpb,iy) = mpc(ia,ngpb-1,iy)
	END DO	
END DO		

!compute MPDs
DO ia = 1,ngpa
	DO iy = 1,ngpy
		DO ib = 1,ngpb-1
			lmpd_b(ia,ib,iy) = (d(ia,ib+1,iy)-d(ia,ib,iy)) / dbgrid(ib)
		END DO
		lmpd_b(ia,ngpb,iy) = lmpd_b(ia,ngpb-1,iy)
	END DO	
END DO		

DO ib = 1,ngpb
	DO iy = 1,ngpy
		DO ia = 1,ngpa-1
			lmpd_a(ia,ib,iy) = (d(ia+1,ib,iy)-d(ia,ib,iy)) / dagrid(ia)
		END DO
		lmpd_a(ngpa,ib,iy) = lmpd_a(ngpa-1,ib,iy)
	END DO	
END DO		


!vectors of useful objectys
!$OMP PARALLEL DO PRIVATE(ia,ib,iy)
DO iaby = 1,naby
	ia = afromaby(iaby)
	ib = bfromaby(iaby)
	iy = yfromaby(iaby)
	
	lmpcvec(abfromab(ia,ib),iy) = mpc(ia,ib,iy)	
	lmpdavec(abfromab(ia,ib),iy) = lmpd_a(ia,ib,iy)
	lmpdbvec(abfromab(ia,ib),iy) = lmpd_b(ia,ib,iy)
	lmargadjcostvec(abfromab(ia,ib),iy) = adjcostfn1(d(ia,ib,iy), agrid(ia))
	leffdiscvec(abfromab(ia,ib),iy) = lmpcvec(abfromab(ia,ib),iy) + (1.0 + lmargadjcostvec(abfromab(ia,ib),iy))*lmpdbvec(abfromab(ia,ib),iy) - lmpdavec(abfromab(ia,ib),iy)
	lmucvec(abfromab(ia,ib),iy) = utilfn1(c(ia,ib,iy))
	lbadjvec(abfromab(ia,ib),iy) = bgrid(ib) * (mpc(ia,ib,iy)/c(ia,ib,iy)) *lmucvec(abfromab(ia,ib),iy)
END DO
!$OMP END PARALLEL DO

!SUBSTITUTION EFFECT: ONE ASSET: effective discounting of xi (see notes) and decay leta
!have to iterate because need to solve for all y points simultaneously

!initialize iterations without y transitions
DO iy = 1,ngpy
	lA = AUCSR(iy)
	lA%val = -lA%val	
	ldiag = leta + rho - rb +lmpcvec(:,iy) 
	CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
	AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 
END DO

!$OMP PARALLEL DO PRIVATE(lvec)
DO iy = 1,ngpy
	lvec = lmucvec(:,iy)
	lsubeff1ass(:,iy) = AUMF(iy) .umfpack. lvec
END DO
!$OMP END PARALLEL DO 

!iterate to convergence
DO iy = 1,ngpy
	lA = AUCSR(iy)
	lA%val = -lA%val	
	ldiag = leta + rho - rb +lmpcvec(:,iy) - ymarkovdiag(iy,iy)	
	CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
	AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 
END DO

it = 1
ldiff = 100.0
lsubeff1ass1 = lsubeff1ass
DO WHILE (ldiff>discMPCtol .and. it<maxiter_discMPC)
	
	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		lvec = lmucvec(:,iy) + MATMUL(ymarkovoff(iy,:),TRANSPOSE(lsubeff1ass1(:,:)))
		lsubeff1ass(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 
	
	ldiff = maxval(abs(lsubeff1ass-lsubeff1ass1))
	lsubeff1ass1 = lsubeff1ass
	IF (Display>=2) write(*,*) 'Sub effect 1 asset iteration ',it,' max change: ',ldiff
	it = it+1

END DO

 
!SUBSTITUTION EFFECT: TWO ASSET: effective discounting of (MPC + (1+ chi'(d)*MPD_B-MPD_A) and decay leta
!have to iterate because need to solve for all y points simultaneously

!initialize iterations without y transitions
DO iy = 1,ngpy
	lA = AUCSR(iy)
	lA%val = -lA%val	
	ldiag = leta + rho - rb +leffdiscvec(:,iy)
	CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
	AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 
END DO

!$OMP PARALLEL DO PRIVATE(lvec)
DO iy = 1,ngpy
	lvec = lmucvec(:,iy)
	lsubeff2ass(:,iy) = AUMF(iy) .umfpack. lvec
END DO
!$OMP END PARALLEL DO 

!iterate to convergence
DO iy = 1,ngpy
	lA = AUCSR(iy)
	lA%val = -lA%val	
	ldiag = leta+ rho - rb  +leffdiscvec(:,iy) - ymarkovdiag(iy,iy)	
	CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
	AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 
END DO

it = 1
ldiff = 100.0
lsubeff2ass1 = lsubeff2ass
DO WHILE (ldiff>discMPCtol .and. it<maxiter_discMPC)
	
	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		lvec = lmucvec(:,iy) + MATMUL(ymarkovoff(iy,:),TRANSPOSE(lsubeff2ass1(:,:)))
		lsubeff2ass(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 
	
	ldiff = maxval(abs(lsubeff2ass-lsubeff2ass1))
	lsubeff2ass1 = lsubeff2ass
	IF (Display>=2) write(*,*) 'Sub effect 2 asset iteration ',it,' max change: ',ldiff
	it = it+1

END DO


!WEALTH EFFECT: ONE ASSET: effective discounting of xi and decay leta
!have to iterate because need to solve for all y points simultaneously

!initialize iterations without y transitions
DO iy = 1,ngpy
	lA = AUCSR(iy)
	lA%val = -lA%val	
	ldiag = leta + rho - rb +lmpcvec(:,iy)
	CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
	AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 
END DO

!$OMP PARALLEL DO PRIVATE(lvec)
DO iy = 1,ngpy
	lvec = lbadjvec(:,iy)
	lwealtheff1ass(:,iy) = AUMF(iy) .umfpack. lvec
END DO
!$OMP END PARALLEL DO 

!iterate to convergence
DO iy = 1,ngpy
	lA = AUCSR(iy)
	lA%val = -lA%val	
	ldiag = leta + rho - rb +lmpcvec(:,iy) - ymarkovdiag(iy,iy)	
	CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
	AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 
END DO

it = 1
ldiff = 100.0
lwealtheff1ass1 = lwealtheff1ass
DO WHILE (ldiff>discMPCtol .and. it<maxiter_discMPC)
	
	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		lvec = lbadjvec(:,iy) + MATMUL(ymarkovoff(iy,:),TRANSPOSE(lwealtheff1ass1(:,:)))
		lwealtheff1ass(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 
	
	ldiff = maxval(abs(lwealtheff1ass-lwealtheff1ass1))
	lwealtheff1ass1 = lwealtheff1ass
	IF (Display>=2) write(*,*) 'Sub effect 1 asset iteration ',it,' max change: ',ldiff
	it = it+1

END DO

 
!WEALTH EFFECT: TWO ASSET: effective discounting of (MPC + (1+ chi'(d)*MPD_B-MPD_A) and decay leta
!have to iterate because need to solve for all y points simultaneously

!initialize iterations without y transitions
DO iy = 1,ngpy
	lA = AUCSR(iy)
	lA%val = -lA%val	
	ldiag = leta+ rho - rb +leffdiscvec(:,iy)
	CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
	AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 
END DO

!$OMP PARALLEL DO PRIVATE(lvec)
DO iy = 1,ngpy
	lvec = lbadjvec(:,iy)
	lwealtheff2ass(:,iy) = AUMF(iy) .umfpack. lvec
END DO
!$OMP END PARALLEL DO 

!iterate to convergence
DO iy = 1,ngpy
	lA = AUCSR(iy)
	lA%val = -lA%val	
	ldiag = leta + rho - rb +leffdiscvec(:,iy) - ymarkovdiag(iy,iy)	
	CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
	AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 
END DO

it = 1
ldiff = 100.0
lwealtheff2ass1 = lwealtheff2ass
DO WHILE (ldiff>discMPCtol .and. it<maxiter_discMPC)
	
	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		lvec = lbadjvec(:,iy) + MATMUL(ymarkovoff(iy,:),TRANSPOSE(lwealtheff2ass1(:,:)))
		lwealtheff2ass(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 
	
	ldiff = maxval(abs(lwealtheff2ass-lwealtheff2ass1))
	lwealtheff2ass1 = lwealtheff2ass
	IF (Display>=2) write(*,*) 'Sub effect 2 asset iteration ',it,' max change: ',ldiff
	it = it+1

END DO


!save functions
!$OMP PARALLEL DO
DO iaby = 1,naby
	subeff1ass(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lsubeff1ass(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
	subeff2ass(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lsubeff2ass(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
	wealtheff1ass(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lwealtheff1ass(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
	wealtheff2ass(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lwealtheff2ass(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
END DO
!$OMP END PARALLEL DO




END SUBROUTINE DiscountedMPC