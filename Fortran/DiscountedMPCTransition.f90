SUBROUTINE DiscountedMPCTransition

USE Parameters
USE Globals
USE Procedures
USE mUMFPACK

IMPLICIT NONE

INTEGER			:: iaby,it,iw(naby),ia,ib,iy
REAL(8) 		:: ldrb
REAL(8), DIMENSION(nab) 	:: ldiag,lvec
REAL(8), DIMENSION(ngpa,ngpb,ngpy) 	:: lmpd_a,lmpd_b
REAL(8), DIMENSION(nab,ngpy) 	:: lsubeff1ass,lsubeff1ass1,lsubeff2ass,lsubeff2ass1,lwealtheff1ass,lwealtheff1ass1,lwealtheff2ass,lwealtheff2ass1
REAL(8), DIMENSION(nab,ngpy) 	:: lmpcvec,lmpdavec,lmpdbvec,lmargadjcostvec,lbadjvec,leffdiscvec,lmucvec


TYPE(tCSR_di), DIMENSION(ngpy) 	:: AUMF     !umfpack type
TYPE(CSR)						:: lA

IF(Display>=1) WRITE(*,*) " Solving discounted MPCs backward "

DO it = Ttransition,1,-1
	IF(Display>=2) WRITE(*,*) "   Solving discounted MPCs backward at time: ",it

	delta = deltatransvec(it)
	
	rho = irfpointer%equmSTICKY(it)%rho
	rb = irfpointer%equmSTICKY(it)%rb
	ldrb = irfpointer%equmSTICKY(it)%rb - equmINITSS%rb
	c = irfpointer%solnSTICKY(it)%c
	d = irfpointer%solnSTICKY(it)%d
	AUCSR = irfpointer%solnSTICKY(it)%AU
	
	IF(it==Ttransition) THEN
		subeff1ass = solnFINALSS%subeff1ass
		subeff2ass = solnFINALSS%subeff2ass
		wealtheff1ass = solnFINALSS%wealtheff1ass
		wealtheff2ass = solnFINALSS%wealtheff2ass
	END IF	
	
	IF(it<Ttransition) THEN
		subeff1ass = irfpointer%solnSTICKY(it+1)%subeff1ass
		subeff2ass = irfpointer%solnSTICKY(it+1)%subeff2ass
		wealtheff1ass = irfpointer%solnSTICKY(it+1)%wealtheff1ass
		wealtheff2ass = irfpointer%solnSTICKY(it+1)%wealtheff2ass
	END IF

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
		lmucvec(abfromab(ia,ib),iy) = utilfn1(c(ia,ib,iy));
		lbadjvec(abfromab(ia,ib),iy) = bgrid(ib) * (mpc(ia,ib,iy)/c(ia,ib,iy)) * lmucvec(abfromab(ia,ib),iy)
		
		lsubeff1ass1(abfromab(ia,ib),iy) = subeff1ass(ia,ib,iy)	
		lsubeff2ass1(abfromab(ia,ib),iy) = subeff2ass(ia,ib,iy)	
		lwealtheff1ass1(abfromab(ia,ib),iy) = wealtheff1ass(ia,ib,iy)	
		lwealtheff2ass1(abfromab(ia,ib),iy) = wealtheff2ass(ia,ib,iy)	
	END DO
	!$OMP END PARALLEL DO
	
	
	
	!SUBSTITUTION EFFECT: ONE ASSET
	DO iy = 1,ngpy
		lA = AUCSR(iy)
		lA%val = -delta*lA%val	
		ldiag = 1.0 + delta*(rho-rb +lmpcvec(:,iy)) - delta*ymarkovdiag(iy,iy)
		CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
		AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 		
	END DO

	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		lvec = delta*lmucvec(:,iy)*ldrb + lsubeff1ass1(:,iy) + delta*MATMUL(ymarkovoff(iy,:),TRANSPOSE(lsubeff1ass1(:,:)))
		lsubeff1ass(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 
	

	!SUBSTITUTION EFFECT: TWO ASSET
	DO iy = 1,ngpy
		lA = AUCSR(iy)
		lA%val = -delta*lA%val	
		ldiag = 1.0 + delta*(rho-rb+leffdiscvec(:,iy)) - delta*ymarkovdiag(iy,iy)
		CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
		AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 		
	END DO

	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		lvec = delta*lmucvec(:,iy)*ldrb + lsubeff2ass1(:,iy) + delta*MATMUL(ymarkovoff(iy,:),TRANSPOSE(lsubeff2ass1(:,:)))
		lsubeff2ass(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 
	
	
	!WEALTH EFFECT: ONE ASSET
	DO iy = 1,ngpy
		lA = AUCSR(iy)
		lA%val = -delta*lA%val	
		ldiag = 1.0 + delta*(rho-rb+lmpcvec(:,iy)) - delta*ymarkovdiag(iy,iy)
		CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
		AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 		
	END DO

	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		lvec = delta*lbadjvec(:,iy)*ldrb + lwealtheff1ass1(:,iy) + delta*MATMUL(ymarkovoff(iy,:),TRANSPOSE(lwealtheff1ass1(:,:)))
		lwealtheff1ass(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 
	

	!WEALTH EFFECT: TWO ASSET: effective discounting of MPC and decay leta
	DO iy = 1,ngpy
		lA = AUCSR(iy)
		lA%val = -delta*lA%val	
		ldiag = 1.0 + delta*(rho-rb +leffdiscvec(:,iy)) - delta*ymarkovdiag(iy,iy)
		CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
		AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 		
	END DO

	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		lvec = delta*lbadjvec(:,iy)*ldrb + lwealtheff2ass1(:,iy) + delta*MATMUL(ymarkovoff(iy,:),TRANSPOSE(lwealtheff2ass1(:,:)))
		lwealtheff2ass(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 
	
	
	
	!save functions
	!$OMP PARALLEL DO
	DO iaby = 1,naby
		subeff1ass(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lsubeff1ass(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
		subeff2ass(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lsubeff2ass(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
		wealtheff1ass(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lwealtheff1ass(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
		wealtheff2ass(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lwealtheff2ass(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
	END DO
	!$OMP END PARALLEL DO
	
	
	irfpointer%solnSTICKY(it)%mpc = mpc
	irfpointer%solnSTICKY(it)%subeff1ass  = subeff1ass
	irfpointer%solnSTICKY(it)%subeff2ass  = subeff2ass
	irfpointer%solnSTICKY(it)%wealtheff1ass  = wealtheff1ass
	irfpointer%solnSTICKY(it)%wealtheff2ass  = wealtheff2ass
END DO

END SUBROUTINE DiscountedMPCTransition