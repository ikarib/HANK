SUBROUTINE IterateBellman

USE Parameters
USE Globals
USE Procedures
USE mUMFPACK

IMPLICIT NONE

INTEGER			:: iaby,ii
REAL(8)			:: lVdiff,ltau,ltau0,ladrift,lc,lh,llabdisutil

delta = deltass
	
!set drifts
ltau = 15.0
ltau0 = (ra+PerfectAnnuityMarkets*deathrate)*(agrid(ngpa)*0.999)**(1.0-ltau)
adrift = ra*agrid + PerfectAnnuityMarkets*deathrate*agrid - ltau0*(agrid**ltau)
bdrift = MERGE((rb+PerfectAnnuityMarkets*deathrate)*bgrid,(rborr+PerfectAnnuityMarkets*deathrate)*bgrid,bgrid>0.0)

IF (Borrowing==1 .and. bgrid(1) < -lumptransfer/(rborr+PerfectAnnuityMarkets*deathrate)) THEN
	write(*,*) 'Warning: natural borrowing limit violated'
END IF


!Initial Guess
IF	(	(EquilibriumR==0 .and. calibrating==.false.) &
 	.or. (EquilibriumR==1 .and. neqmiter<=3 .and. calibrating==.false.) &
 	.or. (CalibrateDiscountRate==1 .and. neqmiter<=3 .and. calibrating==.false.) &
 	.or. (CalibrateRhoAtInitialGuess==1 .and. neqmiter<=3 .and. calibrating==.false.) &
	.or. (calibrating==.true. .and. ImposeEqumInCalibration==1  .and. neqmiter==1 ) &
	.or. (calibrating==.true. .and. ImposeEqumInCalibration==0 ))  THEN
	
	!$OMP PARALLEL DO PRIVATE (ladrift,lh,lc,llabdisutil)
	DO iaby = 1,naby
		IF(LaborSupplyGHH==1) 	lh = (netwage*ygrid(yfromaby(iaby)))**frisch
		IF(LaborSupplySep==1) 	lh = 1.0/3.0
		IF(NoLaborSupply==1) 	lh = 1.0
  		ladrift = (ra + PerfectAnnuityMarkets*deathrate)*agrid(afromaby(iaby))
		lc = netwage*ygrid(yfromaby(iaby))*lh + lumptransfer + (rb+PerfectAnnuityMarkets*deathrate)*bgrid(bfromaby(iaby))
		llabdisutil = lh**(1.0+1.0/frisch)/(1.0+1.0/frisch)	
		IF(NoLaborSupply==1 .or. LaborSupplyGHH==1)V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = utilfn(lc)  / (rho+deathrate)
		IF(LaborSupplySep==1)V(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = (utilfn(lc) - chi*llabdisutil) / (rho+deathrate)
	END DO
	!$OMP END PARALLEL DO
	
END IF



ii = 1
lVdiff = 1.0
DO WHILE (ii<=maxiter .and. lVdiff>Vtol)

	CALL HJBUpdate

	!check for convergence
	lVdiff = maxval(abs(Vnew-V))
	IF(Display>=2) write(*,*) "Iteration: ",ii," max V change: ",lVdiff
	
	!update
	V = Vnew
	ii = ii+1
	
END DO



END SUBROUTINE IterateBellman