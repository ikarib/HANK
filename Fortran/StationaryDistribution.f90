SUBROUTINE StationaryDistribution

USE Parameters
USE Globals
USE Procedures
USE mUMFPACK

IMPLICIT NONE

INTEGER			:: ia,ib,iy,iaby,iab,it,iw(naby)
INTEGER, ALLOCATABLE	:: iwk(:)
REAL(8) 		:: ldiff
REAL(8), DIMENSION(nab) 	:: ldavec,ldbvec,ldiag,lgmat
REAL(8), DIMENSION(ngpa,ngpb) 	:: labdist
REAL(8), DIMENSION(ngpy,ngpy) 	:: leye,lmat
REAL(8), DIMENSION(nab,ngpy) 	:: lgmat1
TYPE(tCSR_di), DIMENSION(ngpy) 	 :: AUMF     !umfpack type

!identity matrix
leye = 0.0; DO iy = 1,ngpy
	leye(iy,iy) = 1.0
END DO


!$OMP PARALLEL DO
DO iab = 1,nab
	ldavec(iab) = adelta(afromab(iab))
	ldbvec(iab) = bdelta(bfromab(iab))
END DO
!$OMP END PARALLEL DO


IF	((EquilibriumR==0 .and. calibrating==.false.) &
  	.or. (EquilibriumR==1 .and. neqmiter<=100 .and. calibrating==.false.) &
	.or. (CalibrateDiscountRate==1 .and. neqmiter<=100 .and. calibrating==.false.) &
 	.or. (CalibrateRhoAtInitialGuess==1 .and. neqmiter<=100 .and. calibrating==.false.) &
	.or. (calibrating==.true. .and. ImposeEqumInCalibration==1  .and. neqmiter==1 ) &
	.or. (calibrating==.true. .and. ImposeEqumInCalibration==0 ) ) THEN

	gmat = 0.0

	!$OMP PARALLEL DO PRIVATE(iab,ia,ib)
	DO iy = 1,ngpy
		DO iab = 1,nab
			ia = afromab(iab)
			ib = bfromab(iab)
			IF(deathrate==0.0) THEN
				IF(Borrowing==0 .and. ia==1 .and. ib==2) gmat(iab,iy) = ydist(iy)
				IF(Borrowing==1 .and. ia==1 .and. ib==ngpbNEG+2) gmat(iab,iy) = ydist(iy)
!				IF(Borrowing==1 .and. ia==2 .and. ib==ngpbNEG+2) gmat(iab,:) = ydist
			ELSE IF(deathrate>0.0) THEN
				IF(Borrowing==0 .and. ia==1 .and. ib==1) gmat(iab,iy) = ydist(iy)
				IF(Borrowing==1 .and. ia==1 .and. ib==ngpbNEG+1) gmat(iab,iy) = ydist(iy)
			END IF

		END DO
 		gmat(:,iy) = ydist(iy)*gmat(:,iy)/SUM(gmat(:,iy)*ldavec*ldbvec)
	END DO
	!$OMP END PARALLEL DO
END IF

!$OMP PARALLEL DO PRIVATE(iwk,iw,ldiag)
DO iy = 1,ngpy
	IF(ALLOCATED(iwk)) DEALLOCATE(iwk)
	
	IF(ALLOCATED(BCSR(iy)%val)) DEALLOCATE(BCSR(iy)%val)
	IF(ALLOCATED(BCSR(iy)%col)) DEALLOCATE(BCSR(iy)%col)
	IF(ALLOCATED(BCSR(iy)%row)) DEALLOCATE(BCSR(iy)%row)
	
	ALLOCATE(iwk(AUCSR(iy)%nz))
	ALLOCATE(BCSR(iy)%val(AUCSR(iy)%nz), BCSR(iy)%col(AUCSR(iy)%nz), BCSR(iy)%row(AUCSR(iy)%n+1))
	
	!transpose AU matrix
	BCSR(iy) = AUCSR(iy)
	CALL csrcsc (nab, 1, 1, AUCSR(iy)%val, AUCSR(iy)%col, AUCSR(iy)%row, BCSR(iy)%val, BCSR(iy)%col, BCSR(iy)%row )

	!adjust A' matrix for non-linearly spaced grids
	ldiag = 1.0/(ldavec*ldbvec)
	CALL diamua (nab, 1, BCSR(iy)%val, BCSR(iy)%col, BCSR(iy)%row, ldiag,BCSR(iy)%val, BCSR(iy)%col, BCSR(iy)%row )
	ldiag = ldavec*ldbvec
	CALL amudia (nab, 1, BCSR(iy)%val, BCSR(iy)%col, BCSR(iy)%row, ldiag,BCSR(iy)%val, BCSR(iy)%col, BCSR(iy)%row )
	

 	BCSR(iy)%val = -deltakfe*BCSR(iy)%val

	ldiag = 1.0  - deltakfe*ymarkovdiag(iy,iy) + deltakfe*deathrate
	CALL apldia (nab, 0, BCSR(iy)%val, BCSR(iy)%col, BCSR(iy)%row, ldiag, BCSR(iy)%val, BCSR(iy)%col, BCSR(iy)%row, iw )
	
	AUMF(iy) = tCSR_di(BCSR(iy)%row-1,BCSR(iy)%col-1,BCSR(iy)%val) 
	
	
END DO
!$OMP END PARALLEL DO


lmat = leye + deltakfe*TRANSPOSE(ymarkovoff) 

ldiff = 1.0
it = 1

DO WHILE (ldiff>KFEtol .and. it<maxiterKFE)

	!sweep over y
	!$OMP PARALLEL DO PRIVATE(lgmat)
	DO iy = 1,ngpy	

 		lgmat = MATMUL(gmat(:,:),lmat(iy,:))
 		lgmat = lgmat + deltakfe*deathrate*MERGE(SUM(gmat(:,iy)*ldavec*ldbvec),0.0_8,afromab==1 .and. bfromab==ngpbNEG+1 )/(ldavec*ldbvec)
		lgmat1(:,iy) = AUMF(iy) .umfpack. lgmat		
		
	END DO
	!$OMP END PARALLEL DO
		
	ldiff = maxval(abs(gmat-lgmat1))
	gmat = lgmat1
	it = it+1
	IF (Display>=2) write(*,*) 'KFE iteration ',it,' max g change: ',ldiff
END DO

gmat = MERGE(0.0_8,gmat,abs(gmat)<1.0e-50_8)

!$OMP PARALLEL DO PRIVATE(ia,ib,iy,iab)	
DO iaby = 1,naby
	ia = afromaby(iaby)
	ib = bfromaby(iaby)
	iy = yfromaby(iaby)
	iab = 	abfromab(ia,ib)
	gjoint(ia,ib,iy) = gmat(iab,iy)
	gvec(iaby) = gmat(iab,iy)
END DO
!$OMP END PARALLEL DO

!marginal distributions
DO ia = 1,ngpa
DO iy = 1,ngpy
	gamarg(ia,iy) = SUM(gjoint(ia,:,iy)*bdelta)
END DO
END DO

DO ib = 1,ngpb
DO iy = 1,ngpy
	gbmarg(ib,iy) = SUM(gjoint(:,ib,iy)*adelta)
END DO
END DO

!ab joint distribution
DO ib = 1,ngpb
	gabmarg(:,ib) = SUM(gjoint(:,ib,:),DIM=2)*adelta*bdelta(ib)
END DO

!ab cumulative distribution
labdist(1,:) = gabmarg(1,:)
DO ia =2,ngpa
	labdist(ia,:) = labdist(ia-1,:)+ gabmarg(ia,:)
END DO
gabcum(:,1) = labdist(:,1)
DO ib = 2,ngpb
	gabcum(:,ib) = gabcum(:,ib-1)+ labdist(:,ib)
END DO


END SUBROUTINE StationaryDistribution