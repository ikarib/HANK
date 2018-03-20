SUBROUTINE CumulativeConsumption

USE Parameters
USE Globals
USE Procedures
USE mUMFPACK

IMPLICIT NONE

INTEGER			:: iy,iaby,it,iw(naby),nsteps
REAL(8) 		:: ldelta
REAL(8), DIMENSION(nab) 	:: ldiag,lvec
REAL(8), DIMENSION(nab,ngpy) 	:: lcvec,lccumvec,lccumvec1,ldvec,ldcumvec,ldcumvec1
TYPE(tCSR_di), DIMENSION(ngpy) 	:: AUMF     !umfpack type
TYPE(CSR)						:: lA

IF (Display>=1) write(*,*) "Computing cumulative consumption for MPC"

!$OMP PARALLEL DO
DO iaby = 1,naby
	lcvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby)) = c(afromaby(iaby),bfromaby(iaby),yfromaby(iaby))
	ldvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby)) = d(afromaby(iaby),bfromaby(iaby),yfromaby(iaby))
END DO
!$OMP END PARALLEL DO

ldelta = deltacumcon

DO iy = 1,ngpy
	lA = AUCSR(iy)
	lA%val = -ldelta*lA%val
	ldiag = 1.0 - ldelta*ymarkovdiag(iy,iy)	
	CALL apldia (nab, 0, lA%val, lA%col, lA%row, ldiag, lA%val, lA%col, lA%row, iw )
	AUMF(iy) = tCSR_di(lA%row-1,lA%col-1,lA%val) 
END DO

!one quarter cumulative expected consumption
nsteps = nint(1.0/ldelta)
lccumvec = 0.0
ldcumvec = 0.0
DO it = 1,nsteps

	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		!consumption	
		!consumption	
		lvec = lccumvec(:,iy) + ldelta*lcvec(:,iy) + ldelta*MATMUL(ymarkovoff(iy,:),TRANSPOSE(lccumvec(:,:)))
		lccumvec1(:,iy) = AUMF(iy) .umfpack. lvec
		!deposits
		lvec = ldcumvec(:,iy) + ldelta*ldvec(:,iy) + ldelta*MATMUL(ymarkovoff(iy,:),TRANSPOSE(ldcumvec(:,:)))
		ldcumvec1(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 

	lccumvec = lccumvec1
	ldcumvec = ldcumvec1
END DO

!$OMP PARALLEL DO
DO iaby = 1,naby
	ccum1(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lccumvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
	dcum1(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = ldcumvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
END DO
!$OMP END PARALLEL DO


!two quarter cumulative expected consumption
nsteps = nint(2.0/ldelta)
lccumvec = 0.0
ldcumvec = 0.0
DO it = 1,nsteps

	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		!consumption	
		lvec = lccumvec(:,iy) + ldelta*lcvec(:,iy) + ldelta*MATMUL(ymarkovoff(iy,:),TRANSPOSE(lccumvec(:,:)))
		lccumvec1(:,iy) = AUMF(iy) .umfpack. lvec
		!deposits
		lvec = ldcumvec(:,iy) + ldelta*ldvec(:,iy) + ldelta*MATMUL(ymarkovoff(iy,:),TRANSPOSE(ldcumvec(:,:)))
		ldcumvec1(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 

	lccumvec = lccumvec1
	ldcumvec = ldcumvec1
END DO


!$OMP PARALLEL DO
DO iaby = 1,naby
	ccum2(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lccumvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
	dcum2(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = ldcumvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
END DO
!$OMP END PARALLEL DO


!four quarter cumulative expected consumption
nsteps = nint(4.0/ldelta)
lccumvec = 0.0
ldcumvec = 0.0
DO it = 1,nsteps

	!$OMP PARALLEL DO PRIVATE(lvec)
	DO iy = 1,ngpy
		!consumption	
		lvec = lccumvec(:,iy) + ldelta*lcvec(:,iy) + ldelta*MATMUL(ymarkovoff(iy,:),TRANSPOSE(lccumvec(:,:)))
		lccumvec1(:,iy) = AUMF(iy) .umfpack. lvec
		!deposits
		lvec = ldcumvec(:,iy) + ldelta*ldvec(:,iy) + ldelta*MATMUL(ymarkovoff(iy,:),TRANSPOSE(ldcumvec(:,:)))
		ldcumvec1(:,iy) = AUMF(iy) .umfpack. lvec
	END DO
	!$OMP END PARALLEL DO 

	lccumvec = lccumvec1
	ldcumvec = ldcumvec1
END DO


!$OMP PARALLEL DO
DO iaby = 1,naby
	ccum4(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = lccumvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
	dcum4(afromaby(iaby),bfromaby(iaby),yfromaby(iaby)) = ldcumvec(abfromab(afromaby(iaby),bfromaby(iaby)),yfromaby(iaby))
END DO
!$OMP END PARALLEL DO



cumINITSS = CumulativePolicyType(ccum1,ccum2,ccum4,dcum1,dcum2,dcum4)

END SUBROUTINE CumulativeConsumption