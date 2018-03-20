SUBROUTINE CombinedProcess

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER	:: iy,iy1,iy2,iyn,iyn1,iyn2,ii,imin(1),it,ierr
INTEGER :: iorder(ngpy1*ngpy2)
REAL(8), DIMENSION(ngpy1*ngpy2) :: lygrid_combined,lydist_combined,lydist_combined_update
REAL(8), DIMENSION(ngpy1*ngpy2,ngpy1*ngpy2) :: lytrans_qu_combined,lymarkov_combined,lmat,lmatinv,leye
REAL(8) 	:: lbeta,lerr


iy = 0
DO iy1 = 1,ngpy1
	DO iy2 = 1,ngpy2
		iy = iy +1
		lygrid_combined(iy) = y1grid(iy1) + y2grid(iy2)
		lydist_combined(iy) = y1dist(iy1)*y2dist(iy2)
		
		iyn = 0
		DO iyn1 = 1,ngpy1
			DO iyn2 = 1,ngpy2
				iyn = iyn +1
		
 				lytrans_qu_combined(iy,iyn) = y1trans_qu(iy1,iyn1)*y2trans_qu(iy2,iyn2)
				IF(iy1==iyn1 .and. iy2==iyn2) lymarkov_combined(iy,iyn) = y1markov(iy1,iyn1) + y2markov(iy2,iyn2)
				IF(iy1==iyn1 .and. iy2.ne.iyn2) lymarkov_combined(iy,iyn) = y2markov(iy2,iyn2) 
				IF(iy1.ne.iyn1 .and. iy2==iyn2) lymarkov_combined(iy,iyn) = y1markov(iy1,iyn1)
				IF(iy1.ne.iyn1 .and. iy2.ne.iyn2) lymarkov_combined(iy,iyn) = 0.0
			END DO
		END DO
	END DO
END DO


!sort into ascending order
imin = MINLOC(ygrid_combined)
iorder(1) = imin(1)

DO ii = 2,ngpy1*ngpy2
	imin = MINLOC(lygrid_combined, MASK = lygrid_combined>lygrid_combined(iorder(ii-1)))
	iorder(ii) = imin(1)
END DO

DO iy = 1,ngpy1*ngpy2
	ygrid_combined(iy) = lygrid_combined(iorder(iy))
	ydist_combined(iy) = lydist_combined(iorder(iy))
	DO iy2 = 1,ngpy1*ngpy2
		ytrans_qu_combined(iy,iy2) = lytrans_qu_combined(iorder(iy),iorder(iy2))
		ymarkov_combined(iy,iy2) = lymarkov_combined(iorder(iy),iorder(iy2))
	END DO
END DO

! 
! !compute ergodic distribution
! lydist_combined = 0.1
! ii = ((ngpy1+1)/2-1)*ngpy2 + (ngpy2+1)/2
! lydist_combined(ii) = 0.9
! 
! lbeta = 1.0e2
! !identity matrix
! leye = 0.0; DO iy = 1,ngpy1*ngpy2
! 	leye(iy,iy) = 1.0
! END DO
! 
! lmat = leye - lbeta*ymarkov_combined
! CALL InvertMatrix(lmat, lmatinv, ngpy1*ngpy2, ierr)
! it = 1
! lerr = 1.0
! DO WHILE (lerr>1.0e-14 .and. it<1000)
! 	lydist_combined_update = MATMUL(lydist_combined,lmatinv)
! 	
! 	lydist_combined_update = MERGE(0.0_8,lydist_combined_update,abs(lydist_combined_update)<1.0e-20_8)
! 	lydist_combined_update = lydist_combined_update/SUM(lydist_combined_update)
! 	
! 	lerr = MAXVAL(ABS(lydist_combined_update-lydist_combined))
! 	lydist_combined = lydist_combined_update
! 	write(*,*) 'it ',it,' ', lerr
! 	it = it+1	
! END DO
! ydist_combined = lydist_combined


END SUBROUTINE CombinedProcess





