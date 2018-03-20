SUBROUTINE TransitionMatrix(lparam,ly1grid,ly1markov,ly1dist,ly2grid,ly2markov,ly2dist)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN)			:: lparam(4)
REAL(8), INTENT(OUT)		:: ly1grid(ngpy1),ly1markov(ngpy1,ngpy1),ly1dist(ngpy1)
REAL(8), INTENT(OUT)		:: ly2grid(ngpy2),ly2markov(ngpy2,ngpy2),ly2dist(ngpy2)
INTEGER 	:: iy,iy1,iy2,it,ierr,ij,ii(2)
REAL(8) 	:: lerr,ldelta,lp(2)
REAL(8), DIMENSION(ngpy1,ngpy1) 	:: leye1,ljump1,ldrift1,lmat1,lmatinv1
REAL(8), DIMENSION(ngpy2,ngpy2) 	:: leye2,ljump2,ldrift2,lmat2,lmatinv2
REAL(8), DIMENSION(ngpy1) 	:: ly1dist_update
REAL(8), DIMENSION(ngpy2) 	:: ly2dist_update
REAL(8), DIMENSION(ngpy1-1) 	:: ldy1grid
REAL(8), DIMENSION(ngpy2-1) 	:: ldy2grid

!identity matrix
leye1 = 0.0; DO iy = 1,ngpy1
	leye1(iy,iy) = 1.0
END DO
leye2 = 0.0; DO iy = 1,ngpy2
	leye2(iy,iy) = 1.0
END DO

!extract parametersye
y1width = lparam(1)
y2width = lparam(2)
y1gridpar = lparam(3)
y2gridpar = lparam(4)


!construct grid: y1
ly1grid(1) = -y1width/2.0
ly1grid(ngpy1) = y1width/2.0
CALL PowerSpacedGrid ((ngpy1+1)/2,y1gridpar,0.0_8,ly1grid(ngpy1),ly1grid((ngpy1+1)/2:ngpy1))
DO iy = 1,(ngpy1+1)/2-1
	ly1grid(iy) = ly1grid(1) + (ly1grid(ngpy1) - ly1grid(ngpy1+1-iy))
END DO
ldy1grid = ly1grid(2:ngpy1)-ly1grid(1:ngpy1-1)

!construct grid: y2
ly2grid(1) = -y2width/2.0
ly2grid(ngpy2) = y2width/2.0
CALL PowerSpacedGrid ((ngpy2+1)/2,y2gridpar,0.0_8,ly2grid(ngpy2),ly2grid((ngpy2+1)/2:ngpy2))
DO iy = 1,(ngpy2+1)/2-1
	ly2grid(iy) = ly2grid(1) + (ly2grid(ngpy2) - ly2grid(ngpy2+1-iy))
END DO
! CALL SymmetricPowerSpacedGrid (ngpy2,y2gridpar,y2width/2.0,0.0_8,ly2grid)
ldy2grid = ly2grid(2:ngpy2)-ly2grid(1:ngpy2-1)
	
IF(AddPointsCloseZeroY2==1) THEN
	! replace bottom and top 2 points with unweighted mean
	ly2grid(1) = SUM(ly2grid(1:2))/2.0
	ly2grid(ngpy2) = SUM(ly2grid(ngpy2-1:ngpy2))/2.0
	
	!move points out
	IF(ngpy2>= 7) THEN
		ly2grid(2:(ngpy2+1)/2-2) = ly2grid(3:(ngpy2+1)/2-1)
		ly2grid((ngpy2+1)/2+2:ngpy2-1) = ly2grid((ngpy2+1)/2+1:ngpy2-2)
	END IF

	! add middle points
	ly2grid((ngpy2+1)/2-1) = SUM(ly2grid((ngpy2+1)/2-1:(ngpy2+1)/2))/2.0
	ly2grid((ngpy2+1)/2+1) = SUM(ly2grid((ngpy2+1)/2:(ngpy2+1)/2+1))/2.0
	
	ldy2grid = ly2grid(2:ngpy2)-ly2grid(1:ngpy2-1)
	
END IF

!construct matrices for jump processes
DO iy = 1,ngpy1
	DO iy2 = 1,ngpy1
		IF (RestrictJumpDomain==0) THEN
			IF(iy2==1) THEN		
				ljump1(iy,iy2) = NormalCDF(y1grid(iy2)+0.5*ldy1grid(iy2),rho1*y1grid(iy),sigma1)
			ELSE IF(iy2>1 .and. iy2<ngpy1) THEN
				ljump1(iy,iy2) = NormalCDF(y1grid(iy2)+0.5*ldy1grid(iy2),rho1*y1grid(iy),sigma1) - NormalCDF(y1grid(iy2)-0.5*ldy1grid(iy2-1),rho1*y1grid(iy),sigma1)
			ELSE IF(iy2==ngpy1) THEN
				ljump1(iy,iy2) = 1.0 - NormalCDF(y1grid(iy2)-0.5*ldy1grid(iy2-1),rho1*y1grid(iy),sigma1)
			END IF
		ELSE IF (RestrictJumpDomain==1) THEN
			IF(iy2==1) THEN		
				ljump1(iy,iy2) = TruncatedNormalCDF(y1grid(iy2)+0.5*ldy1grid(iy2),rho1*y1grid(iy),sigma1,y1grid(1),y1grid(ngpy1))
			ELSE IF(iy2>1 .and. iy2<ngpy1) THEN
				ljump1(iy,iy2) = TruncatedNormalCDF(y1grid(iy2)+0.5*ldy1grid(iy2),rho1*y1grid(iy),sigma1,y1grid(1),y1grid(ngpy1)) - TruncatedNormalCDF(y1grid(iy2)-0.5*ldy1grid(iy2-1),rho1*y1grid(iy),sigma1,y1grid(1),y1grid(ngpy1))
			ELSE IF(iy2==ngpy1) THEN
				ljump1(iy,iy2) = 1.0 - TruncatedNormalCDF(y1grid(iy2)-0.5*ldy1grid(iy2-1),rho1*y1grid(iy),sigma1,y1grid(1),y1grid(ngpy1))
			END IF
		END IF
			
	END DO
	ljump1(iy,:) = ljump1(iy,:) / SUM(ljump1(iy,:))
END DO
ljump1 = lambda1*(ljump1-leye1)

DO iy = 1,ngpy2
	DO iy2 = 1,ngpy2
		IF (RestrictJumpDomain==0) THEN
			IF(iy2==1) THEN		
				ljump2(iy,iy2) = NormalCDF(y2grid(iy2)+0.5*ldy2grid(iy2),rho2*y2grid(iy),sigma2)
			ELSE IF(iy2>1 .and. iy2<ngpy2) THEN
				ljump2(iy,iy2) = NormalCDF(y2grid(iy2)+0.5*ldy2grid(iy2),rho2*y2grid(iy),sigma2) - NormalCDF(y2grid(iy2)-0.5*ldy2grid(iy2-1),rho2*y2grid(iy),sigma2)
			ELSE IF(iy2==ngpy2) THEN
				ljump2(iy,iy2) = 1.0 - NormalCDF(y2grid(iy2)-0.5*ldy2grid(iy2-1),rho2*y2grid(iy),sigma2)
			END IF
		ELSE IF (RestrictJumpDomain==1) THEN
			IF(iy2==1) THEN		
				ljump2(iy,iy2) = TruncatedNormalCDF(y2grid(iy2)+0.5*ldy2grid(iy2),rho2*y2grid(iy),sigma2,y2grid(1),y2grid(ngpy2))
			ELSE IF(iy2>1 .and. iy2<ngpy2) THEN
				ljump2(iy,iy2) = TruncatedNormalCDF(y2grid(iy2)+0.5*ldy2grid(iy2),rho2*y2grid(iy),sigma2,y2grid(1),y2grid(ngpy2)) - TruncatedNormalCDF(y2grid(iy2)-0.5*ldy2grid(iy2-1),rho2*y2grid(iy),sigma2,y2grid(1),y2grid(ngpy2))
			ELSE IF(iy2==ngpy2) THEN
				ljump2(iy,iy2) = 1.0 - TruncatedNormalCDF(y2grid(iy2)-0.5*ldy2grid(iy2-1),rho2*y2grid(iy),sigma2,y2grid(1),y2grid(ngpy2))
			END IF
		END IF
		
	END DO
	ljump2(iy,:) = ljump2(iy,:) / SUM(ljump2(iy,:))
END DO
ljump2 = lambda2*(ljump2-leye2)

!construct matric for drift process
ldrift1 = 0.0
DO iy1 = 1,ngpy1 !rows

	IF (y1grid(iy1) .ne. 0.0) THEN		
		CALL FindLinProb1 ((1.0-beta1*deltaforapprox)*y1grid(iy1),y1grid,ii,lp)
		IF(DriftPointsVersion==2) THEN
			ldrift1(iy1,iy1) = -1.0
			IF (y1grid(iy1)<0.0) THEN
				ldrift1(iy1,iy1) = ldrift1(iy1,iy1) + (y1grid(ii(2)) - (1.0-beta1)*y1grid(iy1))/(y1grid(ii(2))-y1grid(iy1))
				ldrift1(iy1,ii(2)) = (-y1grid(iy1) + (1.0-beta1)*y1grid(iy1))/(y1grid(ii(2))-y1grid(iy1))
			ELSE IF (y1grid(iy1)>0.0) THEN
				ldrift1(iy1,ii(1)) = (y1grid(iy1) - (1.0-beta1)*y1grid(iy1))/(y1grid(iy1)-y1grid(ii(1)))
				ldrift1(iy1,iy1) = ldrift1(iy1,iy1) + (-y1grid(ii(1)) + (1.0-beta1)*y1grid(iy1))/(y1grid(iy1)-y1grid(ii(1)))
			END IF
		ELSE IF(DriftPointsVersion==3) THEN
			ldrift1(iy1,iy1) = -1.0
			ldrift1(iy1,ii(1)) = ldrift1(iy1,ii(1)) + (y1grid(ii(2)) - (1.0-beta1)*y1grid(iy1))/(y1grid(ii(2))-y1grid(ii(1)))
			ldrift1(iy1,ii(2)) = ldrift1(iy1,ii(2)) + (-y1grid(ii(1)) + (1.0-beta1)*y1grid(iy1))/(y1grid(ii(2))-y1grid(ii(1)))
		END IF
	END IF
END DO	
! 
! DO iy1 = 1,ngpy1 !rows
! 	write(*,*) ldrift1(iy1,:)
! END DO
! stop

ldrift2 = 0.0
DO iy1 = 1,ngpy2 !rows

	IF (y2grid(iy1) .ne. 0.0) THEN		
		CALL FindLinProb1 ((1.0-beta2*deltaforapprox)*y2grid(iy1),y2grid,ii,lp)
		IF(DriftPointsVersion==2) THEN
			ldrift2(iy1,iy1) = -1.0
			IF (y2grid(iy1)<0.0) THEN
				ldrift2(iy1,iy1) = ldrift2(iy1,iy1) + (y2grid(ii(2)) - (1.0-beta2)*y2grid(iy1))/(y2grid(ii(2))-y2grid(iy1))
				ldrift2(iy1,ii(2)) = (-y2grid(iy1) + (1.0-beta2)*y2grid(iy1))/(y2grid(ii(2))-y2grid(iy1))
			ELSE IF (y2grid(iy1)>0.0) THEN
				ldrift2(iy1,ii(1)) = (y2grid(iy1) - (1.0-beta2)*y2grid(iy1))/(y2grid(iy1)-y2grid(ii(1)))
				ldrift2(iy1,iy1) = ldrift2(iy1,iy1) + (-y2grid(ii(1)) + (1.0-beta2)*y2grid(iy1))/(y2grid(iy1)-y2grid(ii(1)))
			END IF
		ELSE IF(DriftPointsVersion==3) THEN
			ldrift2(iy1,iy1) = -1.0
			ldrift2(iy1,ii(1)) = ldrift2(iy1,ii(1)) + (y2grid(ii(2)) - (1.0-beta2)*y2grid(iy1))/(y2grid(ii(2))-y2grid(ii(1)))
			ldrift2(iy1,ii(2)) = ldrift2(iy1,ii(2)) + (-y2grid(ii(1)) + (1.0-beta2)*y2grid(iy1))/(y2grid(ii(2))-y2grid(ii(1)))
		END IF
	END IF
END DO

!total transistion matrix
ly1markov = ljump1 + ldrift1
ly2markov = ljump2 + ldrift2



!compute ergodic distribution using implicit method
ly1dist = 0.0
ly1dist((ngpy1+1)/2) = 1.0

ldelta = 1.0 !1.0e6
lmat1 = leye1 - ldelta*ly1markov
CALL InvertMatrix(lmat1, lmatinv1, ngpy1, ierr)
it = 1
lerr = 1.0
DO WHILE (lerr>1.0e-14 .and. it<1000)
	ly1dist_update = MATMUL(ly1dist,lmatinv1)
	
	ly1dist_update = MERGE(0.0_8,ly1dist_update,abs(ly1dist_update)<1.0e-20_8)
	ly1dist_update = ly1dist_update/SUM(ly1dist_update)
	
	lerr = MAXVAL(ABS(ly1dist_update-ly1dist))
	ly1dist = ly1dist_update
	it = it+1	
END DO

IF(Display>=1) write(*,*) 'stationary dist: ', it, ' iters, error ',lerr
IF(Display>=1) write(*,*) 'sum y1dist', SUM(ly1dist)

ly2dist = 0.0
ly2dist((ngpy2+1)/2) = 1.0

ldelta = 1.0 !1.0e6
lmat2 = leye2 - ldelta*ly2markov
CALL InvertMatrix(lmat2, lmatinv2, ngpy2, ierr)
it = 1
lerr = 1.0
DO WHILE (lerr>1.0e-14 .and. it<1000)
	ly2dist_update = MATMUL(ly2dist,lmatinv2)
	
	ly2dist_update = MERGE(0.0_8,ly2dist_update,abs(ly2dist_update)<1.0e-20_8)
	ly2dist_update = ly2dist_update/SUM(ly2dist_update)
	
	lerr = MAXVAL(ABS(ly2dist_update-ly2dist))
	ly2dist = ly2dist_update
	it = it+1	
END DO

IF(Display>=1) write(*,*) 'stationary dist: ', it, ' iters, error ',lerr
IF(Display>=1) write(*,*) 'sum y2dist', SUM(ly2dist)
IF(Display>=2) THEN
	DO iy = 1,ngpy1
		write(*,*) 'sum y1markov row ', SUM(ly1markov(iy,:))
	END DO

	DO iy = 1,ngpy2
		write(*,*) 'sum y2markov row ', SUM(ly2markov(iy,:))
	END DO
END IF

! STOP
END SUBROUTINE TransitionMatrix