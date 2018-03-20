SUBROUTINE Simulate

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER	:: in,it,iy
REAL(8)	:: cumtvec(Tsim),leye1(ngpy1,ngpy1),leye2(ngpy2,ngpy2)

!identity matrix
leye1 = 0.0; DO iy = 1,ngpy1
	leye1(iy,iy) = 1.0
END DO
leye2 = 0.0; DO iy = 1,ngpy2
	leye2(iy,iy) = 1.0
END DO

!discretized transition matrix
y1trans = dt*y1markov +leye1
y2trans = dt*y2markov +leye2

if(any(y1trans<0.0)) then
	write(*,*) 'y1trans<0'
	write(*,*) 'y1markov',y1markov
	write(*,*) 'y1trans',y1trans
	pause
end if
if(any(y2trans<0.0)) then
	write(*,*) 'y2trans<0'
	write(*,*) 'y2markov',y2markov
	write(*,*) 'y2trans',y2trans
	pause
end if

!quarterly transition matrix
y1trans_qu = y1trans
y2trans_qu = y2trans
DO it = 1,FLOOR(1.0/dt)-1
	y1trans_qu = MATMUL(y1trans_qu,y1trans)
	y2trans_qu = MATMUL(y2trans_qu,y2trans)
END DO	

!cumulative time vector (in quarters)
DO it = 1,Tsim
	cumtvec(it) = it*dt
END DO

!$OMP PARALLEL DO PRIVATE(it)
DO in = 1,nsim
	
	IF(Display>=2) write(*,*) 'simulating person ',in, ' of ',nsim
	
	!simulate income path in dt increments
	CALL DiscreteDist1(y1simI(in,1),ngpy1,y1dist,y1rand(in,1))
	CALL DiscreteDist1(y2simI(in,1),ngpy2,y2dist,y2rand(in,1))
	ysim(in,1) = y1grid(y1simI(in,1)) + y2grid(y2simI(in,1))
	
	DO it = 2,Tsim
		CALL DiscreteDist1(y1simI(in,it),ngpy1,y1trans(y1simI(in,it-1),:),y1rand(in,it))
		CALL DiscreteDist1(y2simI(in,it),ngpy2,y2trans(y2simI(in,it-1),:),y2rand(in,it))
		ysim(in,it) = y1grid(y1simI(in,it)) + y2grid(y2simI(in,it))
	END DO

	!aggregate to annual income
	ylevsim(in,:) = exp(ysim(in,:))
	DO it = 1,5
		yannlevsim(in,it) = SUM(ylevsim(in,:),MASK = cumtvec>4.0*real(it-1) .and. cumtvec<=4.0*real(it))
	END DO
	yannsim(in,:) = log(yannlevsim(in,:))

END DO
!$OMP END PARALLEL DO

END SUBROUTINE Simulate
