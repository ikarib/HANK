SUBROUTINE Simulate

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER			:: in,it,it1,itN
REAL(8)			:: lssvar1,lssvar2

IF(rho1 .ne. 1.0) THEN
	IF(delta1 == 0.0) lssvar1 = (sigma1**2) / (1.0 -rho1**2)
	IF(delta1 .ne. 0.0) lssvar1 = lambda1*(sigma1**2) / (2.0*delta1 + lambda1*(1.0 -rho1**2))
ELSE IF (rho1==1.0) THEN	
	lssvar1 = (sigma1**2) / (1.0 -0.99**2)
END IF
	
IF(Include2ndProcess==1) THEN
	IF(delta2 == 0.0) lssvar2 = (sigma2**2) / (1.0 -rho2**2)
	IF(delta2 .ne. 0.0) lssvar2 = lambda2*(sigma2**2) / (2.0*delta2 + lambda2*(1.0 -rho2**2))
ELSE
	lssvar2 = 0.0
END IF

!$OMP PARALLEL DO PRIVATE(it)
DO in = 1,nsim
	
	IF(Display==1) write(*,*) 'simulating person ',in, ' of ',nsim
	
	!draw initial from normal distribution with same mean and variance
	IF(UseNormalDist==1) THEN
		y1sim(in,1) = sqrt(lssvar1)*y1rand(in,1)
		y2sim(in,1) = sqrt(lssvar2)*y2rand(in,1)
	ELSE IF(UseDoubleParetoDist==1) THEN
		y1sim(in,1) = DoubleParetoInverseCDF(y1rand(in,1),0.0_8,zeta1P) 
		y2sim(in,1) = DoubleParetoInverseCDF(y2rand(in,1),0.0_8,zeta2P) 
	END IF	
	ysim(in,1) = y1sim(in,1) + y2sim(in,1)
	
	!simulate income path in dt increments
	DO it = 1,Tsim-1
		CALL DiscreteDist1(y1jumpI(in,it),2,(/1.0-dt*lambda1, dt*lambda1/),y1jumprand(in,it))
		y1jumpI(in,it) = y1jumpI(in,it)-1
		IF(y1jumpI(in,it)==1) THEN
			IF(UseNormalDist==1) y1sim(in,it+1) = rho1*y1sim(in,it) + sigma1*y1rand(in,it+1)
			IF(UseDoubleParetoDist==1) y1sim(in,it+1) = DoubleParetoInverseCDF(y1rand(in,it+1),rho1*y1sim(in,it) ,zeta1P) 
		ELSE
			IF(AdditiveDrift==0) y1sim(in,it+1) = (1.0 - dt*delta1)*y1sim(in,it)
			IF(AdditiveDrift==1) y1sim(in,it+1) = y1sim(in,it) - dt*delta1*abs(y1sim(in,it))/y1sim(in,it)
		END IF
		
		IF(Include2ndProcess==1) THEN
		
			CALL DiscreteDist1(y2jumpI(in,it),2,(/1.0-dt*lambda2, dt*lambda2/),y2jumprand(in,it))
			y2jumpI(in,it) = y2jumpI(in,it)-1
			IF(y2jumpI(in,it)==1) THEN
				IF(UseNormalDist==1) y2sim(in,it+1) = rho2*y2sim(in,it) + sigma2*y2rand(in,it+1)
				IF(UseDoubleParetoDist==1) y2sim(in,it+1) = DoubleParetoInverseCDF(y2rand(in,it+1),rho2*y2sim(in,it) ,zeta2P) 
			ELSE
				IF(AdditiveDrift==0) y2sim(in,it+1) = (1.0 - dt*delta2)*y2sim(in,it)
				IF(AdditiveDrift==1) y2sim(in,it+1) = y2sim(in,it) - dt*delta2*abs(y2sim(in,it))/y2sim(in,it)
			
			END IF
		ELSE
			y2sim(in,it+1) = 0.0
		END IF
			
		ysim(in,it+1) = y1sim(in,it+1) + y2sim(in,it+1)
	END DO
	ylevsim(in,:) = exp(ysim(in,:))

	!aggregate to annual income
	DO it = 1,Tann
		it1 = Tburn + 1 + FLOOR(4.0/dt)*(it-1)
		itN = it1-1 + FLOOR(4.0/dt)
		yannlevsim(in,it) = SUM(ylevsim(in,it1:itN))
	END DO
	yannsim(in,:) = log(yannlevsim(in,:))

END DO
!$OMP END PARALLEL DO

END SUBROUTINE Simulate
