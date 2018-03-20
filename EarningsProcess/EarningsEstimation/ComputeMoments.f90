SUBROUTINE ComputeMoments

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

!central moments: logs
muy = SUM(yannsim(:,1))/real(nsim)
mu2y = SUM((yannsim(:,1)-muy)**2)/real(nsim)
mu3y = SUM((yannsim(:,1)-muy)**3)/real(nsim)
mu4y = SUM((yannsim(:,1)-muy)**4)/real(nsim)

!standardised moments: logs
IF (mu2y>0.0) THEN
	gam3y = mu3y/(mu2y**1.5)
	gam4y = mu4y/(mu2y**2)
ELSE
	gam3y = 0.0
	gam4y = 0.0
END IF

!central moments: logs
muylev = SUM(yannlevsim(:,1))/real(nsim)
mu2ylev = SUM((yannlevsim(:,1)-muylev)**2)/real(nsim)
mu3ylev = SUM((yannlevsim(:,1)-muylev)**3)/real(nsim)
mu4ylev = SUM((yannlevsim(:,1)-muylev)**4)/real(nsim)

!standardised moments: logs
IF (mu2ylev>0.0) THEN
	gam3ylev = mu3ylev/(mu2ylev**1.5)
	gam4ylev = mu4ylev/(mu2ylev**2)
ELSE
	gam3ylev = 0.0
	gam4ylev = 0.0
END IF

!central moments: 1 year log changes
mudy1 = SUM(yannsim(:,2)-yannsim(:,1))/real(nsim)
mu2dy1 = SUM((yannsim(:,2)-yannsim(:,1)-mudy1)**2)/real(nsim)
mu3dy1 = SUM((yannsim(:,2)-yannsim(:,1)-mudy1)**3)/real(nsim)
mu4dy1 = SUM((yannsim(:,2)-yannsim(:,1)-mudy1)**4)/real(nsim)

!standardised moments: 1 year log changes
IF (mu2dy1>0.0) THEN
	gam3dy1 = mu3dy1/(mu2dy1**1.5)
	gam4dy1 = mu4dy1/(mu2dy1**2)
ELSE
	gam3dy1 = 0.0
	gam4dy1 = 0.0
END IF

!central moments: 5 year log changes
mudy5 = SUM(yannsim(:,5)-yannsim(:,1))/real(nsim)
mu2dy5 = SUM((yannsim(:,5)-yannsim(:,1)-mudy5)**2)/real(nsim)
mu3dy5 = SUM((yannsim(:,5)-yannsim(:,1)-mudy5)**3)/real(nsim)
mu4dy5 = SUM((yannsim(:,5)-yannsim(:,1)-mudy5)**4)/real(nsim)

!standardised moments: 5 year log changes
IF (mu2dy5>0.0) THEN
	gam3dy5 = mu3dy5/(mu2dy5**1.5)
	gam4dy5 = mu4dy5/(mu2dy5**2)
ELSE
	gam3dy5 = 0.0
	gam4dy5 = 0.0
END IF

!fraction 1 year log changes in ranges
fracdy1less5 = COUNT(ABS(yannsim(:,2)-yannsim(:,1)) < 0.05)/real(nsim)
fracdy1less10 = COUNT(ABS(yannsim(:,2)-yannsim(:,1)) < 0.1)/real(nsim)
fracdy1less20 = COUNT(ABS(yannsim(:,2)-yannsim(:,1)) < 0.2)/real(nsim)
fracdy1less50 = COUNT(ABS(yannsim(:,2)-yannsim(:,1)) < 0.5)/real(nsim)



END SUBROUTINE ComputeMoments