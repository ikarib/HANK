SUBROUTINE SaveOutput

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE


IF (SaveSimulations==1) THEN
	OPEN(3, FILE = trim(OutputDir) // 'yannsim.txt', STATUS = 'replace'); CALL WriteMatrixCSV(3,nsim,5,yannsim)	
END IF

OPEN(3, FILE = trim(OutputDir) // 'parameters.txt', STATUS = 'replace')
	write(3,*) 'lambda1 ',lambda1
	write(3,*) 'zeta1P ',zeta1P
	write(3,*) 'zeta1N ',zeta1N
	write(3,*) 'lambda2 ',lambda2
	write(3,*) 'zeta2P ',zeta2P
	write(3,*) 'zeta2N ',zeta2N
	write(3,*) 'rho1 ',rho1	
	write(3,*) 'rho2 ',rho2	
	write(3,*) 'sigma1 ',sigma1	
	write(3,*) 'sigma2 ',sigma2	
	write(3,*) 'delta1 ',delta1	
	write(3,*) 'delta2 ',delta2	
CLOSE(3)	

OPEN(3, FILE = trim(OutputDir) // 'moments.txt', STATUS = 'replace')
	write(3,*) 'muy ',muy
	write(3,*) 'mu2y ',mu2y
	write(3,*) 'mu3y ',mu3y
	write(3,*) 'mu4y ',mu4y
	write(3,*) 'gam3y ',gam3y
	write(3,*) 'gam4y ',gam4y
	
	write(3,*) 'muylev ',muylev
	write(3,*) 'mu2ylev ',mu2ylev
	write(3,*) 'mu3ylev ',mu3ylev
	write(3,*) 'mu4ylev ',mu4ylev
	write(3,*) 'gam3ylev ',gam3ylev
	write(3,*) 'gam4ylev ',gam4ylev
	
	write(3,*) 'mu2dy1 ',mu2dy1
	write(3,*) 'mu3dy1 ',mu3dy1
	write(3,*) 'mu4dy1 ',mu4dy1
	write(3,*) 'gam3dy1 ',gam3dy1
	write(3,*) 'gam4dy1 ',gam4dy1
	
	write(3,*) 'mu2dy5 ',mu2dy5
	write(3,*) 'mu3dy5 ',mu3dy5
	write(3,*) 'mu4dy5 ',mu4dy5
	write(3,*) 'gam3dy5 ',gam3dy5
	write(3,*) 'gam4dy5 ',gam4dy5	
	
	write(3,*) 'fracdy1less5 ',fracdy1less5
	write(3,*) 'fracdy1less10 ',fracdy1less10
	write(3,*) 'fracdy1less20 ',fracdy1less20
	write(3,*) 'fracdy1less50 ',fracdy1less50
	
CLOSE(3)

!construct combined matrix

END SUBROUTINE SaveOutput