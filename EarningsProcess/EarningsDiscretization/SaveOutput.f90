SUBROUTINE SaveOutput

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE


OPEN(3, FILE = trim(OutputDir) // 'y1grid.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy1,1,y1grid)
OPEN(3, FILE = trim(OutputDir) // 'y1dist.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy1,1,y1dist)
OPEN(3, FILE = trim(OutputDir) // 'y1markov.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy1,ngpy1,y1markov)
OPEN(3, FILE = trim(OutputDir) // 'y1trans.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy1,ngpy1,y1trans)
OPEN(3, FILE = trim(OutputDir) // 'y1trans_qu.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy1,ngpy1,y1trans_qu)

OPEN(3, FILE = trim(OutputDir) // 'y2grid.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy2,1,y2grid)
OPEN(3, FILE = trim(OutputDir) // 'y2dist.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy2,1,y2dist)
OPEN(3, FILE = trim(OutputDir) // 'y2markov.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy2,ngpy2,y2markov)
OPEN(3, FILE = trim(OutputDir) // 'y2trans.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy2,ngpy2,y2trans)
OPEN(3, FILE = trim(OutputDir) // 'y2trans_qu.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy2,ngpy2,y2trans_qu)


OPEN(3, FILE = trim(OutputDir) // 'ygrid_combined.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy1*ngpy2,1,ygrid_combined)
OPEN(3, FILE = trim(OutputDir) // 'ydist_combined.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy1*ngpy2,1,ydist_combined)
OPEN(3, FILE = trim(OutputDir) // 'ymarkov_combined.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy1*ngpy2,ngpy1*ngpy2,ymarkov_combined)
OPEN(3, FILE = trim(OutputDir) // 'ytrans_qu_combined.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpy1*ngpy2,ngpy1*ngpy2,ytrans_qu_combined)


IF (SaveSimulations==1) THEN
	OPEN(3, FILE = trim(OutputDir) // 'yannsim.txt', STATUS = 'replace'); CALL WriteMatrixCSV(3,nsim,5,yannsim)	
END IF

OPEN(3, FILE = trim(OutputDir) // 'parameters.txt', STATUS = 'replace')
	write(3,*) 'ngpy1 ',ngpy1
	write(3,*) 'ngpy2 ',ngpy2
	write(3,*) 'y1width ',y1width
	write(3,*) 'y2width ',y2width
	write(3,*) 'beta1 ',beta1
	write(3,*) 'sigma1 ',sigma1
	write(3,*) 'lambda1 ',lambda1
	write(3,*) 'beta2 ',beta2
	write(3,*) 'sigma2 ',sigma2
	write(3,*) 'lambda2 ',lambda2
	write(3,*) 'y1gridpar ',y1gridpar	
	write(3,*) 'y2gridpar ',y2gridpar	
	write(3,*) 'rho1 ',rho1	
	write(3,*) 'rho2 ',rho2	
	write(3,*) 'sigma1 ',sigma1	
	write(3,*) 'sigma2 ',sigma2	
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
	
CLOSE(3)

!construct combined matrix

END SUBROUTINE SaveOutput