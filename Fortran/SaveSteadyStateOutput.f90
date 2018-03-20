SUBROUTINE SaveSteadyStateOutput

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE
INTEGER		:: iy
CHARACTER	:: lstring*80


IF (Display>=1) write(*,*) 'Saving output to disk'

!grids
!OPEN(3, FILE = trim(OutputDir) // 'agrid.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,ngpa,1,agrid)
!OPEN(3, FILE = trim(OutputDir) // 'bgrid.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,ngpb,1,bgrid)
!OPEN(3, FILE = trim(OutputDir) // 'ygrid.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,ngpy,1,ygrid)
!OPEN(3, FILE = trim(OutputDir) // 'ydist.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,ngpy,1,ydist)
!OPEN(3, FILE = trim(OutputDir) // 'adelta.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,ngpa,1,adelta)
!OPEN(3, FILE = trim(OutputDir) // 'bdelta.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,ngpb,1,bdelta)

!initial steady state summary stats
OPEN(3, FILE = trim(OutputDir) // 'InitialSteadyStateParameters.txt', STATUS = 'replace')
	write(3,*) 'gam ',gam
	write(3,*) 'rho ',rho
	write(3,*) 'deathrate ',deathrate
	write(3,*) 'kappa0_w ',kappa0_w
	write(3,*) 'kappa1_w ',kappa1_w
	write(3,*) 'kappa2_w ',kappa2_w
	write(3,*) 'kappa3 ',kappa3
	write(3,*) 'kappa4_w ',kappa4_w
	write(3,*) 'kappa0_d ',kappa0_d
	write(3,*) 'kappa1_d ',kappa1_d
	write(3,*) 'kappa2_d ',kappa2_d
	write(3,*) 'kappa4_d ',kappa4_d
	write(3,*) 'corptax ',corptax

	write(3,*) 'ra ',equmINITSS%ra
	write(3,*) 'rb ',equmINITSS%rb
	write(3,*) 'rborr ',equmINITSS%rborr
	write(3,*) 'rcapital ',equmINITSS%rcapital
	write(3,*) 'wage ',equmINITSS%wage
	write(3,*) 'netwage ',equmINITSS%netwage
	write(3,*) 'bond ',equmINITSS%bond
	write(3,*) 'capital ',equmINITSS%capital
	write(3,*) 'equity ',equmINITSS%equity
	write(3,*) 'labor ',equmINITSS%labor
	write(3,*) 'output ',equmINITSS%output
	write(3,*) 'investment ',equmINITSS%investment
	write(3,*) 'govexp ',equmINITSS%govexp
	write(3,*) 'lumptransfer ',equmINITSS%lumptransfer
	write(3,*) 'labtax ',equmINITSS%labtax
	write(3,*) 'taxrev ',equmINITSS%taxrev
	write(3,*) 'govbond ',equmINITSS%govbond
	write(3,*) 'worldbond ',equmINITSS%worldbond
	write(3,*) 'fundbond ',equmINITSS%fundbond
	write(3,*) 'KYratio ',equmINITSS%KYratio
	write(3,*) 'KNratio ',equmINITSS%KNratio
	write(3,*) 'mc ',equmINITSS%mc
	write(3,*) 'rb ',equmINITSS%rb
	write(3,*) 'tfp ',equmINITSS%tfp
	write(3,*) 'pi ',equmINITSS%pi
	write(3,*) 'rnom ',equmINITSS%rnom
	write(3,*) 'priceadjust ',equmINITSS%priceadjust
	write(3,*) 'profit ',equmINITSS%profit
	write(3,*) 'dividend ',equmINITSS%dividend
	write(3,*) 'divrate ',equmINITSS%divrate

	write(3,*) 'borrwedge ',equmINITSS%borrwedge
	write(3,*) 'rho ',equmINITSS%rho
	write(3,*) 'fundlev ',equmINITSS%fundlev
	write(3,*) 'deprec ',deprec
	
	write(3,*) 'Ea ',statsINITSS%Ea
	write(3,*) 'Eb ',statsINITSS%Eb
	write(3,*) 'Ec ',statsINITSS%Ec
	write(3,*) 'Ehours ',statsINITSS%Ehours
	write(3,*) 'Elabor ',statsINITSS%Elabor
	write(3,*) 'Ed ',statsINITSS%Ed
	write(3,*) 'Ewage ',statsINITSS%Ewage
	write(3,*) 'Enetlabinc ',statsINITSS%Enetlabinc
	write(3,*) 'Egrosslabinc ',statsINITSS%Egrosslabinc
	write(3,*) 'Enetprofinc ',statsINITSS%Enetprofinc
	write(3,*) 'Egrossprofinc ',statsINITSS%Egrossprofinc
	write(3,*) 'Einc ',statsINITSS%Einc
	write(3,*) 'Enw ',statsINITSS%Enw
	write(3,*) 'FRACa0 ',statsINITSS%FRACa0
	write(3,*) 'FRACb0 ',statsINITSS%FRACb0
	write(3,*) 'FRACb0a0 ',statsINITSS%FRACb0a0
	write(3,*) 'FRACnw0 ',statsINITSS%FRACnw0
	write(3,*) 'FRACb0aP ',statsINITSS%FRACb0aP
	write(3,*) 'FRACa0close ',statsINITSS%FRACa0close
	write(3,*) 'FRACb0close ',statsINITSS%FRACb0close
	write(3,*) 'FRACb0a0close ',statsINITSS%FRACb0a0close
	write(3,*) 'FRACnw0close ',statsINITSS%FRACnw0close
	write(3,*) 'FRACbN ',statsINITSS%FRACbN
	write(3,*) 'EbN ',statsINITSS%EbN
	write(3,*) 'EbP ',statsINITSS%EbP
	write(3,*) 'Eadjcost ',statsINITSS%Eadjcost
	write(3,*) 'GINIa ',statsINITSS%GINIa
	write(3,*) 'GINIb ',statsINITSS%GINIb
	write(3,*) 'GINInw ',statsINITSS%GINInw
	write(3,*) 'GINIc ',statsINITSS%GINIc
	write(3,*) 'GINIinc ',statsINITSS%GINIinc
	write(3,*) 'Ec_bN ',statsINITSS%Ec_bN
	write(3,*) 'Ec_b0close ',statsINITSS%Ec_b0close
	write(3,*) 'Ec_b0far ',statsINITSS%Ec_b0far
	


CLOSE(3)

!OPEN(3, FILE = trim(OutputDir) // 'c.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,reshape(c,(/nab,ngpy/))) 
!OPEN(3, FILE = trim(OutputDir) // 'h.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,reshape(h,(/nab,ngpy/))) 
!OPEN(3, FILE = trim(OutputDir) // 's.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,reshape(s,(/nab,ngpy/))) 
!OPEN(3, FILE = trim(OutputDir) // 'd.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,reshape(d,(/nab,ngpy/))) 
!OPEN(3, FILE = trim(OutputDir) // 'u.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,reshape(u,(/nab,ngpy/))) 
!OPEN(3, FILE = trim(OutputDir) // 'bdot.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,reshape(bdot,(/nab,ngpy/))) 
!OPEN(3, FILE = trim(OutputDir) // 'V.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,reshape(V,(/nab,ngpy/))) 
!OPEN(3, FILE = trim(OutputDir) // 'Vnew.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,reshape(Vnew,(/nab,ngpy/))) 
!OPEN(3, FILE = trim(OutputDir) // 'ACOOrow.txt', STATUS = 'replace'); 
!OPEN(4, FILE = trim(OutputDir) // 'ACOOcol.txt', STATUS = 'replace'); 
!OPEN(5, FILE = trim(OutputDir) // 'ACOOval.txt', STATUS = 'replace'); 
!OPEN(6, FILE = trim(OutputDir) // 'AUCOOrow.txt', STATUS = 'replace'); 
!OPEN(7, FILE = trim(OutputDir) // 'AUCOOcol.txt', STATUS = 'replace'); 
!OPEN(8, FILE = trim(OutputDir) // 'AUCOOval.txt', STATUS = 'replace'); 
!OPEN(9, FILE = trim(OutputDir) // 'BCSRrow.txt', STATUS = 'replace'); 
!OPEN(10, FILE = trim(OutputDir) // 'BCSRcol.txt', STATUS = 'replace'); 
!OPEN(11, FILE = trim(OutputDir) // 'BCSRval.txt', STATUS = 'replace'); 
!DO iy = 1,ngpy
!	WRITE(3,'(10000G26.15)') ACOO(iy)%row
!	WRITE(4,'(10000G26.15)') ACOO(iy)%col
!	WRITE(5,'(10000G26.15)') ACOO(iy)%val
!	WRITE(6,'(6000G26.15)') AUCOO(iy)%row
!	WRITE(7,'(6000G26.15)') AUCOO(iy)%col
!	WRITE(8,'(6000G26.15)') AUCOO(iy)%val
!	WRITE(9,'(2001G26.15)') BCSR(iy)%row
!	WRITE(10,'(10000G26.15)') BCSR(iy)%col
!	WRITE(11,'(10000G26.15)') BCSR(iy)%val
!END DO
!CLOSE(3)
!CLOSE(4)
!CLOSE(5)
!CLOSE(6)
!CLOSE(7)
!CLOSE(8)
!CLOSE(9)
!CLOSE(10)
!CLOSE(11)
!OPEN(3, FILE = trim(OutputDir) // 'gmat.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,reshape(gmat,(/nab,ngpy/))) 

!OPEN(3, FILE = trim(OutputDir) // 'lccumvec.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,lccumvec) 
!OPEN(3, FILE = trim(OutputDir) // 'ldcumvec.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,ldcumvec) 
!CALL EXIT(0)

!OPEN(1, FILE = trim(OutputDir) // 'AUCOOrow.txt', STATUS = 'replace'); 
!OPEN(2, FILE = trim(OutputDir) // 'AUCOOcol.txt', STATUS = 'replace'); 
!OPEN(3, FILE = trim(OutputDir) // 'AUCOOval.txt', STATUS = 'replace'); 
!DO iy = 1,ngpy
!	AUCOO(iy)%row = 0
!	AUCOO(iy)%col = 0
!	AUCOO(iy)%val = 0.0_8
!	AUCOO(iy)%nz = 0
!	CALL csrcoo (nab, 3, AUCSR(iy)%nz, AUCSR(iy)%val, AUCSR(iy)%col, AUCSR(iy)%row, AUCOO(iy)%nz, AUCOO(iy)%val, AUCOO(iy)%row, AUCOO(iy)%col, ierr )
!	IF (ierr) THEN
!		WRITE(*,*) "error=", ierr
!	END IF
!	WRITE(1,'(6000G26.15)') AUCOO(iy)%row
!	WRITE(2,'(6000G26.15)') AUCOO(iy)%col
!	WRITE(3,'(6000G26.15)') AUCOO(iy)%val
!END DO
!CLOSE(1)
!CLOSE(2)
!CLOSE(3)

!OPEN(3, FILE = trim(OutputDir) // 'lsubeff1ass.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,lsubeff1ass) 
!OPEN(3, FILE = trim(OutputDir) // 'lsubeff2ass.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,lsubeff2ass) 
!OPEN(3, FILE = trim(OutputDir) // 'lwealtheff1ass.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,lwealtheff1ass) 
!OPEN(3, FILE = trim(OutputDir) // 'lwealtheff2ass.txt', STATUS = 'replace'); 
!CALL WriteMatrix(3,nab,ngpy,lwealtheff2ass) 


!initial steady state distributions and policy functions
CALL system ("mkdir -p " // trim(OutputDir) // "INITSS")

OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'PERCa.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,11,statsINITSS%PERCa)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'PERCb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,11,statsINITSS%PERCb)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'PERCnw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,11,statsINITSS%PERCnw)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'PERCc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,11,statsINITSS%PERCc)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'PERCinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,11,statsINITSS%PERCinc)

OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'Ea_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,4,statsINITSS%Ea_nwQ)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'Eb_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,4,statsINITSS%Eb_nwQ)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'Ec_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,4,statsINITSS%Ec_nwQ)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'Einc_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,4,statsINITSS%Einc_nwQ)

OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'Ea_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,4,statsINITSS%Ea_incQ)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'Eb_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,4,statsINITSS%Eb_incQ)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'Ec_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,4,statsINITSS%Ec_incQ)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'Einc_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,1,4,statsINITSS%Einc_incQ)

DO iy = 1,ngpy
	IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
	IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy

	OPEN(3, FILE = trim(OutputDir) // "INITSS/" // 'V_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%V(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'dep_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%d(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'con_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%c(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'sav_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%s(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'hour_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%h(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'bdot_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%bdot(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'gjoint_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%gjoint(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/" // 'B_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%B(iy)%val)

	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'ccum1_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,cumINITSS%ccum1(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'ccum2_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,cumINITSS%ccum2(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'ccum4_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,cumINITSS%ccum4(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'dcum1_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,cumINITSS%dcum1(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'dcum2_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,cumINITSS%dcum2(:,:,iy))
	OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'dcum4_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,cumINITSS%dcum4(:,:,iy))

	IF(ComputeDiscountedMPC==1) THEN
		OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'mpc_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%mpc(:,:,iy))
		OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'subeff1ass_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%subeff1ass(:,:,iy))
		OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'subeff2ass_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%subeff2ass(:,:,iy))
		OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'wealtheff1ass_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%wealtheff1ass(:,:,iy))
		OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'wealtheff2ass_INITSS_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,solnINITSS%wealtheff2ass(:,:,iy))
	END IF	

END DO

OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'gamarg_INITSS.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpy,solnINITSS%gamarg)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'gbmarg_INITSS.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpb,ngpy,solnINITSS%gbmarg)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'gabmarg_INITSS.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,gabmarg)
OPEN(3, FILE = trim(OutputDir) // "INITSS/"  // 'gabcum_INITSS.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,gabcum)


END SUBROUTINE SaveSteadyStateOutput
