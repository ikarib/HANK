SUBROUTINE SaveIRFOutput

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE
INTEGER		:: iy,it,ip,iq
CHARACTER	:: lstring*80,lt*80
REAL(8)		:: lmat(Ttransition,11),lmat4(Ttransition,4)

!sticky price transition
IF(SolveStickyPriceTransition==1) THEN
	CALL system ("mkdir -p " // trim(OutputDirIRF) // "STICKY")	

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'ra.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%ra)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'rcapital.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%rcapital)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'wage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%wage)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'netwage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%netwage)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'KYratio.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%KYratio)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'KNratio.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%KNratio)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'mc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%mc)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'rb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%rb)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'rborr.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%rborr)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'tfp.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%tfp)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'pi.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%pi)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'rnom.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%rnom)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'gap.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%gap)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'capital.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%capital)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'bond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%bond)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'labor.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%labor)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'output.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%output)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'investment.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%investment)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'taxrev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%taxrev)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'govexp.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%govexp)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'govbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%govbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'worldbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%worldbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'borrwedge.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%borrwedge)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'rho.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%rho)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'mpshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%mpshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'prefshock.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%prefshock)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'elast.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%elast)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'fundlev.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%fundlev)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'gam.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%gam)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'fundbond.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%fundbond)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'priceadjust.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%priceadjust)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'profit.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%profit)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'dividend.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%dividend)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'divrate.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%divrate)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'lumptransfer.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%lumptransfer)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'equity.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%equity)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'caputil.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%caputil)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'deprec.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%deprec)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'tfpadj.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%tfpadj)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'illassetdrop.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%equmSTICKY%illassetdrop)

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ea.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ea)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Eb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Eb)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ec)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Elabor.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Elabor)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ed.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ed)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ewage.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ewage)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Enetlabinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Enetlabinc)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Egrosslabinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Egrosslabinc)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Einc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Einc)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ehours.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ehours)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Enw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Enw)

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACa0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACa0)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACa0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACa0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACb0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACb0)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACb0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACb0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACb0a0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACb0a0)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACb0aP.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACb0aP)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACbN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACbN)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACnw0.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACnw0)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'FRACnw0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%FRACnw0close)

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'EbN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%EbN)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'EbP.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%EbP)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Eadjcost.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Eadjcost)

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'GINIa.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%GINIa)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'GINIb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%GINIb)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'GINInw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%GINInw)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'GINIc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%GINIc)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'GINIinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%GINIinc)

	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_bN.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ec_bN)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_b0close.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ec_b0close)
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_b0far.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,1,irfsave%statsSTICKY%Ec_b0far)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsSTICKY%PERCa(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'PERCa.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsSTICKY%PERCb(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'PERCb.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsSTICKY%PERCnw(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'PERCnw.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsSTICKY%PERCc(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'PERCc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO ip = 1,11;lmat(:,ip) = irfsave%statsSTICKY%PERCinc(ip); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'PERCinc.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,11,lmat)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Ea_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ea_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Eb_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Eb_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Ec_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Einc_nwQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Einc_nwQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Ea_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ea_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Eb_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Eb_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Ec_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Ec_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	DO iq = 1,4;lmat4(:,iq) = irfsave%statsSTICKY%Einc_incQ(iq); END DO	
	OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'Einc_incQ.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,Ttransition,4,lmat4)

	IF(SaveTime1PolicyFns==1) THEN
		DO iy = 1,ngpy
			IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
			IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy
			
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'V_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%V(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'dep_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%d(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'con_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%c(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'sav_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%s(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'hour_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%h(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'bdot_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%bdot(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'gjoint_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%gjoint(:,:,iy))
			IF(ComputeDiscountedMPC==1) THEN
				OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'mpc_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%mpc(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'subeff1ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%subeff1ass(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'subeff2ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%subeff2ass(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'wealtheff1ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%wealtheff1ass(:,:,iy))
				OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/"  // 'wealtheff2ass_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%solnSTICKY(1)%wealtheff2ass(:,:,iy))
			END IF
		END DO
	END IF
	
	IF (SaveCumPolicyFnsIRF==1) THEN
		DO iy = 1,ngpy
			IF (iy<10) WRITE(UNIT=lstring, FMT='(I1)') iy
			IF (iy >= 10) WRITE(UNIT=lstring, FMT='(I2)') iy
			
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'ccum1_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%ccum1(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'ccum2_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%ccum2(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'ccum4_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%ccum4(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'dcum1_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%dcum1(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'dcum2_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%dcum2(:,:,iy))
			OPEN(3, FILE = trim(OutputDirIRF) // "STICKY/" // 'dcum4_T1_y' // trim(lstring) //'.txt', STATUS = 'replace'); CALL WriteMatrixExpon(3,ngpa,ngpb,irfsave%cumSTICKY%dcum4(:,:,iy))
		END DO
		
	END IF	
	

END IF


END SUBROUTINE