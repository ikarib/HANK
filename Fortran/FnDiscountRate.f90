REAL(8) FUNCTION FnDiscountRate(lrhoT)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lrhoT
REAL(8)				:: lKNratio,lKYratio

rho = -log(logistic(lrhoT))

IF (Display>=2) THEN
	WRITE(*,*) '*******************************************'
	WRITE(*,*) ' ITERATION : 			',neqmiter
	WRITE(*,*) ' rho guess: 			',rho
	WRITE(*,*) '  target KY ratio: ',KYratio
	WRITE(*,*) '  target KN ratio: ',KNratio
END IF

CALL Grids
CALL IterateBellman
CALL StationaryDistribution
CALL DistributionStatistics

IF(DividendFundLumpSum==0) capital = Ea/(1.0-fundlev)
IF(DividendFundLumpSum==1) THEN
	IF(DistributeProfitsInProportion==0) capital = Ea/(1.0-fundlev + (1.0-mc)*(1.0-corptax)/(ra*KYratio))
	IF(DistributeProfitsInProportion==1) capital = Ea/(1.0-fundlev + (1.0-mc)*(1.0-corptax)*profdistfrac/(ra*KYratio))
END IF

labor = Elabor
equity = Ea - (1.0-fundlev)*capital

lKNratio = capital / labor
lKYratio = (lKNratio**(1.0-alpha)) / tfp

IF(Display>=1) write(*,"(A,I2,A,E11.4,A,E11.4)") '  Rho iter ',neqmiter, ', rho ',rho, ', K/N err',lKNratio/KNratio - 1.0

FnDiscountRate = lKNratio/KNratio - 1.0


IF(CalibrateDiscountRate==1) OPEN(3, FILE = trim(OutputDir) // 'DiscountRateCalibration.txt', ACCESS = 'append')
IF(CalibrateRhoAtInitialGuess==1) OPEN(3, FILE = trim(OutputDir) // 'DiscountRateAtInitialGuess.txt', ACCESS = 'append')
WRITE(3,*) '*******************************************'
WRITE(3,*) ' ITERATION : 			',neqmiter
WRITE(3,*) ' rho guess: 			',rho
WRITE(3,*) '  target KY ratio: ',KYratio
WRITE(3,*) '  implied KY ratio: ',lKYratio
WRITE(3,*) '  target KN ratio: ',KNratio
WRITE(3,*) '  implied KN ratio: ',lKNratio
WRITE(3,*) '  relative error: ',FnDiscountRate
CLOSE(3)

neqmiter = neqmiter+1

!implied aggregate statistics
bond = Eb
investment = deprec*capital
priceadjust = 0.0
profit = (1.0-mc)*capital/KYratio - priceadjust
dividend = profit*(1.0-corptax)
IF(DistributeProfitsInProportion==1) dividend = profdistfrac*dividend 
output = tfp*(capital**alpha)*(labor**(1.0-alpha))
fundbond = -capital*fundlev
bondelast = bondelastrelgdp*output
caputil 	= 1.0
tfpadj = ((tfp**(1.0+utilelast)) * (mc*alpha/rcapital)**(alpha*utilelast))**(1.0/utilelastalpha)
taxrev = labtax*wage*labor - lumptransfer + corptax*profit
IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) taxrev = taxrev + labtax*(1.0-profdistfrac)*profit*(1.0-corptax)

IF(GovBondResidualZeroWorld==0) THEN
	govbond = -ssdebttogdp*output
	govexp = taxrev + rb*govbond 
	worldbond = -bond-govbond-fundbond
ELSE IF(GovBondResidualZeroWorld==1) THEN
	worldbond = 0.0
	govbond = -bond-worldbond-fundbond
	govexp = taxrev + rb*govbond		
END IF

END FUNCTION FnDiscountRate
