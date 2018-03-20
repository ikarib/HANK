SUBROUTINE SolveSteadyStateEqum

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8) 	:: ldiffKN, lKNratio,lKYratio,lwage,lrcapital,lstepKN,ldiffprof,lprofit

OPEN(3, FILE = trim(OutputDir) // 'SteadyStateEqumIteration.txt', STATUS = 'replace'); CLOSE(3)

converged = .false.
neqmiter = 1
IF(calibrating == .false.) lstepKN = stepequmss
IF(calibrating == .true. ) lstepKN = 0.2*stepequmss
ldiffKN = 1.0

IF(OneAssetNoCapital==0) THEN
	DO WHILE (neqmiter<=maxiterequmss .and. ldiffKN>tolequmss )

		IF (Display>=2) THEN
			WRITE(*,*) '*******************************************'
			WRITE(*,*) ' ITERATION : 			',neqmiter
			WRITE(*,*) ' r guess: 			',rcapital
			WRITE(*,*) '  implied ra: 		',ra
			WRITE(*,*) '  implied borr rate: ',rborr
			WRITE(*,*) '  implied wage: ',wage
			WRITE(*,*) '  implied KY ratio: ',KYratio
			WRITE(*,*) '  implied KN firm: 	',KNratio
		END IF

		IF(initialSS == .true.) CALL Grids
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
		lwage = mc*(1.0-alpha)* tfp * (lKNratio**alpha)
		lrcapital = mc*alpha/lKYratio

		ldiffKN= abs(lKNratio/KNratio - 1.0)
		IF (Display>=1) write(*,"(A,I2,A,E11.4)") ' Steady state equm iter ',neqmiter, ', K/N error',lKNratio/KNratio - 1.0

		OPEN(3, FILE = trim(OutputDir) // 'SteadyStateEqumIteration.txt', ACCESS = 'append')
		WRITE(3,*) '*******************************************'
		WRITE(3,*) ' ITERATION : 			',neqmiter
		WRITE(3,*) ' r guess: 	',rcapital
		WRITE(3,*) '  implied ra: 	',ra
		WRITE(3,*) '  implied borr rate: ',rborr	
		WRITE(3,*) '  implied wage: ',wage
		WRITE(3,*) '  implied KY ratio: ',KYratio
		WRITE(3,*) '  actual KY ratio: ',lKYratio
		WRITE(3,*) '  implied KN firm: 	',KNratio
		WRITE(3,*) '  implied KN hh: 	',lKNratio
		WRITE(3,*) '  relative error: ',lKNratio/KNratio - 1.0
		CLOSE(3)
		
		!update KN ratio
		IF (neqmiter<=maxiterequmss .and. ldiffKN>tolequmss ) THEN
			KNratio = (1.0-lstepKN)*KNratio +lstepKN*lKNratio
		ELSE
			KNratio = lKNratio
		END IF

		KYratio = (KNratio**(1.0-alpha)) / tfp
		profit = (1.0-mc)*capital/KYratio
		rcapital = mc*alpha/KYratio
		wage = mc*(1.0-alpha)*tfp*(KNratio**alpha)
		netwage = (1.0-labtax)*wage
		IF(DividendFundLumpSum==1) divrate = 0.0
		IF(DividendFundLumpSum==0) divrate =  (1.0-corptax)*(1.0-mc)/KYratio !outside of steady state include price adjustments
		IF(DistributeProfitsInProportion==1) divrate =  profdistfrac*divrate
		ra = (rcapital - deprec + divrate - fundlev*rb)/(1.0-fundlev)
	
		neqmiter = neqmiter+1 

	END DO

ELSE IF(OneAssetNoCapital==1) THEN

	IF(DistributeProfitsInProportion==0) THEN
		IF(initialSS == .true.)CALL Grids
		CALL IterateBellman
		CALL StationaryDistribution
		CALL DistributionStatistics
		labor = Elabor
		profit = (1.0-mc)*tfp*labor
	
	ELSE IF(DistributeProfitsInProportion==1) THEN  	!iterate on profits
		ldiffprof=1
		DO WHILE (neqmiter<=maxiterequmss .and. ldiffprof>tolequmss )
			IF(initialSS == .true.)CALL Grids
			CALL IterateBellman
			CALL StationaryDistribution
			CALL DistributionStatistics
			labor = Elabor
			lprofit = (1.0-mc)*tfp*labor
			
			ldiffprof= abs(lprofit/profit - 1.0)			
			IF (Display>=1) write(*,"(A,I2,A,E11.4)") ' Steady state equm iter ',neqmiter, ', Profit error',lprofit/profit - 1.0
			
			profit = lprofit
			neqmiter = neqmiter+1 
			
		END DO	
	END IF
	
	capital = 0.0
	equity = 0.0
	KYratio = 0.0
	KNratio = 0.0
END IF

bond = Eb
investment = deprec*capital
priceadjust = 0.0
dividend = profit*(1.0-corptax)
IF(DistributeProfitsInProportion==1) dividend = profdistfrac*dividend 
output = tfp*(capital**alpha)*(labor**(1.0-alpha))
fundbond = -capital*fundlev
bondelast = bondelastrelgdp*output
caputil 	= 1.0
IF(OneAssetNoCapital==0) tfpadj = ((tfp**(1.0+utilelast)) * (mc*alpha/rcapital)**(alpha*utilelast))**(1.0/utilelastalpha)
IF(OneAssetNoCapital==1) tfpadj = tfp
taxrev = labtax*wage*labor - lumptransfer + corptax*profit
IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) taxrev = taxrev + labtax*(1.0-profdistfrac)*profit*(1.0-corptax)
illassetdrop = 1.0

IF(GovBondResidualZeroWorld==0) THEN
	govbond = -ssdebttogdp*output
	govexp = taxrev + rb*govbond 
	worldbond = -bond-govbond-fundbond
ELSE IF(GovBondResidualZeroWorld==1) THEN
	worldbond = 0.0
	govbond = -bond-worldbond-fundbond
	govexp = taxrev + rb*govbond		
END IF


END SUBROUTINE SolveSteadyStateEqum

