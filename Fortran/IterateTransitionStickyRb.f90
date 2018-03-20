SUBROUTINE IterateTransitionStickyRb

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: it,ii,itfg
REAL(8) 	:: lldK,ldiffB,ldiffK,lminmargcost,lpvgovbc,lpvlumpincr,linitlumpincr
REAL(8), DIMENSION(Ttransition) :: lcapital,lcapital1,lbond,lfirmdiscount,lrb,lrb1,lfundbond,lworldbond,lrgov

iteratingtransition = .true.

lminmargcost = 0.01

IF(Display>=1 .and. stickytransition==.true.) write(*,*)' Solving for sticky price transition'

!forward guidance time
itfg = MINLOC(cumdeltatrans, 1, MASK = cumdeltatrans>=ForwardGuideShockQtrs)


!guess capital demand and liquid return

!construct sequence of guesses of capital: assume log linear in capital (constant if a temporary transition)
equmTRANS(1)%capital = equmINITSS%capital
lldK = (log(equmFINALSS%capital)  - log(equmTRANS(1)%capital)) / real(Ttransition)
DO it = 2,Ttransition
	equmTRANS(it)%capital = equmTRANS(1)%capital  *exp(lldK*it)
END DO

!construct sequence of guesses of Rb
equmTRANS(:)%pi = equmINITSS%pi

IF(forwardguide==.false.) THEN
	equmTRANS(:)%rnom = equmINITSS%rnom +phitaylor*equmTRANS(:)%pi + equmTRANS(:)%mpshock
			
ELSE IF(forwardguide==.true. ) THEN
	equmTRANS(1:itfg-1)%rnom = equmINITSS%rnom +phifg*equmTRANS(1:itfg-1)%pi + equmTRANS(1:itfg-1)%mpshock
	equmTRANS(itfg:Ttransition)%rnom = equmINITSS%rnom +phitaylor*equmTRANS(itfg:Ttransition)%pi + equmTRANS(itfg:Ttransition)%mpshock
END IF	
equmTRANS(:)%rb = equmTRANS(:)%rnom - equmTRANS(:)%pi
equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge

!world bond
equmTRANS(1)%worldbond = equmINITSS%worldbond
DO it = 1,Ttransition-1
	CALL WorldBondFunction2( equmTRANS(it)%rb,equmTRANS(it+1)%worldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
	equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + bondadjust*deltatransvec(it)*(equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)
END DO

!fund bond
equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev


!solve phillips curve backwards for marginal costs
IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb
IF (FirmDiscountRate==3) lfirmdiscount = equmINITSS%ra
IF (FirmDiscountRate==4) lfirmdiscount = equmTRANS(:)%rb
IF (FirmDiscountRate==5) lfirmdiscount = equmINITSS%ra

!marginal costs
!final period of transition
it = Ttransition
equmTRANS(it)%mc = (lfirmdiscount(it) 	- (equmFINALSS%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
										- alpha*(equmFINALSS%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
										- alpha*(equmFINALSS%caputil-equmTRANS(it)%caputil)/(equmTRANS(it)%caputil*deltatransvec(it)) &
										- (1.0-alpha)*(equmFINALSS%labor-equmTRANS(it)%labor)/(equmTRANS(it)%labor*deltatransvec(it)) ) *equmFINALSS%pi * theta/ equmTRANS(it)%elast &
										+ (equmTRANS(it)%elast-1.0)/equmTRANS(it)%elast - ((equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it)) * theta/ equmTRANS(it)%elast
equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)

!solve backwards
DO it = Ttransition-1,1,-1
	equmTRANS(it)%mc = (lfirmdiscount(it) 	- (equmTRANS(it+1)%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
											- alpha*(equmTRANS(it+1)%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
											- alpha*(equmTRANS(it+1)%caputil-equmTRANS(it)%caputil)/(equmTRANS(it)%caputil*deltatransvec(it)) &
											- (1.0-alpha)*(equmTRANS(it+1)%labor-equmTRANS(it)%labor)/(equmTRANS(it)%labor*deltatransvec(it)) ) *equmTRANS(it+1)%pi * theta/ equmTRANS(it)%elast &
											+ (equmTRANS(it)%elast-1.0)/equmTRANS(it)%elast - ((equmTRANS(it+1)%pi-equmTRANS(it)%pi)/deltatransvec(it)) * theta/ equmTRANS(it)%elast
	equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)

END DO

equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
equmTRANS(:)%tfpadj = (equmTRANS(:)%tfp**((1.0+utilelast)/utilelastalpha)) * ((equmTRANS(:)%mc*alpha/equmINITSS%rcapital)**(alpha*utilelast/utilelastalpha))
equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
equmTRANS(:)%wage = equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfpadj * (equmTRANS(:)%KNratio**(alpha/utilelastalpha))
equmTRANS(:)%netwage = (1.0-equmTRANS(:)%labtax)*equmTRANS(:)%wage
equmTRANS(:)%caputil = ((equmTRANS(:)%mc*alpha*equmTRANS(:)%tfp/equmINITSS%rcapital) * equmTRANS(:)%KNratio**(alpha-1.0)) ** (utilelast/utilelastalpha)
equmTRANS(:)%output = equmTRANS(:)%tfpadj * (equmTRANS(:)%capital**(alpha/utilelastalpha)) * (equmTRANS(:)%labor**((1.0-alpha)*(1.0+utilelast)/utilelastalpha))
equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / (equmTRANS(:)%tfp* (equmTRANS(:)%caputil**alpha))
equmTRANS(:)%rcapital = ((equmINITSS%rcapital**utilelast) * equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio ) ** (1.0/(1.0+utilelast))
equmTRANS(:)%priceadjust = (theta/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%capital/equmTRANS(:)%KYratio
equmTRANS(:)%profit = (1.0-equmTRANS(:)%mc)*equmTRANS(:)%capital/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust

equmTRANS(:)%deprec = equmINITSS%deprec + (utilelast*equmINITSS%rcapital/(1.0+ utilelast)) * ((equmTRANS(:)%rcapital/equmINITSS%rcapital)**(1.0+utilelast) -1.0)

!solve backward for investment
it = Ttransition
equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
DO it = Ttransition-1,1,-1
	equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
END DO

!dividends and illiquid return
equmTRANS(:)%dividend 	= equmTRANS(:)%profit*(1.0-corptax)
IF(DistributeProfitsInProportion==1) equmTRANS(:)%dividend  = profdistfrac*equmTRANS(:)%dividend 

IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital


equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)

!value of equity component of investmemt fund
IF(DividendFundLumpSum==0) THEN
	equmTRANS(:)%equity = 0.0
	equmTRANS(:)%illassetdrop = 1.0	
ELSE IF(DividendFundLumpSum==1) THEN
	it = Ttransition
	equmTRANS(it)%equity = (equmFINALSS%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%equity = (equmTRANS(it+1)%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
	END DO
equmTRANS(:)%illassetdrop = ((1.0-equmTRANS(1)%fundlev)*equmTRANS(1)%capital + equmTRANS(1)%equity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
END IF

	
!government budget constraint,expenditures and tax rates
IF (AdjGovBudgetConstraint==1) THEN !adjust spending
	equmTRANS(:)%govbond = equmINITSS%govbond
	equmTRANS(:)%labtax = equmINITSS%labtax
	equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
		
	equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
	
	equmTRANS(:)%govexp = equmTRANS(:)%taxrev + equmTRANS(:)%rb*equmINITSS%govbond

ELSE IF(AdjGovBudgetConstraint==2) THEN  !adjust lump sum taxes
	equmTRANS(:)%govbond = equmINITSS%govbond
	equmTRANS(:)%govexp = equmINITSS%govexp
	equmTRANS(:)%labtax = equmINITSS%labtax
	equmTRANS(:)%taxrev = equmTRANS(:)%govexp - equmTRANS(:)%rb*equmINITSS%govbond
	equmTRANS(:)%lumptransfer = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor + corptax*equmTRANS(:)%profit + equmTRANS(:)%rb*equmINITSS%govbond - equmTRANS(:)%govexp
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)

			
ELSE IF(AdjGovBudgetConstraint==3) THEN !adjust debt
	IF(GovExpConstantFracOutput==0) equmTRANS(:)%govexp = equmINITSS%govexp
	IF(GovExpConstantFracOutput==1) equmTRANS(:)%govexp = equmTRANS(:)%output*equmINITSS%govexp/equmINITSS%output

	equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
	equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)

	!compute required increase in lumptransfer
	lrgov = equmTRANS(:)%rb
	lpvgovbc = equmFINALSS%govbond
	lpvlumpincr = 0.0
	DO it = Ttransition,1,-1
		lpvgovbc = (lpvgovbc + deltatransvec(it)*(equmTRANS(it)%govexp - equmTRANS(it)%taxrev))/(1.0+deltatransvec(it)*lrgov(it))
		IF(cumdeltatrans(it)>=taxincrstart) lpvlumpincr = (lpvlumpincr + deltatransvec(it))/(1.0+deltatransvec(it)*(lrgov(it)+taxincrdecay))
		IF(cumdeltatrans(it)<taxincrstart) lpvlumpincr = lpvlumpincr/(1.0+deltatransvec(it)*lrgov(it))
	END DO	

	linitlumpincr = (equmINITSS%govbond-lpvgovbc) / lpvlumpincr
	DO it = 1,Ttransition
		IF(cumdeltatrans(it)>=taxincrstart) equmTRANS(it)%lumptransfer = equmTRANS(it)%lumptransfer + linitlumpincr*exp(-taxincrdecay*(cumdeltatrans(it)-taxincrstart))
	END DO

	equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
		
	equmTRANS(Ttransition)%govbond = equmFINALSS%govbond
	DO it = Ttransition-1,2,-1
		equmTRANS(it)%govbond = (equmTRANS(it+1)%govbond - deltatransvec(it)*(equmTRANS(it)%taxrev-equmTRANS(it)%govexp)) / (1.0+deltatransvec(it)*lrgov(it))
	END DO
	equmTRANS(1)%govbond = equmINITSS%govbond

	equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + (equmTRANS(:)%rb-lrgov(:))*equmTRANS(:)%govbond
	equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
	
ELSE IF(AdjGovBudgetConstraint==4) THEN  !adjust proportional tax rate
	equmTRANS(:)%govbond = equmINITSS%govbond
	equmTRANS(:)%govexp = equmINITSS%govexp
	equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer
	equmTRANS(:)%taxrev = equmTRANS(:)%govexp - equmTRANS(:)%rb*equmINITSS%govbond

	IF(DistributeProfitsInProportion == 0 .or. TaxHHProfitIncome == 0) equmTRANS(:)%labtax  = (equmTRANS(:)%lumptransfer - corptax*equmTRANS(:)%profit - equmTRANS(:)%rb*equmINITSS%govbond + equmTRANS(:)%govexp) / (equmTRANS(:)%wage*equmTRANS(:)%labor)
	IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%labtax  = (equmTRANS(:)%lumptransfer - corptax*equmTRANS(:)%profit - equmTRANS(:)%rb*equmINITSS%govbond + equmTRANS(:)%govexp) / (equmTRANS(:)%wage*equmTRANS(:)%labor + (1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax))

END IF

equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond

ii = 1	 
ldiffK = 1.0
ldiffB = 1.0
DO WHILE (ii<=maxitertranssticky .and. max(ldiffK,ldiffB)>toltransition )
	!solve for transtion
	CALL Transition
	
	!computed implied equilibrium quantities
	lbond = statsTRANS(:)%Eb
	lcapital = (statsTRANS(:)%Ea - equmTRANS(:)%equity)/ (1.0 - equmTRANS(:)%fundlev)
	IF(ConvergenceRelToOutput==0) THEN
		ldiffK= maxval(abs(lcapital/equmTRANS(:)%capital - 1.0))
		ldiffB= maxval(abs(lbond/equmTRANS(:)%bond - 1.0))
	ELSEIF(ConvergenceRelToOutput==1) THEN
		ldiffK= maxval(abs(lcapital-equmTRANS(:)%capital)/equmINITSS%output)
		ldiffB= maxval(abs(lbond-equmTRANS(:)%bond)/equmINITSS%output)
	END IF
	IF (Display>=1) write(*,"(A,I,A)") '  Transition iter ',ii, ':'
	IF (Display>=1) write(*,"(A,E10.3,A,E10.3,A,E10.3)") '   K err',ldiffK, ',  B err',ldiffB
	IF (Display>=1) write(*,*) '   household bond',lbond(2), ',  target bond',equmTRANS(2)%bond
	
	!update capital and interest rate
	IF (ii<maxitertranssticky .and. max(ldiffK,ldiffB)>toltransition ) THEN
		CALL PartialUpdate(Ttransition-1,stepstickytransK,equmTRANS(2:Ttransition)%capital,lcapital(2:Ttransition),lcapital1(2:Ttransition))
		lcapital1(1)  = equmTRANS(1)%capital
		equmTRANS(:)%capital = lcapital1
		
		lfundbond = -lcapital*equmTRANS(:)%fundlev
		lworldbond = -lbond - equmTRANS(:)%govbond - lfundbond

		it = Ttransition
		CALL WorldBondInverse2( (equmFINALSS%worldbond-lworldbond(it))/(bondadjust*deltatransvec(it)) + lworldbond(it) &
							,lrb(it),equmINITSS%worldbond,equmINITSS%rb,bondelast)
		DO it = Ttransition-1,1,-1
			CALL WorldBondInverse2( (lworldbond(it+1)-lworldbond(it))/(bondadjust*deltatransvec(it)) + lworldbond(it) &
								,lrb(it),equmINITSS%worldbond,equmINITSS%rb,bondelast)
		END DO

		CALL PartialUpdate(Ttransition,stepstickytransB,equmTRANS(:)%rb,lrb,lrb1)
		equmTRANS(:)%rb = lrb1
			

	ElSE
		!run distribution stats with full
		iteratingtransition = .false.
		CALL Transition
		equmTRANS(:)%capital = lcapital
		equmTRANS(:)%bond = lbond
		equmTRANS(:)%rb = lrb
		
	END IF
	
	equmTRANS(:)%rborr = equmTRANS(:)%rb + equmTRANS(:)%borrwedge
	
	!world bond
	equmTRANS(1)%worldbond = equmINITSS%worldbond
	DO it = 1,Ttransition-1
		CALL WorldBondFunction2( equmTRANS(it)%rb,equmTRANS(it+1)%worldbond,equmINITSS%worldbond,equmINITSS%rb,bondelast)
		equmTRANS(it+1)%worldbond = equmTRANS(it)%worldbond + bondadjust*deltatransvec(it)*(equmTRANS(it+1)%worldbond-equmTRANS(it)%worldbond)
	END DO

	!fund bond
	equmTRANS(:)%fundbond = -equmTRANS(:)%capital*equmTRANS(:)%fundlev

	!inflation and nominal interest rates
	IF(forwardguide==.false.) THEN
		equmTRANS(:)%pi = (equmTRANS(:)%rb - equmINITSS%rnom - equmTRANS(:)%mpshock) / (phitaylor-1.0) !taylor rule
		equmTRANS(:)%rnom = equmTRANS(:)%rb + equmTRANS(:)%pi !fisher equn
		
	ELSE IF(forwardguide==.true.) THEN

		equmTRANS(1:itfg-1)%pi = (equmTRANS(1:itfg-1)%rb - equmINITSS%rnom - equmTRANS(1:itfg-1)%mpshock) / (phifg-1.0) !taylor rule
		equmTRANS(1:itfg-1)%rnom = equmTRANS(1:itfg-1)%rb + equmTRANS(1:itfg-1)%pi !fisher equn
		
		equmTRANS(itfg:Ttransition)%pi = (equmTRANS(itfg:Ttransition)%rb - equmINITSS%rnom - equmTRANS(itfg:Ttransition)%mpshock) / (phitaylor-1.0) !taylor rule
		equmTRANS(itfg:Ttransition)%rnom = equmTRANS(itfg:Ttransition)%rb + equmTRANS(itfg:Ttransition)%pi !fisher equn
			
	END IF		



	!solve phillips curve backwards for marginal costs
	IF (FirmDiscountRate==1) lfirmdiscount = equmTRANS(:)%rho
	IF (FirmDiscountRate==2) lfirmdiscount = equmINITSS%rb
	IF (FirmDiscountRate==3) lfirmdiscount = equmINITSS%ra
	IF (FirmDiscountRate==4) lfirmdiscount = equmTRANS(:)%rb
	IF (FirmDiscountRate==5) lfirmdiscount = equmTRANS(:)%ra

	!labor
	equmTRANS(:)%labor = statsTRANS(:)%Elabor

	!marginal costs
	!final period of transition
	it = Ttransition
	equmTRANS(it)%mc = (lfirmdiscount(it) 	- (equmFINALSS%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
											- alpha*(equmFINALSS%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
											- alpha*(equmFINALSS%caputil-equmTRANS(it)%caputil)/(equmTRANS(it)%caputil*deltatransvec(it)) &
											- (1.0-alpha)*(equmFINALSS%labor-equmTRANS(it)%labor)/(equmTRANS(it)%labor*deltatransvec(it)) ) *equmFINALSS%pi * theta/ equmTRANS(it)%elast &
											+ (equmTRANS(it)%elast-1.0)/equmTRANS(it)%elast - ((equmFINALSS%pi-equmTRANS(it)%pi)/deltatransvec(it)) * theta/ equmTRANS(it)%elast
	equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)

	!solve backwards
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%mc = (lfirmdiscount(it) 	- (equmTRANS(it+1)%tfp-equmTRANS(it)%tfp)/(equmTRANS(it)%tfp*deltatransvec(it)) &
												- alpha*(equmTRANS(it+1)%capital-equmTRANS(it)%capital)/(equmTRANS(it)%capital*deltatransvec(it)) &
												- alpha*(equmTRANS(it+1)%caputil-equmTRANS(it)%caputil)/(equmTRANS(it)%caputil*deltatransvec(it)) &
												- (1.0-alpha)*(equmTRANS(it+1)%labor-equmTRANS(it)%labor)/(equmTRANS(it)%labor*deltatransvec(it)) ) *equmTRANS(it+1)%pi * theta/ equmTRANS(it)%elast &
												+ (equmTRANS(it)%elast-1.0)/equmTRANS(it)%elast - ((equmTRANS(it+1)%pi-equmTRANS(it)%pi)/deltatransvec(it)) * theta/ equmTRANS(it)%elast
		equmTRANS(it)%mc = max(lminmargcost,equmTRANS(it)%mc)

	END DO
 

	equmTRANS(:)%gap = equmTRANS(:)%elast*equmTRANS(:)%mc / (equmTRANS(:)%elast-1.0) - 1.0
	equmTRANS(:)%tfpadj = (equmTRANS(:)%tfp**((1.0+utilelast)/utilelastalpha)) * ((equmTRANS(:)%mc*alpha/equmINITSS%rcapital)**(alpha*utilelast/utilelastalpha))
	equmTRANS(:)%KNratio = equmTRANS(:)%capital/equmTRANS(:)%labor
	equmTRANS(:)%wage = equmTRANS(:)%mc*(1.0-alpha)* equmTRANS(:)%tfpadj * (equmTRANS(:)%KNratio**(alpha/utilelastalpha))
	equmTRANS(:)%netwage = (1.0-equmTRANS(:)%labtax)*equmTRANS(:)%wage
	equmTRANS(:)%caputil = ((equmTRANS(:)%mc*alpha*equmTRANS(:)%tfp/equmINITSS%rcapital) * equmTRANS(:)%KNratio**(alpha-1.0)) ** (utilelast/utilelastalpha)
	equmTRANS(:)%output = equmTRANS(:)%tfpadj * (equmTRANS(:)%capital**(alpha/utilelastalpha)) * (equmTRANS(:)%labor**((1.0-alpha)*(1.0+utilelast)/utilelastalpha))
	equmTRANS(:)%KYratio = (equmTRANS(:)%KNratio**(1.0-alpha)) / (equmTRANS(:)%tfp* (equmTRANS(:)%caputil**alpha))
	equmTRANS(:)%rcapital = ((equmINITSS%rcapital**utilelast) * equmTRANS(:)%mc * alpha / equmTRANS(:)%KYratio ) ** (1.0/(1.0+utilelast))
	equmTRANS(:)%priceadjust = (theta/2.0)*(equmTRANS(:)%pi**2)*equmTRANS(:)%capital/equmTRANS(:)%KYratio
	equmTRANS(:)%profit = (1.0-equmTRANS(:)%mc)*equmTRANS(:)%capital/equmTRANS(:)%KYratio - equmTRANS(:)%priceadjust

	equmTRANS(:)%deprec = equmINITSS%deprec + (utilelast*equmINITSS%rcapital/(1.0+ utilelast)) * ((equmTRANS(:)%rcapital/equmINITSS%rcapital)**(1.0+utilelast) -1.0)

	!solve backward for investment
	it = Ttransition
	equmTRANS(it)%investment = (equmFINALSS%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
	DO it = Ttransition-1,1,-1
		equmTRANS(it)%investment = (equmTRANS(it+1)%capital-equmTRANS(it)%capital)/deltatransvec(it) + equmTRANS(it)%deprec*equmTRANS(it)%capital
	END DO

	!dividends and illiquid return
	equmTRANS(:)%dividend 	= equmTRANS(:)%profit*(1.0-corptax)
	IF(DistributeProfitsInProportion==1) equmTRANS(:)%dividend  = profdistfrac*equmTRANS(:)%dividend

	IF(DividendFundLumpSum==1) equmTRANS(:)%divrate = 0.0
	IF(DividendFundLumpSum==0) equmTRANS(:)%divrate = equmTRANS(:)%dividend/equmTRANS(:)%capital

	equmTRANS(:)%ra = (equmTRANS(:)%rcapital*equmTRANS(:)%caputil - equmTRANS(:)%deprec + equmTRANS(:)%divrate - equmTRANS(:)%fundlev*equmTRANS(:)%rb) / (1.0-equmTRANS(:)%fundlev)

	!value of equity component of investmemt fund
	IF(DividendFundLumpSum==0) THEN
		equmTRANS(:)%equity = 0.0
		equmTRANS(:)%illassetdrop = 1.0	
	ELSE IF(DividendFundLumpSum==1) THEN
		it = Ttransition
		equmTRANS(it)%equity = (equmFINALSS%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
		DO it = Ttransition-1,1,-1
			equmTRANS(it)%equity = (equmTRANS(it+1)%equity + equmTRANS(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
		END DO
		equmTRANS(:)%illassetdrop = ((1.0-equmTRANS(1)%fundlev)*equmTRANS(1)%capital + equmTRANS(1)%equity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
	END IF

	!government budget constraint,expenditures and tax rates
	IF (AdjGovBudgetConstraint==1) THEN !adjust spending
		equmTRANS(:)%govbond = equmINITSS%govbond
		equmTRANS(:)%labtax = equmINITSS%labtax
		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
			
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
	
		equmTRANS(:)%govexp = equmTRANS(:)%taxrev + equmTRANS(:)%rb*equmINITSS%govbond

	ELSE IF(AdjGovBudgetConstraint==2) THEN  !adjust lump sum taxes
		equmTRANS(:)%govbond = equmINITSS%govbond
		equmTRANS(:)%govexp = equmINITSS%govexp
		equmTRANS(:)%labtax = equmINITSS%labtax
		equmTRANS(:)%taxrev = equmTRANS(:)%govexp - equmTRANS(:)%rb*equmINITSS%govbond
		equmTRANS(:)%lumptransfer = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor + corptax*equmTRANS(:)%profit + equmTRANS(:)%rb*equmINITSS%govbond - equmTRANS(:)%govexp
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
					
	ELSE IF(AdjGovBudgetConstraint==3) THEN !adjust debt
		IF(GovExpConstantFracOutput==0) equmTRANS(:)%govexp = equmINITSS%govexp
		IF(GovExpConstantFracOutput==1) equmTRANS(:)%govexp = equmTRANS(:)%output*equmINITSS%govexp/equmINITSS%output

		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer		
		
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)

		!compute required increase in lumptransfer
		lrgov = equmTRANS(:)%rb
		lpvgovbc = equmFINALSS%govbond
		lpvlumpincr = 0.0
		DO it = Ttransition,1,-1
			lpvgovbc = (lpvgovbc + deltatransvec(it)*(equmTRANS(it)%govexp - equmTRANS(it)%taxrev))/(1.0+deltatransvec(it)*lrgov(it))
			IF(cumdeltatrans(it)>=taxincrstart) lpvlumpincr = (lpvlumpincr + deltatransvec(it))/(1.0+deltatransvec(it)*(lrgov(it)+taxincrdecay))
			IF(cumdeltatrans(it)<taxincrstart) lpvlumpincr = lpvlumpincr/(1.0+deltatransvec(it)*lrgov(it))
		END DO	

		linitlumpincr = (equmINITSS%govbond-lpvgovbc) / lpvlumpincr
		DO it = 1,Ttransition
			IF(cumdeltatrans(it)>=taxincrstart) equmTRANS(it)%lumptransfer = equmTRANS(it)%lumptransfer + linitlumpincr*exp(-taxincrdecay*(cumdeltatrans(it)-taxincrstart))
		END DO

		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)
		
		equmTRANS(Ttransition)%govbond = equmFINALSS%govbond
		DO it = Ttransition-1,2,-1
			equmTRANS(it)%govbond = (equmTRANS(it+1)%govbond - deltatransvec(it)*(equmTRANS(it)%taxrev-equmTRANS(it)%govexp)) / (1.0+deltatransvec(it)*lrgov(it))
		END DO
		equmTRANS(1)%govbond = equmINITSS%govbond

		equmTRANS(:)%lumptransfer = equmTRANS(:)%lumptransfer + (equmTRANS(:)%rb-lrgov(:))*equmTRANS(:)%govbond
		equmTRANS(:)%taxrev = equmTRANS(:)%labtax*equmTRANS(:)%wage*equmTRANS(:)%labor - equmTRANS(:)%lumptransfer + corptax*equmTRANS(:)%profit
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%taxrev = equmTRANS(:)%taxrev + equmTRANS(:)%labtax*(1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax)

	ELSE IF(AdjGovBudgetConstraint==4) THEN  !adjust proportional tax rate
		equmTRANS(:)%govbond = equmINITSS%govbond
		equmTRANS(:)%govexp = equmINITSS%govexp
		equmTRANS(:)%lumptransfer = equmINITSS%lumptransfer
		equmTRANS(:)%taxrev = equmTRANS(:)%govexp - equmTRANS(:)%rb*equmINITSS%govbond

		IF(DistributeProfitsInProportion == 0 .or. TaxHHProfitIncome == 0) equmTRANS(:)%labtax  = (equmTRANS(:)%lumptransfer - corptax*equmTRANS(:)%profit - equmTRANS(:)%rb*equmINITSS%govbond + equmTRANS(:)%govexp) / (equmTRANS(:)%wage*equmTRANS(:)%labor)
		IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) equmTRANS(:)%labtax  = (equmTRANS(:)%lumptransfer - corptax*equmTRANS(:)%profit - equmTRANS(:)%rb*equmINITSS%govbond + equmTRANS(:)%govexp) / (equmTRANS(:)%wage*equmTRANS(:)%labor + (1.0-profdistfrac)*equmTRANS(:)%profit*(1.0-corptax))
	

	END IF


	!household bonds
	equmTRANS(:)%bond = -equmTRANS(:)%worldbond - equmTRANS(:)%govbond - equmTRANS(:)%fundbond
	
	ii = ii+1	
END DO

IF(stickytransition==.true.) THEN
	irfpointer%equmSTICKY = equmTRANS
	irfpointer%statsSTICKY = statsTRANS
	irfpointer%solnSTICKY = solnTRANS
END IF


END SUBROUTINE IterateTransitionStickyRb