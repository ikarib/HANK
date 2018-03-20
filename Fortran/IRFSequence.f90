SUBROUTINE IRFSequence(lIRFDir)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

CHARACTER, INTENT(IN) :: lIRFDir*20
CHARACTER	:: lstring*80
INTEGER	:: it,ipe
REAL(8) :: lequity

irfpointer => irfstruct(0)		

IF(SolveStickyPriceTransition==1) THEN
	stickytransition = .true.
	IF(OneAssetNoCapital==0) CALL IterateTransitionStickyRb
	IF(OneAssetNoCapital==1) CALL IterateTransOneAssetStickyRb
	IF (SaveCumPolicyFnsIRF==1) CALL CumulativeConsTransition
	IF(ComputeDiscountedMPC==1) CALL DiscountedMPCTransition		
	stickytransition = .false.		
END IF	

irfsave => irfstruct(0)
OutputDirIRF = trim(OutputDir) // "IRF_" // trim(lIRFDir) // "/NOFS/"
CALL SaveIRFOutput

IF (DoPriceExperiments==1) THEN
	irfpointer => irfpriceexp
	
	DO ipe = 1,15
				
		IF(SolveStickyPriceTransition==1) THEN
			
			write(*,*) "Solving for sticky price transition without ZLB, price experiment ", ipe
		
			DO it = 1,Ttransition
				equmTRANS(it) = equmINITSS
			END DO

			SELECT CASE(ipe)
				CASE(1) !change wage only
					equmTRANS%wage = irfstruct(0)%equmSTICKY%wage
					equmTRANS%netwage = irfstruct(0)%equmSTICKY%netwage
					equmTRANS%labtax = irfstruct(0)%equmSTICKY%labtax					

				CASE(2) !only change profits
					equmTRANS%profit = irfstruct(0)%equmSTICKY%profit

				CASE(3) !only change profits and wage
					equmTRANS%profit = irfstruct(0)%equmSTICKY%profit
					equmTRANS%wage = irfstruct(0)%equmSTICKY%wage
					equmTRANS%netwage = irfstruct(0)%equmSTICKY%netwage
					equmTRANS%labtax = irfstruct(0)%equmSTICKY%labtax					
				
				CASE(4) !only change rb, rborr only
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr

				CASE(5) !only change ra only
					equmTRANS%ra = irfstruct(0)%equmSTICKY%ra

				CASE(6) !only change illiquid asset drop
					equmTRANS%illassetdrop  = irfstruct(0)%equmSTICKY%illassetdrop

				CASE(7) !only change transfers
					equmTRANS%lumptransfer = irfstruct(0)%equmSTICKY%lumptransfer


				CASE(8) !change all 6
					equmTRANS%wage = irfstruct(0)%equmSTICKY%wage
					equmTRANS%netwage = irfstruct(0)%equmSTICKY%netwage
					equmTRANS%labtax = irfstruct(0)%equmSTICKY%labtax					
					equmTRANS%lumptransfer = irfstruct(0)%equmSTICKY%lumptransfer
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr
					equmTRANS%ra = irfstruct(0)%equmSTICKY%ra
					equmTRANS%illassetdrop  = irfstruct(0)%equmSTICKY%illassetdrop
					equmTRANS%profit = irfstruct(0)%equmSTICKY%profit

 				CASE(9) !change lump transfer by direct effect from government interest payments
					equmTRANS%lumptransfer = equmINITSS%lumptransfer + (irfstruct(0)%equmSTICKY%rb - equmTRANS%rb)*equmINITSS%govbond

 				CASE(10) !change lump transfer, not including direct effect from govt interest payments
					equmTRANS%lumptransfer = irfstruct(0)%equmSTICKY%lumptransfer - (irfstruct(0)%equmSTICKY%rb - equmTRANS%rb)*equmINITSS%govbond

				CASE(11) !change rb, rborr, and change ra by same amount as rb
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr
					equmTRANS%ra = equmINITSS%ra + (irfstruct(0)%equmSTICKY%rb - equmINITSS%rb)

				CASE(12) !change ra, and change rb, rborr by same amount as ra
					equmTRANS%ra = irfstruct(0)%equmSTICKY%ra
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb + (irfstruct(0)%equmSTICKY%ra - equmINITSS%ra)
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr + (irfstruct(0)%equmSTICKY%ra - equmINITSS%ra)

				CASE(13) !change only proportional labor tax
					equmTRANS%labtax = irfstruct(0)%equmSTICKY%labtax					
					equmTRANS%netwage = equmINITSS%wage*(1.0-irfstruct(0)%equmSTICKY%labtax)

				CASE(14) !change rb, rborr, and change ra by same amount as rb, and discount eqm profits at implied ra
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr
					equmTRANS%ra = equmINITSS%ra + (irfstruct(0)%equmSTICKY%rb - equmINITSS%rb)
					IF(DividendFundLumpSum==0) THEN
						equmTRANS(:)%illassetdrop = 1.0	
					ELSE IF(DividendFundLumpSum==1) THEN
						it = Ttransition
						lequity = (equmFINALSS%equity + irfstruct(0)%equmSTICKY(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
						DO it = Ttransition-1,1,-1
							lequity = (lequity + irfstruct(0)%equmSTICKY(it)%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
						END DO
						equmTRANS(:)%illassetdrop = ((1.0-irfstruct(0)%equmSTICKY(1)%fundlev) * irfstruct(0)%equmSTICKY(1)%capital + lequity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
					END IF					
					
				CASE(15) !change rb, rborr, and change ra by same amount as rb, and discount initimake al profits at implied ra
					equmTRANS%rb = irfstruct(0)%equmSTICKY%rb
					equmTRANS%rborr = irfstruct(0)%equmSTICKY%rborr
					equmTRANS%ra = equmINITSS%ra + (irfstruct(0)%equmSTICKY%rb - equmINITSS%rb)
					IF(DividendFundLumpSum==0) THEN
						equmTRANS(:)%illassetdrop = 1.0	
					ELSE IF(DividendFundLumpSum==1) THEN
						it = Ttransition
						lequity = (equmFINALSS%equity + equmINITSS%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
						DO it = Ttransition-1,1,-1
							lequity = (lequity + equmINITSS%dividend*deltatransvec(it)) / (1.0+deltatransvec(it)*equmTRANS(it)%ra)		
						END DO
						equmTRANS(:)%illassetdrop = ((1.0-irfstruct(0)%equmSTICKY(1)%fundlev) * irfstruct(0)%equmSTICKY(1)%capital + lequity) / ((1.0-equmINITSS%fundlev)*equmINITSS%capital + equmINITSS%equity)
					END IF					
					
			END SELECT	
			CALL Transition

			IF (SaveCumPolicyFnsIRF==1) CALL CumulativeConsTransition

			irfpointer%equmSTICKY = equmTRANS
			irfpointer%statsSTICKY = statsTRANS
			irfpointer%solnSTICKY = solnTRANS
			irfpointer%cumSTICKY = CumulativePolicyType(ccum1,ccum2,ccum4,dcum1,dcum2,dcum4)

			IF(ComputeDiscountedMPC==1) CALL DiscountedMPCTransition	

		END IF

		irfsave => irfpriceexp
		IF(ipe<10) WRITE(UNIT=lstring, FMT='(I1)') ipe
		IF(ipe>=10) WRITE(UNIT=lstring, FMT='(I2)') ipe
		OutputDirIRF = trim(OutputDir) // "IRF_" // trim(lIRFDir) // "/PE"// trim(lstring) // "/"		
		CALL SaveIRFOutput
		
	END DO
	
END IF



END SUBROUTINE IRFSequence