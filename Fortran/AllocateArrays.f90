SUBROUTINE AllocateArrays

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER	:: it,ifs

write(*,*) 'Allocating arrays for transition'

DO it = 1,Ttransition
 	CALL AllocateSolutionType(solnTRANS(it))
END DO

DO it = 1,Ttransition
	ifs = 0
	IF(SolveStickyPriceTransition==1) CALL AllocateSolutionType(irfstruct(ifs)%solnSTICKY(it))

	IF (DoPriceExperiments==1) THEN
		IF(SolveStickyPriceTransition==1) CALL AllocateSolutionType(irfpriceexp%solnSTICKY(it))
	END IF
	
END DO	

ifs = 0
IF (SaveCumPolicyFnsIRF==1)	THEN
	IF(SolveStickyPriceTransition==1) CALL AllocateCumulativePolicyType(irfstruct(ifs)%cumSTICKY)

	IF (DoPriceExperiments==1) THEN
		IF(SolveStickyPriceTransition==1) CALL AllocateCumulativePolicyType(irfpriceexp%cumSTICKY)
	END IF
END IF	


write(*,*) 'Finished allocating arrays'

END SUBROUTINE