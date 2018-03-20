SUBROUTINE OptimalConsumption(lVb,lc,lh,ls,lHc)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8),INTENT(IN)	:: lVb
REAL(8),INTENT(OUT)	:: lc,lh,ls,lHc
REAL(8)	:: llabdisutil,lhmax,lhmin
INTEGER	:: iflag
REAL(8), EXTERNAL :: FnHoursBC

!set lVb = -999.9 to use budget constraint rather than FOC for consumption


IF (lVb > -999.0) THEN !not at stationary point or limit
	IF(NoLaborSupply==1) THEN
		lc = utilfn1inv(lVb)
		lh = 1.0
		ls = gbdrift + lh*gnetwage -lc
		IF(lc>0.0) THEN
			lHc = utilfn(lc) + lVb*ls
		ELSE
			lHc = -1.0e12
		END IF

	ELSE IF	(LaborSupplySep==1) THEN
		IF(ScaleDisutilityIdio==0) lh = (gnetwage*lVb/chi)**(frisch)
		IF(ScaleDisutilityIdio==1) lh = ((gnetwage/gidioprod)*lVb/chi)**(frisch)
		IF(ImposeMaxHours==1)lh = min(lh, 1.0_8)
		IF(ScaleDisutilityIdio==0) llabdisutil = lh**(1.0+1.0/frisch)/(1.0+1.0/frisch)
		IF(ScaleDisutilityIdio==1) llabdisutil = gidioprod * lh**(1.0+1.0/frisch)/(1.0+1.0/frisch)
		lc = utilfn1inv(lVb)
		ls = gbdrift + lh*gnetwage -lc
		IF(lc>0.0) THEN
			lHc = utilfn(lc) - chi*llabdisutil + lVb*ls
		ELSE
			lHc = -1.0e12
		END IF
		
	ELSE IF	(LaborSupplyGHH==1) THEN
		IF(ScaleDisutilityIdio==0) lh = (gnetwage/chi)**(frisch)
		IF(ScaleDisutilityIdio==1) lh = ((gnetwage/gidioprod)/chi)**(frisch)
		IF(ImposeMaxHours==1) lh = min(lh, 1.0_8)
		IF(ScaleDisutilityIdio==0) llabdisutil = lh**(1.0+1.0/frisch)/(1.0+1.0/frisch)
		IF(ScaleDisutilityIdio==1) llabdisutil = gidioprod * lh**(1.0+1.0/frisch)/(1.0+1.0/frisch)
		lc = utilfn1inv(lVb) + chi*llabdisutil
		ls = gbdrift + lh*gnetwage -lc
		IF(lc-chi*llabdisutil>0.0) THEN
			lHc = utilfn(lc-chi*llabdisutil) + lVb*ls
		ELSE
			lHc = -1.0e12
		END IF
		
	END IF	
	
ELSEIF (lVb <= -999.0) THEN !at stationary point or limit
	IF(NoLaborSupply==1) THEN
		lh = 1.0
		lc = gbdrift + lh*gnetwage
		ls = 0.0
		IF(lc>0.0) THEN
			lHc = utilfn(lc)
		ELSE
			lHc = -1.0e12
		END IF

	ELSE IF	(LaborSupplySep==1) THEN		
		ls = 0.0
		lhmin = max(0.0_8,-gbdrift/gnetwage+1.0e-5)
		IF(ImposeMaxHours==1) lhmax = 1.0
		IF(ImposeMaxHours==0) lhmax = 100.0
		IF (FnHoursBC(lhmin) >= 0.0) lh = lhmin
		IF (FnHoursBC(lhmax) <= 0.0) lh = lhmax
		IF (FnHoursBC(lhmin)<0.0 .and. FnHoursBC(lhmax)>0.0) THEN		
			CALL rtsec(FnHoursBC,lhmin,lhmax,1.0e-8_8,lh,iflag)
			IF(iflag<0) CALL rtbis(FnHoursBC,lhmin,lhmax,1.0e-8_8,1.0e-10_8,lh) 
		END IF
		IF(ScaleDisutilityIdio==0) llabdisutil = lh**(1.0+1.0/frisch)/(1.0+1.0/frisch)
		IF(ScaleDisutilityIdio==1) llabdisutil = gidioprod * lh**(1.0+1.0/frisch)/(1.0+1.0/frisch)
		
		lc = gbdrift + lh*gnetwage
		IF(lc>0.0) THEN
			lHc = utilfn(lc) - chi*llabdisutil
		ELSE
			lHc = -1.0e12
		END IF
		
		
	ELSE IF	(LaborSupplyGHH==1) THEN
		
		ls = 0.0
		IF(ScaleDisutilityIdio==0) lh = (gnetwage/chi)**(frisch)
		IF(ScaleDisutilityIdio==1) lh = ((gnetwage/gidioprod)/chi)**(frisch)
		
		IF(ImposeMaxHours==1) lh = min(lh, 1.0_8)
		IF(ScaleDisutilityIdio==0) llabdisutil = lh**(1.0+1.0/frisch)/(1.0+1.0/frisch)
		IF(ScaleDisutilityIdio==1) llabdisutil = gidioprod * lh**(1.0+1.0/frisch)/(1.0+1.0/frisch)
		
		lc = gbdrift + lh*gnetwage
		
		IF(lc-chi*llabdisutil>0.0) THEN
			lHc = utilfn(lc-chi*llabdisutil)
		ELSE
			lHc = -1.0e12
		END IF
		
	END IF	
	
	
	
END IF		



END SUBROUTINE OptimalConsumption