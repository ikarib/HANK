REAL(8) FUNCTION  FnHoursBC(lh)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lh

IF(LaborSupplySep==1) THEN
	IF(ScaleDisutilityIdio==0) FnHoursBC = chi*(lh**(1.0/frisch)) - utilfn1(gbdrift + lh*gnetwage) * gnetwage
	IF(ScaleDisutilityIdio==1) FnHoursBC = chi*(lh**(1.0/frisch)) - utilfn1(gbdrift + lh*gnetwage) * (gnetwage/gidioprod)
END IF

END FUNCTION FnHoursBC