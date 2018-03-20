REAL(8) FUNCTION  FnCapitalEquity(lcapital)

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

REAL(8), INTENT(IN) :: lcapital

IF(DistributeProfitsInProportion==0) FnCapitalEquity = (1.0-fundlev)*(lcapital*elast*ra + (1.0-corptax)*tfp*(lcapital**alpha)*(labor**(1.0-alpha))) - Ea*elast*ra 
IF(DistributeProfitsInProportion==1) FnCapitalEquity = (1.0-fundlev)*(lcapital*elast*ra + profdistfrac*(1.0-corptax)*tfp*(lcapital**alpha)*(labor**(1.0-alpha))) - Ea*elast*ra 

END FUNCTION FnCapitalEquity