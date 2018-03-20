MODULE Parameters
IMPLICIT NONE
SAVE

integer, parameter	:: ngpy1	= 3	 !must be odd
integer, parameter	:: ngpy2	= 11
integer, parameter 	:: nsim = 50000
real(8), parameter  :: dt 	= 0.25					!time step in quarters 
integer, parameter 	:: Tsim = FLOOR(20.0/dt)+1		!only need 20 quarters

END MODULE Parameters
