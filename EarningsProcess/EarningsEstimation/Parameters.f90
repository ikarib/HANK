MODULE Parameters
IMPLICIT NONE
SAVE

integer, parameter 	:: nsim = 5000 !50000
real(8), parameter  :: dt = 0.25 !0.005					!time step in quarters 
integer, parameter 	:: Tburn = FLOOR(100.0/dt)+1		
integer, parameter 	:: Tsim = Tburn + FLOOR(20.0/dt)+1		!only need 20 quarters
integer, parameter 	:: Tann = (Tsim-Tburn)*dt/4

END MODULE Parameters
