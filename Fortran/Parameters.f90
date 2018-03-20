MODULE Parameters
IMPLICIT NONE
SAVE

!OPTIONS
integer,parameter  ::  TwoPointWageProcess 		= 0	!set ngpy=2
integer,parameter  ::  Borrowing				= 1

!GRID SIZES
integer,parameter :: ngpa		= 40 !2		    !grid for illiquid assets
integer,parameter :: ngpbPOS	= 40 !60		    !grid for liquid assets, >=0 range
integer,parameter :: ngpbNEG	= 10		    !grid for liquid assets, <0 range only active if Borrowing==1
integer,parameter :: ngpb		= ngpbPOS + Borrowing*ngpbNEG
integer,parameter :: ngpy		= 33 !30
integer,parameter :: naby		= ngpa*ngpb*ngpy
integer,parameter :: nab		= ngpa*ngpb

!PARAMETERS FOR GRID CONSTRUCTION
real(8), parameter   :: agridparam = 0.15 !0.9 !0.15		!for a: approaches linear as goes to 1, approaches L shaped as goes to 0
real(8), parameter   :: bgridparam = 0.35 !0.25 !0.35		!for b pos: approaches linear as goes to 1, approaches L shaped as goes to 0
real(8), parameter   :: bgridparamNEG = 0.4		!for b neg: approaches linear as goes to 1, approaches L shaped as goes to 0
real(8), parameter   :: amax  = 2000.0 !0.1 !2000.0			!multiple of quarterly output
real(8), parameter   :: bmax  = 40.0 !500.0 !40.0

!OTHER PARAMETERS
real(8), parameter 	 :: cmin = 1.0e-5	!minimum consumption for natural borrowing limit
real(8), parameter 	 :: dmax = 1.0e10	!maximum deposit rate, for numerical stability while converging
real(8), parameter 	 :: facc = 1.0e-10 !1.0e-6

integer, parameter   :: Ttransition = 200 !no. time steps for the transition (each step is can be a different number of time units)


END MODULE Parameters
