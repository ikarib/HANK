SUBROUTINE MakeGuess

USE Parameters
USE Globals

IMPLICIT NONE

y1widthguess 	= 2.0*sigma1
y1widthmin		= 1.0*sigma1
y1widthmax 		= 10.0*sigma1

y1gridparguess 	= 0.7	!1 for linear, 0 for all points in middle, defaults to 1 if not estimated and DefaultY1GridPar1=1
y1gridmin 		= 0.3 !0.1

y2widthguess 	= 2.0*sigma2
y2widthmin		= 1.0*sigma2
y2widthmax 		= 10.0*sigma2

y2gridparguess 	= 0.7	!1 for linear, 0 for all points in middle, defaults to 1 if not estimated and DefaultY2GridPar1=1
y2gridmin 		= 0.3 !0.1


END SUBROUTINE MakeGuess