MODULE Procedures

USE Parameters
USE Globals

IMPLICIT NONE
CONTAINS

!-------------------------------------------------------

REAL(8) FUNCTION DoubleParetoPDF(lx,lmu,lzeta) 
REAL(8), INTENT(IN)     :: lx,lmu,lzeta

IF(lx<=lmu) THEN
	DoubleParetoPDF = (lzeta/2.0)*exp(lzeta*(lx-lmu))
ELSE IF(lx>lmu) THEN
	DoubleParetoPDF = (lzeta/2.0)*exp(-lzeta*(lx-lmu))
END IF

END FUNCTION DoubleParetoPDF

!-------------------------------------------------------

REAL(8) FUNCTION DoubleParetoCDF(lx,lmu,lzeta) 
REAL(8), INTENT(IN)     :: lx,lmu,lzeta

IF(lx<=lmu) THEN
	DoubleParetoCDF = 0.5*exp(lzeta*(lx-lmu))
ELSE IF(lx>lmu) THEN
	DoubleParetoCDF = 1.0 - 0.5*exp(-lzeta*(lx-lmu))
END IF

END FUNCTION DoubleParetoCDF

!-------------------------------------------------------

REAL(8) FUNCTION DoubleParetoInverseCDF(lu,lmu,lzeta) 
REAL(8), INTENT(IN)     :: lu,lmu,lzeta

IF(lu<=0.5) THEN
	DoubleParetoInverseCDF = lmu + (1.0/lzeta)*log(2.0*lu)
ELSE IF(lu>0.5) THEN
	DoubleParetoInverseCDF = lmu - (1.0/lzeta)*log(2.0*(1.0-lu)) 
END IF

END FUNCTION DoubleParetoInverseCDF

!-------------------------------------------------------

REAL(8) FUNCTION DoubleParetoSkewPDF(lx,lmu,lzetaP,lzetaN) 
REAL(8), INTENT(IN)     :: lx,lmu,lzetaP,lzetaN
REAL(8)					:: lc

!lzetaN<0<lzetaP
lc = 1.0/(1.0/lzetaP-1.0/lzetaN)

IF(lx<=lmu) THEN
	DoubleParetoSkewPDF = lc*exp(-lzetaN*(lx-lmu))
ELSE IF(lx>lmu) THEN
	DoubleParetoSkewPDF = lc*exp(-lzetaP*(lx-lmu))
END IF

END FUNCTION DoubleParetoSkewPDF

!-------------------------------------------------------

REAL(8) FUNCTION DoubleParetoSkewCDF(lx,lmu,lzetaP,lzetaN) 
REAL(8), INTENT(IN)     :: lx,lmu,lzetaP,lzetaN
REAL(8)					:: lc

!lzetaN<0<lzetaP
lc = 1.0/(1.0/lzetaP-1.0/lzetaN)

IF(lx<=lmu) THEN
	DoubleParetoSkewCDF = -(lc/lzetaN)*exp(-lzetaN*(lx-lmu))
ELSE IF(lx>lmu) THEN
	DoubleParetoSkewCDF = -(lc/lzetaN)-(lc/lzetaP)*(exp(-lzetaP*(lx-lmu)) -1.0)
END IF

END FUNCTION DoubleParetoSkewCDF
!-------------------------------------------------------

REAL(8) FUNCTION TruncatedDoubleParetoSkewCDF(lx,lmu,lzetaP,lzetaN,lxL,lxU) 
REAL(8), INTENT(IN)     :: lx,lmu,lzetaP,lzetaN,lxL,lxU
REAL(8)					:: lFL,lFU,lFx

lFL = DoubleParetoSkewCDF(lxL,lmu,lzetaP,lzetaN) 
lFU = DoubleParetoSkewCDF(lxU,lmu,lzetaP,lzetaN) 
lFx = DoubleParetoSkewCDF(lx,lmu,lzetaP,lzetaN) 

TruncatedDoubleParetoSkewCDF = (lFx-lFL)/(lFU-lFL)

END FUNCTION TruncatedDoubleParetoSkewCDF

!-------------------------------------------------------

REAL(8) FUNCTION DoubleParetoSkewPDFv2(lx,lmu,lzetaP,lzetaN) 
REAL(8), INTENT(IN)     :: lx,lmu,lzetaP,lzetaN
REAL(8)					:: lcN,lcP

!lzetaN<0<lzetaP
lcN = -lzetaN/2.0
lcP = lzetaP/2.0

IF(lx<=lmu) THEN
	DoubleParetoSkewPDFv2 = lcN*exp(-lzetaN*(lx-lmu))
ELSE IF(lx>lmu) THEN
	DoubleParetoSkewPDFv2 = lcP*exp(-lzetaP*(lx-lmu))
END IF

END FUNCTION DoubleParetoSkewPDFv2

!-------------------------------------------------------

REAL(8) FUNCTION DoubleParetoSkewCDFv2(lx,lmu,lzetaP,lzetaN) 
REAL(8), INTENT(IN)     :: lx,lmu,lzetaP,lzetaN
REAL(8)					:: lcN,lcP

!lzetaN<0<lzetaP
lcN = -lzetaN/2.0
lcP = lzetaP/2.0

IF(lx<=lmu) THEN
	DoubleParetoSkewCDFv2 = 0.5*exp(-lzetaN*(lx-lmu))
ELSE IF(lx>lmu) THEN
	DoubleParetoSkewCDFv2 = 0.5 -0.5*(exp(-lzetaP*(lx-lmu))-1.0)
END IF

END FUNCTION DoubleParetoSkewCDFv2
!-------------------------------------------------------

REAL(8) FUNCTION TruncatedDoubleParetoSkewCDFv2(lx,lmu,lzetaP,lzetaN,lxL,lxU) 
REAL(8), INTENT(IN)     :: lx,lmu,lzetaP,lzetaN,lxL,lxU
REAL(8)					:: lFL,lFU,lFx

lFL = DoubleParetoSkewCDFv2(lxL,lmu,lzetaP,lzetaN) 
lFU = DoubleParetoSkewCDFv2(lxU,lmu,lzetaP,lzetaN) 
lFx = DoubleParetoSkewCDFv2(lx,lmu,lzetaP,lzetaN) 

TruncatedDoubleParetoSkewCDFv2 = (lFx-lFL)/(lFU-lFL)

END FUNCTION TruncatedDoubleParetoSkewCDFv2

!-------------------------------------------------------
REAL(8) FUNCTION NormalCDF(lx,lmu,lsigma) 
REAL(8), INTENT(IN)     :: lx,lmu,lsigma
REAL(8) 	:: lc,ltemp

CALL cumnor((lx-lmu)/lsigma,lc,ltemp)
NormalCDF = lc

END FUNCTION NormalCDF
!-------------------------------------------------------

REAL(8) FUNCTION TruncatedNormalCDF(lx,lmu,lsigma,lxL,lxU) 
REAL(8), INTENT(IN)     :: lx,lmu,lsigma,lxL,lxU
REAL(8)					:: lFL,lFU,lFx

lFL = NormalCDF(lxL,lmu,lsigma) 
lFU = NormalCDF(lxU,lmu,lsigma) 
lFx = NormalCDF(lx,lmu,lsigma) 

TruncatedNormalCDF = (lFx-lFL)/(lFU-lFL)

END FUNCTION TruncatedNormalCDF

!-------------------------------------------------------
!-------------------------------------------------------

REAL(8) FUNCTION logistic(lx) 
REAL(8), INTENT(IN)     :: lx

logistic = 1.0/(1.0+dexp(-lx))

END FUNCTION logistic
!-------------------------------------------------------
REAL(8) FUNCTION invlogistic(lx) 
REAL(8), INTENT(IN)     :: lx

if (lx<=0.0 .or. lx>=1) then
	write(*,*) 'error in invlogistic'
	return
else
	invlogistic = -dlog(1.0/lx - 1.0)
end if

END FUNCTION invlogistic

!-----------------------------------------------

SUBROUTINE LinInterp (n,x,y,ni,xi,yi)
!this does linear interpolation of (x,y) at points xi
!requires x to be sorted in ascending order
!extrapolates out of range
IMPLICIT NONE
INTEGER, INTENT(in)		:: n,ni
REAL(8), INTENT(in)	:: x(:),y(:),xi(:)
REAL(8), INTENT(out)	:: yi(:)
REAL(8), DIMENSION(ni)	:: xL,xH,yL,yH
INTEGER					:: i,locL(ni)


DO i = 1,ni

	LocL(i) = MAXLOC(x,1,MASK=xi(i)>x)
	
	IF (xi(i)<=x(1)) THEN
		LocL(i) = 1
	END IF

	IF (LocL(i)>=n) THEN 
		LocL(i) = n-1
	END IF


	xL(i) = x(locL(i))
	xH(i) = x(locL(i)+1)
	yL(i) = y(locL(i))
	yH(i) = y(locL(i)+1)
	
	yi(i) = yL(i) + (xi(i)-xL(i))*((yH(i)-yL(i))/(xH(i)-xL(i)))

END DO


END SUBROUTINE LinInterp

!-----------------------------------------------

SUBROUTINE LinInterp1 (n,x,y,xi,yi)
!this does linear interpolation of (x,y) at points only point,xi
!requires x to be sorted in ascending order
!extrapolates out of range
IMPLICIT NONE
INTEGER, INTENT(in)		:: n
REAL(8), INTENT(in)	:: x(:),y(:),xi
REAL(8), INTENT(out)	:: yi
REAL(8)	            :: xL,xH,yL,yH
INTEGER					:: locL

LocL = MAXLOC(x,1,MASK=xi>x)

IF (xi<=x(1)) THEN
	LocL = 1
END IF

IF (LocL>=n) THEN 
	LocL = n-1
END IF

xL  = x(locL)
xH  = x(locL +1)
yL  = y(locL)
yH  = y(locL +1)

yi  = yL  + (xi -xL )*((yH -yL )/(xH -xL ))

END SUBROUTINE LinInterp1

!-----------------------------------------------

SUBROUTINE BiLinInterp1 (nx,x,ny,y,f,xi,yi,fi)
!this does linear interpolation of f(x,y) at points (xi,yi)
!requires x and y to be sorted in ascending order
!extrapolates out of range
IMPLICIT NONE
INTEGER, INTENT(in)		:: nx,ny
REAL(8), INTENT(in)	:: x(:),y(:),f(:,:),xi,yi
REAL(8), INTENT(out)	:: fi
REAL(8)	            :: xL,xH,yL,yH,fLL,fHH,fLH,fHL,dxdy
INTEGER					:: xlocL,ylocL

xlocL = MAXLOC(x,1,MASK=xi>x)
ylocL = MAXLOC(y,1,MASK=yi>y)

IF (xi<=x(1)) THEN
	xlocL = 1
END IF

IF (xLocL>=nx) THEN 
	xLocL = nx-1
END IF

IF (yi<=y(1)) THEN
	ylocL = 1
END IF

IF (yLocL>=ny) THEN 
	yLocL = ny-1
END IF

xL  = x(xlocL)
xH  = x(xlocL +1)
yL  = y(ylocL)
yH  = y(ylocL +1)
fLL = f(xlocL,ylocL)
fLH = f(xlocL,ylocL+1)
fHL = f(xlocL+1,ylocL)
fHH = f(xlocL+1,ylocL+1)

dxdy = (xH-xL)*(yH-yL)
fi = fLL*(xH-xi)*(yH-yi)/(dxdy) + fHL*(xi-xL)*(yH-yi)/(dxdy) + fLH*(xH-xi)*(yi-yL)/(dxdy) + fHH*(xi-xL)*(yi-yL)/(dxdy)


END SUBROUTINE BiLinInterp1

!-----------------------------------------------

SUBROUTINE Bi45DegInterp1 (nx,x,ny,y,f,xi,yi,fi)
!this does linear interpolation of f(x,y) at points (xi,yi)
!requires x and y to be sorted in ascending order
!extrapolates out of range
!interpolates along the 45 degree line xi + yi = K

IMPLICIT NONE
INTEGER, INTENT(in)		:: nx,ny
REAL(8), INTENT(in)	:: x(:),y(:),f(:,:),xi,yi
REAL(8), INTENT(out)	:: fi
REAL(8)	            :: xL,xH,yL,yH,fLL,fHH,fLH,fHL,K,ficheck,lyi
REAL(8)	            :: x1,y1,f1,x2,y2,f2
INTEGER					:: xlocL,ylocL,check

 check = 0
lyi = yi


K = xi+lyi

xlocL = MAXLOC(x,1,MASK=xi>x)
ylocL = MAXLOC(y,1,MASK=lyi>y)

IF (xi<=x(1)) THEN
	xlocL = 1
END IF

IF (xLocL>=nx) THEN 
	xLocL = nx-1
END IF

IF (lyi<=y(1)) THEN
	ylocL = 1
END IF

IF (yLocL>=ny) THEN 
	yLocL = ny-1
END IF

xL  = x(xlocL)
xH  = x(xlocL +1)
yL  = y(ylocL)
yH  = y(ylocL +1)
fLL = f(xlocL,ylocL)
fLH = f(xlocL,ylocL+1)
fHL = f(xlocL+1,ylocL)
fHH = f(xlocL+1,ylocL+1)
!check corners
IF(xi==xL .and. lyi==yL) THEN
	fi = fLL
	return
END IF
IF(xi==xH .and. lyi==yH) THEN
	fi = fHH
	return
END IF



!DA or AB
IF((K-xL<=yH .and. K-xL>=yL) .or. xL+yL>K) THEN !DA
	x1 = xL
	y1 = K-xL
	f1  = fLL  + (y1 -yL )*((fLH -fLL )/(yH -yL ))
ELSE IF ((K-yH<=xH .and. K-yH>=xL) .or. xH+yH<K) THEN!AB
	x1 = K-yH
	y1 = yH
	f1  = fLH  + (x1 -xL )*((fHH -fLH )/(xH -xL ))
ELSE
	IF(Display>=2) write(*,*) 'problem 1 in Bi45DegInterp1'
END IF


!BC or CD
IF((K-xH<=yH .and. K-xH>=yL) .or. xH +yH<K) THEN !BC
	x2 = xH
	y2 = K-xH
	f2  = fHL  + (y2 -yL )*((fHH -fHL )/(yH -yL ))
ELSE IF ((K-yL<=xH .and. K-yL>=xL) .or. xL+yL>K ) THEN!CD
	x2 = K-yL
	y2 = yL
	f2  = fLL  + (x2 -xL )*((fHL -fLL )/(xH -xL ))
ELSE
	IF(Display>=2) write(*,*) 'problem 2 in Bi45DegInterp1'
END IF

!interpolate (x1,y1,f1) and (x2,y2,f2) at (xi,yi)
fi  = f1  + (xi -x1 )*((f2 -f1 )/(x2 -x1 ))

!check
ficheck = f1  + (lyi -y1 )*((f2 -f1 )/(y2 -y1 ))

IF (abs(fi-ficheck)>1.0e-10) THEN
	IF(Display>=2) write(*,*) 'problem check in Bi45DegInterp1: ', abs(fi-ficheck)
END IF


END SUBROUTINE Bi45DegInterp1

!-----------------------------------------------

SUBROUTINE BiLinInterp (nx,x,ny,y,f,xi,yi,fi)
!this does linear interpolation of f(x,y) at points (xi,yi)
!requires x and y to be sorted in ascending order
!xi and yi are vectors of the same length ni
!extrapolates out of range
IMPLICIT NONE
INTEGER, INTENT(in)		:: nx,ny
REAL(8), INTENT(in)	:: x(:),y(:),f(:,:),xi(:),yi(:)
REAL(8), INTENT(out)	:: fi(:)
REAL(8)	            :: xL,xH,yL,yH,fLL,fHH,fLH,fHL,dxdy
INTEGER					:: xlocL,ylocL,i,ni

ni = size(xi)

DO i = 1,ni
	xlocL = MAXLOC(x,1,MASK=xi(i)>x)
	ylocL = MAXLOC(y,1,MASK=yi(i)>y)

	IF (xi(i)<=x(1)) THEN
		xlocL = 1
	END IF

	IF (xLocL>=nx) THEN 
		xLocL = nx-1
	END IF

	IF (yi(i)<=y(1)) THEN
		ylocL = 1
	END IF

	IF (yLocL>=ny) THEN 
		yLocL = ny-1
	END IF

	xL  = x(xlocL)
	xH  = x(xlocL +1)
	yL  = y(ylocL)
	yH  = y(ylocL +1)
	fLL = f(xlocL,ylocL)
	fLH = f(xlocL,ylocL+1)
	fHL = f(xlocL+1,ylocL)
	fHH = f(xlocL+1,ylocL+1)

	dxdy = (xH-xL)*(yH-yL)
	fi(i) = fLL*(xH-xi(i))*(yH-yi(i))/(dxdy) + fHL*(xi(i)-xL)*(yH-yi(i))/(dxdy) + fLH*(xH-xi(i))*(yi(i)-yL)/(dxdy) + fHH*(xi(i)-xL)*(yi(i)-yL)/(dxdy)
END DO

END SUBROUTINE BiLinInterp

!-----------------------------------------------

SUBROUTINE PowerSpacedGrid (n,k,low,high,y)

!gives a grid spaced between low and high based on the unit interval with a function x^(1/k)
!k = 1 is linear, k = 0 is L-shaped

IMPLICIT NONE

INTEGER, INTENT(in)	:: n
REAL(8), INTENT(in)	:: k,low,high
REAL(8), INTENT(out) :: y(:)
INTEGER				:: i
REAL(8)				:: x(n),z(n)

IF(n<2) THEN
	write(*,*) 'n must be at least 2 to make grids'
	return
END IF

IF (n==2) THEN
	y(1) = low
	y(2) = high
	return
END IF

x(1) = 0.0
x(n) = 1.0
DO i = 2,n-1
	x(i) = (i-1)/real(n-1)
END DO

z = x**(1.0/k)

y = low + (high-low)*z


END SUBROUTINE PowerSpacedGrid

!-----------------------------------------------

SUBROUTINE LinearSpacedGrid (n,low,high,y)

!gives a linear spaced between low and high
IMPLICIT NONE

INTEGER, INTENT(in)	:: n
REAL(8), INTENT(in)	:: low,high
REAL(8), INTENT(out) :: y(:)
INTEGER				:: i
REAL(8)				:: x(n)


IF (n==2) THEN
	y(1) = low
	y(2) = high
	return
END IF

x(1) = 0.0
x(n) = 1.0
DO i = 2,n-1
	x(i) = (i-1)/real(n-1)
END DO

y = low + (high-low)*x


END SUBROUTINE LinearSpacedGrid
!-----------------------------------------------

SUBROUTINE FindLinProb1 (xi,x,y,p)
!this takes in xi, and finds two points in either side of in in x
!and returns the indices of them y and associated probabilities p
IMPLICIT NONE
REAL(8), INTENT(in)	    :: x(:),xi
REAL(8), INTENT(out)	:: p(2)
INTEGER, INTENT(out)    :: y(2)

INTEGER					:: locL,n

n = size(x)
IF (n==1) THEN
	p(1) = 1.0
	p(2) = 0.0
	y(1) = 1
	y(2) = 1
	return
END IF

LocL = MAXLOC(x,1,MASK=xi>x)

IF (xi<=x(1)) THEN
	y(1) = 1
	y(2) = 2
	p(1) = 1.0
	p(2) = 0.0
ELSE IF (LocL>=n) THEN 
	LocL = n-1
    y(1) = n-1
    y(2) = n
    p(1) = 0.0
    p(2) = 1.0
ELSE IF (x(LocL+1) == x(LocL)) THEN
	y(1) = LocL
    y(2) = LocL+1
	p(1) = 0.5
	p(2) = 0.5
ELSE
    y(1) = LocL
    y(2) = LocL+1
    p(2) = (xi - x(LocL))/real(x(LocL+1)-x(LocL))
    p(1) = 1.0 - p(2)
END IF

END SUBROUTINE FindLinProb1

!----------------------------------------------
SUBROUTINE WriteMatrix(f,n1,n2,mat)

INTEGER, INTENT(IN)			:: f,n1,n2
REAL(8),INTENT(in)		:: mat(n1,n2)
CHARACTER		::lstring*80
INTEGER			::i1

WRITE(UNIT=lstring, FMT='(I5)') n2
lstring = '('//trim(lstring) // 'F16.6)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrix

!----------------------------------------------
SUBROUTINE WriteMatrixLong(f,n1,n2,mat)

INTEGER, INTENT(IN)			:: f,n1,n2
REAL(8),INTENT(in)		:: mat(n1,n2)
CHARACTER		::lstring*80
INTEGER			::i1

WRITE(UNIT=lstring, FMT='(I5)') n2
lstring = '('//trim(lstring) // 'F20.14)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrixLong

!----------------------------------------------
SUBROUTINE WriteMatrixInteger(f,n1,n2,mat)

INTEGER, INTENT(IN)			:: f,n1,n2
INTEGER,INTENT(in)		:: mat(n1,n2)
CHARACTER		::lstring*80
INTEGER			::i1

WRITE(UNIT=lstring, FMT='(I5)') n2
lstring = '('//trim(lstring) // 'I16)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrixInteger

!----------------------------------------------
SUBROUTINE WriteMatrixExpon(f,n1,n2,mat)

INTEGER, INTENT(IN)			:: f,n1,n2
REAL(8),INTENT(in)		:: mat(n1,n2)
CHARACTER		::lstring*80
INTEGER			::i1

WRITE(UNIT=lstring, FMT='(I5)') n2
lstring = '('//trim(lstring) // 'E16.8)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrixExpon

!----------------------------------------------

SUBROUTINE WriteMatrixCSV(f,n1,n2,mat)

INTEGER, INTENT(IN)			:: f,n1,n2
REAL(8),INTENT(in)		:: mat(n1,n2)
CHARACTER		::lstring*80
INTEGER			::i1

WRITE(UNIT=lstring, FMT='(I5)') n2-1
lstring = '('//trim(lstring) // '(F20.14,","),F20.14)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrixCSV
!----------------------------------------------

SUBROUTINE MASmooth(n,x,v,d,s)

IMPLICIT NONE

INTEGER, INTENT(IN)					 :: n,d
REAL(8), INTENT(IN), DIMENSION(n) 	 :: x,v
REAL(8), INTENT(OUT), DIMENSION(n) 	 :: s
INTEGER								:: i

!n is length
!x is input vector
!v is vector of time steps
!s is smoothed vector
!d is size of MA filter (number of points either side so d=1 is a 3 point MA filter)

DO i = 1,d
	s(i) = SUM(x(1:i+d)*v(1:i+d))/SUM(v(1:i+d))
END DO
DO i = d+1,n-d
	s(i) = SUM(x(i-d:i+d)*v(i-d:i+d))/SUM(v(i-d:i+d))
END DO
DO i = n-d+1,n
	s(i) = SUM(x(i-d:n)*v(i-d:n))/SUM(v(i-d:n))	
END DO


END SUBROUTINE

!----------------------------------------------
SUBROUTINE DiscreteDist(Nout,Xout,Nin,Pin,lran)
!takes in random numbers and a prob dist over the integers 1 to Nin
!and returns the corresponding values
!Created 11/22/09 - ERS

INTEGER,    INTENT(out) :: Xout(:)
INTEGER,    INTENT(in)  :: Nin,Nout
REAL(8),    INTENT(in)  :: Pin(:)
REAL(8), INTENT(in)  	:: lran(:)
INTEGER			        :: i1,i2

Xout(:) = 0
DO i1 = 1,Nout
    IF ( lran(i1) .le. Pin(1) ) THEN
        Xout(i1) = 1
    ELSE
        i2 = 2
        DO WHILE (i2 .le. Nin)
            IF ( (lran(i1) .le. SUM(Pin(1:i2)) ).and. (lran(i1) > SUM(Pin(1:i2-1)) ) ) THEN
                Xout(i1) = i2
                i2 = Nin+1
            ELSE
                i2 = i2+1
            END IF
        END DO
    END IF
END DO

END SUBROUTINE DiscreteDist

!-----------------------------------------------
SUBROUTINE DiscreteDist1(Xout,Nin,Pin,lran)
!takes in a random number and a prob dist over the integers 1 to Nin
!and returns the corresponding value

INTEGER, INTENT(in)     :: Nin
INTEGER,INTENT(out)     :: Xout
REAL(8), INTENT(in)     :: Pin(Nin),lran
INTEGER                 :: i2

Xout = 0
IF ( lran .le. Pin(1) ) THEN
   Xout = 1
ELSE
   i2 = 2
   DO WHILE (i2 .le. Nin)
       IF ( (lran .le. SUM(Pin(1:i2)) ).and. (lran > SUM(Pin(1:i2-1)) ) ) THEN
           Xout = i2
           i2 = Nin+1
       ELSE
           i2 = i2+1
       END IF
   END DO
END IF

END SUBROUTINE DiscreteDist1


subroutine cumnor ( arg, cum, ccum )

!*****************************************************************************80
!
!! CUMNOR computes the cumulative normal distribution.
!
!  Discussion:
!
!    This function evaluates the normal distribution function:
!
!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo
!
!    This transportable program uses rational functions that
!    theoretically approximate the normal distribution function to
!    at least 18 significant decimal digits.  The accuracy achieved
!    depends on the arithmetic system, the compiler, the intrinsic
!    functions, and proper selection of the machine dependent
!    constants.
!
!  Author: 
!
!    William Cody
!    Mathematics and Computer Science Division
!    Argonne National Laboratory
!    Argonne, IL 60439
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev approximations for the error function,
!    Mathematics of Computation, 
!    1969, pages 631-637.
!
!    William Cody, 
!    Algorithm 715: 
!    SPECFUN - A Portable FORTRAN Package of Special Function Routines 
!    and Test Drivers,
!    ACM Transactions on Mathematical Software,
!    Volume 19, Number 1, 1993, pages 22-32.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the upper limit of integration.
!
!    Output, real ( kind = 8 ) CUM, CCUM, the Normal density CDF and
!    complementary CDF.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) EPS, the argument below which anorm(x) 
!    may be represented by 0.5 and above which  x*x  will not underflow.
!    A conservative value is the largest machine number X
!    such that   1.0D+00 + X = 1.0D+00   to machine precision.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 5 ) :: a = (/ &
    2.2352520354606839287D+00, &
    1.6102823106855587881D+02, &
    1.0676894854603709582D+03, &
    1.8154981253343561249D+04, &
    6.5682337918207449113D-02 /)
  real ( kind = 8 ) arg
  real ( kind = 8 ), parameter, dimension ( 4 ) :: b = (/ &
    4.7202581904688241870D+01, &
    9.7609855173777669322D+02, &
    1.0260932208618978205D+04, &
    4.5507789335026729956D+04 /)
  real ( kind = 8 ), parameter, dimension ( 9 ) :: c = (/ &
    3.9894151208813466764D-01, &
    8.8831497943883759412D+00, &
    9.3506656132177855979D+01, &
    5.9727027639480026226D+02, &
    2.4945375852903726711D+03, &
    6.8481904505362823326D+03, &
    1.1602651437647350124D+04, &
    9.8427148383839780218D+03, &
    1.0765576773720192317D-08 /)
  real ( kind = 8 ) ccum
  real ( kind = 8 ) cum
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    2.2266688044328115691D+01, &
    2.3538790178262499861D+02, &
    1.5193775994075548050D+03, &
    6.4855582982667607550D+03, &
    1.8615571640885098091D+04, &
    3.4900952721145977266D+04, &
    3.8912003286093271411D+04, &
    1.9685429676859990727D+04 /)
  real ( kind = 8 ) del
  real ( kind = 8 ) eps
  integer i
  real ( kind = 8 ), parameter, dimension ( 6 ) :: p = (/ &
    2.1589853405795699D-01, &
    1.274011611602473639D-01, &
    2.2235277870649807D-02, &
    1.421619193227893466D-03, &
    2.9112874951168792D-05, &
    2.307344176494017303D-02 /)
  real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
    1.28426009614491121D+00, &
    4.68238212480865118D-01, &
    6.59881378689285515D-02, &
    3.78239633202758244D-03, &
    7.29751555083966205D-05 /)
  real ( kind = 8 ), parameter :: root32 = 5.656854248D+00
  real ( kind = 8 ), parameter :: sixten = 16.0D+00
  real ( kind = 8 ) temp
  real ( kind = 8 ), parameter :: sqrpi = 3.9894228040143267794D-01
  real ( kind = 8 ), parameter :: thrsh = 0.66291D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xden
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) xsq
!
!  Machine dependent constants
!
  eps = epsilon ( 1.0D+00 ) * 0.5D+00

  x = arg
  y = abs ( x )

  if ( y <= thrsh ) then
!
!  Evaluate  anorm  for  |X| <= 0.66291
!
    if ( eps < y ) then
      xsq = x * x
    else
      xsq = 0.0D+00
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do
    cum = x * ( xnum + a(4) ) / ( xden + b(4) )
    temp = cum
    cum = 0.5D+00 + temp
    ccum = 0.5D+00 - temp
!
!  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
!
  else if ( y <= root32 ) then

    xnum = c(9) * y
    xden = y
    do i = 1, 7
      xnum = ( xnum + c(i) ) * y
      xden = ( xden + d(i) ) * y
    end do
    cum = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = aint ( y * sixten ) / sixten
    del = ( y - xsq ) * ( y + xsq )
    cum = exp ( - xsq * xsq * 0.5D+00 ) * exp ( -del * 0.5D+00 ) * cum
    ccum = 1.0D+00 - cum

    if ( 0.0D+00 < x ) then
      call r8_swap ( cum, ccum )
    end if
!
!  Evaluate ANORM for sqrt(32) < |X|.
!
  else

    cum = 0.0D+00
    xsq = 1.0D+00 / ( x * x )
    xnum = p(6) * xsq
    xden = xsq
    do i = 1, 4
      xnum = ( xnum + p(i) ) * xsq
      xden = ( xden + q(i) ) * xsq
    end do

    cum = xsq * ( xnum + p(5) ) / ( xden + q(5) )
    cum = ( sqrpi - cum ) / y
    xsq = aint ( x * sixten ) / sixten
    del = ( x - xsq ) * ( x + xsq )
    cum = exp ( - xsq * xsq * 0.5D+00 ) &
      * exp ( - del * 0.5D+00 ) * cum
    ccum = 1.0D+00 - cum

    if ( 0.0D+00 < x ) then
      call r8_swap ( cum, ccum )
    end if

  end if

  if ( cum < tiny ( cum ) ) then
    cum = 0.0D+00
  end if

  if ( ccum < tiny ( ccum ) ) then
    ccum = 0.0D+00
  end if

  return
end subroutine cumnor
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP switches two R8's.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end subroutine r8_swap

!!!!!!!!!!!!!!! FROM SparseKit !!!!!!!!!!!!!
subroutine apldia ( nrow, job, a, ja, ia, diag, b, jb, ib, iw )

!*****************************************************************************80
!
!! APLDIA adds a diagonal matrix to a general sparse matrix:  B = A + Diag.
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    The algorithm is in place (b, jb, ib, can be the same as
!    a, ja, ia, on entry). See comments for parameter job.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) JOB, job indicator. Job=0 means get array b only
!    (i.e. assume that a has already been copied into array b,
!    or that algorithm is used in place. ) For all practical
!    puposes enter job=0 for an in-place call and job=1 otherwise.
!    In case there are missing diagonal elements in A,
!    then the option job =0 will be ignored, since the algorithm
!    must modify the data structure (i.e. jb, ib) in this
!    situation.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(NROW), a diagonal matrix.
!
! on return:
!
! b,
! jb,
! ib      = resulting matrix B in compressed sparse row sparse format.
!
!
! iw    = integer ( kind = 4 ) work array of length n. On return iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k = 1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) diag(nrow)
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ib(nrow+1)
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iw(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jb(*)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) ko
  integer ( kind = 4 ) nnz
  logical test
!
!  Copy integer ( kind = 4 ) arrays into B's data structure if required.
!
  if ( job /= 0 ) then
    nnz = ia(nrow+1)-1
    jb(1:nnz) = ja(1:nnz)
    ib(1:nrow+1) = ia(1:nrow+1)
  end if
!
!  Get positions of diagonal elements in data structure.
!
  call diapos ( nrow, ja, ia, iw )
!
!  Count number of holes in diagonal and add DIAG elements to
!  valid diagonal entries.
!
  icount = 0

  do j = 1, nrow

     if ( iw(j) == 0 ) then
        icount = icount + 1
     else
        b(iw(j)) = a(iw(j)) + diag(j)
     end if

  end do
!
!  If no diagonal elements to insert, return.
!
  if ( icount == 0 ) then
    return
  end if
!
!  Shift the nonzero elements if needed, to allow for created
!  diagonal elements.
!
  ko = ib(nrow+1) + icount
!
!  Copy rows backward.
!
  do ii = nrow, 1, -1
!
!  Go through row II.
!
     k1 = ib(ii)
     k2 = ib(ii+1) - 1
     ib(ii+1) = ko
     test = ( iw(ii) == 0 )

     do k = k2, k1, -1

        j = jb(k)

         if ( test .and. j < ii ) then
           test = .false.
           ko = ko - 1
           b(ko) = diag(ii)
           jb(ko) = ii
           iw(ii) = ko
        end if

        ko = ko - 1
        b(ko) = a(k)
        jb(ko) = j

      end do
!
!  The diagonal element has not been added yet.
!
     if ( test ) then
        ko = ko - 1
        b(ko) = diag(ii)
        jb(ko) = ii
        iw(ii) = ko
     end if

  end do

  ib(1) = ko

  return
end subroutine apldia
subroutine coocsr ( nrow, nnz, a, ir, jc, ao, jao, iao )

!*****************************************************************************80
!
!! COOCSR converts COO to CSR.
!
!  Discussion:
!
!    This routine converts a matrix that is stored in COO coordinate format
!    a, ir, jc into a CSR row general sparse ao, jao, iao format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NNZ, the number of nonzero elements.
!
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!         the elements, ir(k) = its row number and jc(k) = its column
!        number. The order of the elements is arbitrary.
!
! on return:
!
! ir       is destroyed
!
!    Output, real AO(*), JAO(*), IAO(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iad
  integer ( kind = 4 ) iao(nrow+1)
  integer ( kind = 4 ) ir(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jao(*)
  integer ( kind = 4 ) jc(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k0
  integer ( kind = 4 ) nnz
  real ( kind = 8 ) x

  iao(1:nrow+1) = 0
!
!  Determine the row lengths.
!
  do k = 1, nnz
    iao(ir(k)) = iao(ir(k)) + 1
  end do
!
!  The starting position of each row.
!
  k = 1
  do j = 1, nrow+1
     k0 = iao(j)
     iao(j) = k
     k = k + k0
  end do
!
!  Go through the structure once more.  Fill in output matrix.
!
  do k = 1, nnz
     i = ir(k)
     j = jc(k)
     x = a(k)
     iad = iao(i)
     ao(iad) = x
     jao(iad) = j
     iao(i) = iad + 1
  end do
!
!  Shift back IAO.
!
  do j = nrow, 1, -1
    iao(j+1) = iao(j)
  end do
  iao(1) = 1

  return
end subroutine coocsr
subroutine diapos ( n, ja, ia, idiag )

!*****************************************************************************80
!
!! DIAPOS returns the positions of the diagonal elements of a sparse matrix.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, JA(*), IA(N+1), the matrix information, (but no values) 
!    in CSR Compressed Sparse Row format.
!
!    Output, integer ( kind = 4 ) IDIAG(N); the I-th entry of IDIAG points to the 
!    diagonal element A(I,I) in the arrays A and JA.  That is,
!    A(IDIAG(I)) = element A(I,I) of matrix A.  If no diagonal element 
!    is found, the entry is set to 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) idiag(n)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k

  idiag(1:n) = 0
!
!  Sweep through the data structure.
!
  do i = 1, n
    do k = ia(i), ia(i+1) -1
      if ( ja(k) == i ) then
        idiag(i) = k
      end if
    end do
  end do

  return
end subroutine diapos
subroutine transp ( nrow, ncol, a, ja, ia, iwk, ierr )

!*****************************************************************************80
!
!! TRANSP carries out in-place transposition routine.
!
!  Discussion:
!
!    This routine transposes a matrix stored in compressed sparse row
!    format.  The transposition is done in place in that the arrays 
!    A, JA, and IA of the transpose are overwritten onto the original arrays.
!
!    If you do not need the transposition to be done in place
!    it is preferrable to use the conversion routine csrcsc
!    (see conversion routines in formats).
!
!    The entries of the output matrix are not sorted (the column
!    indices in each are not in increasing order).  Use CSRCSC
!    if you want them sorted.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Workspace, integer ( kind = 4 ) IWK(*), of the same length as JA.
!
! on return:
!
!
! ncol      = actual row dimension of the transpose of the input matrix.
!         Note that this may be <= the input value for ncol, in
!         case some of the last columns of the input matrix are zero
!         columns. In the case where the actual number of rows found
!         in transp(A) exceeds the input value of ncol, transp will
!         return without completing the transposition. see ierr.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NCOL+1), the transposed
!    matrix in CSR Compressed Sparse Row format.
!
! ierr      = integer ( kind = 4 ). error message. If the number of rows for the
!         transposed matrix exceeds the input value of ncol,
!         then ierr is  set to that number and transp quits.
!         Otherwise ierr is set to 0 (normal return).
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) inext
  integer ( kind = 4 ) init
  integer ( kind = 4 ) iwk(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnz
  real ( kind = 8 ) t
  real ( kind = 8 ) t1

  ierr = 0
  nnz = ia(nrow+1) - 1
!
!  Determine the column dimension.
!
  jcol = 0
  do k = 1, nnz
    jcol = max ( jcol, ja(k) )
  end do

  if ( ncol < jcol ) then
     ierr = jcol
     return
  end if
!
!  Convert to coordinate format.  Use IWK for row indices.
!
  ncol = jcol

  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
      iwk(k) = i
    end do
  end do
!
!  Find pointer array for transpose.
!
  ia(1:ncol+1) = 0

  do k = 1, nnz
    i = ja(k)
    ia(i+1) = ia(i+1) + 1
  end do
  ia(1) = 1

  do i = 1, ncol
    ia(i+1) = ia(i) + ia(i+1)
  end do
!
!  Loop for a cycle in chasing process.
!
  init = 1
  k = 0

 5    continue

  t = a(init)
  i = ja(init)
  j = iwk(init)
  iwk(init) = -1

 6 continue

   k = k + 1
!
!  Current row number is I.  Determine where to go.
!
  l = ia(i)
!
!  Save the chased element.
!
  t1 = a(l)
  inext = ja(l)
!
!  Then occupy its location.
!
  a(l) = t
  ja(l) = j
!
!  Update pointer information for next element to be put in row i.
!
  ia(i) = l + 1
!
!  Determine next element to be chased.
!
  if ( iwk(l) < 0 ) then
    go to 65
  end if

  t = t1
  i = inext
  j = iwk(l)
  iwk(l) = -1

  if ( k < nnz ) then
    go to 6
  end if

  do i = ncol, 1, -1
    ia(i+1) = ia(i)
  end do

  ia(1) = 1

  return

 65   continue

  init = init + 1

  if ( nnz < init ) then

    do i = ncol, 1, -1
      ia(i+1) = ia(i)
    end do

    ia(1) = 1

    return
  end if

  if ( iwk(init) < 0 ) then
    go to 65
  end if
!
!  Restart chasing.
!
  go to 5
end subroutine transp
subroutine csort ( n, a, ja, ia, iwork, values )

!*****************************************************************************80
!
!! CSORT sorts the elements of a CSR matrix.
!
!  Discussion:
!
!    This routine sorts the elements of a CSR matrix (stored in Compressed
!    Sparse Row Format) in increasing order of their column indices within
!    each row. It uses a form of bucket sort with a cost of O(nnz) where
!    nnz = number of nonzero elements.
!
!    Requires an integer ( kind = 4 ) work array of size length 2*nnz.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! iwork = integer ( kind = 4 ) work array of length max ( n+1, 2*nnz )
!         where nnz = 2* (ia(n+1)-ia(1))  ) .
!
! values= logical indicating whether or not the real values a(*) must
!         also be permuted. if (.not. values) then the array a is not
!         touched by csort and can be a dummy array.
!
! on return:
!
! the matrix stored in the structure a, ja, ia is permuted in such a
! way that the column indices are in increasing order within each row.
! iwork(1:nnz) contains the permutation used  to rearrange the elements.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ifirst
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) iwork(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ko
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nnz
  logical values
!
!  Count the number of elements in each column.
!
  iwork(1:n+1) = 0

  do i = 1, n
     do k = ia(i), ia(i+1)-1
        j = ja(k) + 1
        iwork(j) = iwork(j) + 1
     end do
  end do
!
!  Compute pointers from lengths.
!
  iwork(1) = 1

  do i = 1, n
     iwork(i+1) = iwork(i) + iwork(i+1)
  end do
!
!  Get the positions of the nonzero elements in order of columns.
!
  ifirst = ia(1)
  nnz = ia(n+1)-ifirst

  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      next = iwork(j)
      iwork(nnz+next) = k
      iwork(j) = next + 1
    end do
  end do
!
!  Convert to coordinate format.
!
  do i = 1, n
    do k = ia(i), ia(i+1)-1
      iwork(k) = i
    end do
  end do
!
!  Loop to find permutation: for each element find the correct
!  position in (sorted) arrays A, JA.  Record this in IWORK.
!
  do k = 1, nnz
     ko = iwork(nnz+k)
     irow = iwork(ko)
     next = ia(irow)
!
!  The current element should go in next position in row. IWORK
!  records this position.
!
     iwork(ko) = next
     ia(irow) = next + 1
  end do
!
!  Perform an in-place permutation of the arrays.
!
     call ivperm ( nnz, ja(ifirst), iwork )

     if ( values ) then
       call dvperm ( nnz, a(ifirst), iwork )
     end if
!
!  Reshift the pointers of the original matrix back.
!
  do i = n, 1, -1
    ia(i+1) = ia(i)
  end do

  ia(1) = ifirst

  return
end subroutine csort
subroutine dvperm ( n, x, perm )

!*****************************************************************************80
!
!! DVPERM performs an in-place permutation of a real vector.
!
!  Discussion:
!
!    This routine permutes a real vector X using a permutation PERM.
!
!    On return, the vector X satisfies,
!
!      x(perm(j)) :== x(j), j = 1,2,.., n
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of X.
!
!    Input/output, real X(N), the vector to be permuted.
!
!    Input, integer ( kind = 4 ) PERM(N), the permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ii
  integer ( kind = 4 ) init
  integer ( kind = 4 ) k
  integer ( kind = 4 ) next
  integer ( kind = 4 ) perm(n)
  real ( kind = 8 ) tmp
  real ( kind = 8 ) tmp1
  real ( kind = 8 ) x(n)

  init = 1
  tmp = x(init)
  ii = perm(init)
  perm(init)= -perm(init)
  k = 0
!
!  The main loop.
!
 6  continue

   k = k + 1
!
!  Save the chased element.
!
  tmp1 = x(ii)
  x(ii) = tmp
  next = perm(ii)

  if ( next < 0 ) then
    go to 65
  end if
!
!  Test for end.
!
  if ( n < k ) then
    perm(1:n) = -perm(1:n)
    return
  end if

  tmp = tmp1
  perm(ii) = -perm(ii)
  ii = next
!
!  End of the loop.
!
  go to 6
!
!  Reinitialize cycle.
!
 65   continue

  init = init + 1

  if ( n < init ) then 
    perm(1:n) = -perm(1:n)
    return
  end if

  if ( perm(init) < 0 ) then
    go to 65
  end if

  tmp = x(init)
  ii = perm(init)
  perm(init) = -perm(init)
  go to 6

end subroutine dvperm
subroutine ivperm ( n, ix, perm )

!*****************************************************************************80
!
!! IVPERM performs an in-place permutation of an integer ( kind = 4 ) vector.
!
!  Discussion:
!
!    The integer ( kind = 4 ) vector ix is permuted according to the permutation 
!    array perm(*), i.e., on return, the vector x satisfies,
!
!      ix(perm(j)) :== ix(j), j = 1,2,.., n
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vector.
!
!    Input/output, integer ( kind = 4 ) IX(N), the vector to be permuted.
!
!    Input, integer ( kind = 4 ) PERM(N), the permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ii
  integer ( kind = 4 ) init
  integer ( kind = 4 ) ix(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) next
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) tmp
  integer ( kind = 4 ) tmp1

  init = 1
  tmp = ix(init)
  ii = perm(init)
  perm(init)= -perm(init)
  k = 0
!
!  Loop.
!
 6    continue

  k = k + 1
!
!  Save the chased element.
!
  tmp1 = ix(ii)
  ix(ii) = tmp
  next = perm(ii)

  if ( next < 0 ) then
    go to 65
  end if
!
!  Test for end.
!
  if ( n < k ) then
    perm(1:n) = -perm(1:n)
    return
  end if

  tmp = tmp1
  perm(ii) = -perm(ii)
  ii = next
!
!  End of loop.
!
  go to 6
!
!  Reinitilaize cycle.
!
 65   continue

  init = init + 1

  if ( n < init ) then
    perm(1:n) = -perm(1:n)
    return
  end if

  if ( perm(init) < 0 ) then
    go to 65
  end if

  tmp = ix(init)
  ii = perm(init)
  perm(init)=-perm(init)
  go to 6

end subroutine ivperm
subroutine csrcsc ( n, job, ipos, a, ja, ia, ao, jao, iao )
 
!*****************************************************************************80
!
!! CSRCSC converts Compressed Sparse Row to Compressed Sparse Column.
!
!  Discussion:
!
!    This is essentially a transposition operation.  
!
!    It is NOT an in-place algorithm.
!
!    This routine transposes a matrix stored in a, ja, ia format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) JOB, indicates whether or not to fill the values of the
!    matrix AO or only the pattern (IA, and JA).  Enter 1 for yes.
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use
!                call csrcsc (n,1,n+2,a,ja,ia,a,ja,ia(n+2))
!        for any other normal usage, enter ipos=1.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real AO(*), JAO(*), IAO(N+1), the matrix in CSC
!    Compressed Sparse Column format.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) iao(n+1)
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jao(*)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) next
!
!  Compute lengths of rows of A'.
!
  iao(1:n+1) = 0

  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k) + 1
      iao(j) = iao(j) + 1
    end do
  end do
!
!  Compute pointers from lengths.
!
  iao(1) = ipos
  do i = 1, n
    iao(i+1) = iao(i) + iao(i+1)
  end do
!
!  Do the actual copying.
!
  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      next = iao(j)
      if ( job == 1 ) then
        ao(next) = a(k)
      end if
      jao(next) = i
      iao(j) = next + 1
    end do
  end do
!
!  Reshift IAO and leave.
!
  do i = n, 1, -1
    iao(i+1) = iao(i)
  end do
  iao(1) = ipos

  return
end subroutine csrcsc
subroutine amudia ( nrow, job, a, ja, ia, diag, b, jb, ib )

!*****************************************************************************80
!
!! AMUDIA performs the matrix by matrix product B = A * Diag  (in place)
!
!  Discussion:
!
!    The column dimension of A is not needed.
!    The algorithm is "in place", so B can take the place of A.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) JOB, job indicator. Job=0 means get array b
!    only job = 1 means get b, and the integer ( kind = 4 ) arrays ib, jb.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(NROW), the diagonal matrix stored as a vector.
!
!    Output, B(*), JB(*), IB(NROW+1), the resulting matrix B in 
!    compressed sparse row sparse format.
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) diag(nrow)
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ib(nrow+1)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jb(*)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2

  do ii = 1, nrow
!
!  Scale each element.
!
    k1 = ia(ii)
    k2 = ia(ii+1) - 1
    do k = k1, k2
      b(k) = a(k) * diag(ja(k))
    end do

  end do

  if ( job == 0 ) then
    return
  end if

  ib(1) = ia(1)
  do ii = 1, nrow
    ib(ii) = ia(ii)
    do k = ia(ii), ia(ii+1)-1
      jb(k) = ja(k)
    end do
  end do

  return
end subroutine amudia
subroutine diamua ( nrow, job, a, ja, ia, diag, b, jb, ib )

!*****************************************************************************80
!
!! DIAMUA performs the matrix by matrix product B = Diag * A.
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    The algorithm is in-place; that is, B can take the place of A.
!    in this case use job=0.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) JOB, indicates the job to be done.
!    0, means get array B only;
!    1, means get B, and the integer ( kind = 4 ) arrays IB and JB.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(N), a diagonal matrix stored as a vector.
!
!    Output, real B(*), integer ( kind = 4 ) JB(*), 
!    integer ( kind = 4 ) IB(NROW+1), the resulting 
!    matrix B in compressed sparse row sparse format.
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) diag(nrow)
  integer ( kind = 4 ) ia(nrow+1)
  integer ( kind = 4 ) ib(nrow+1)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jb(*)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) scal

  do ii = 1, nrow
!
!  Normalize each row.
!
     k1 = ia(ii)
     k2 = ia(ii+1) - 1
     scal = diag(ii)
     b(k1:k2) = a(k1:k2) * scal

  end do

  if ( job == 0 ) then
    return
  end if

  ib(1) = ia(1)

  do ii = 1, nrow
    ib(ii) = ia(ii)
    do k = ia(ii), ia(ii+1)-1
      jb(k) = ja(k)
    end do
  end do

  return
end subroutine diamua

END MODULE Procedures