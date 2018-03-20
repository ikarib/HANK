MODULE Procedures

USE Parameters
USE Globals
USE mUMFPACK

IMPLICIT NONE
CONTAINS

!-----------------------------------------------

REAL(8) FUNCTION utilfn(lx)
REAL(8), INTENT(in)	::  lx                    

IF(gam==1.0) THEN
	utilfn = prefshock*log(lx)
ELSE
	utilfn = prefshock*(lx**(1.0-gam)-1.0)/(1.0-gam)
! 	utilfn = (lx**(1.0-gam))/(1.0-gam)
END IF

END FUNCTION utilfn

!-----------------------------------------------

REAL(8) FUNCTION utilfninv(lu)
REAL(8), INTENT(in)	::  lu                 

IF(gam==1.0) THEN
	utilfninv = exp(lu/prefshock)
ELSE
	utilfninv = ((1.0-gam)*(lu/prefshock)+1.0)**(1.0/(1.0-gam))
END IF

END FUNCTION utilfninv

!-----------------------------------------------

REAL(8) FUNCTION utilfn1(lx)
REAL(8), INTENT(in)	::  lx                    
	
utilfn1 = prefshock*lx**(-gam)

END FUNCTION utilfn1

!-----------------------------------------------

REAL(8) FUNCTION utilfn1inv(lu)
REAL(8), INTENT(in)	::  lu                
	
utilfn1inv = (lu/prefshock)**(-1.0/gam)

END FUNCTION utilfn1inv


!-----------------------------------------------

REAL(8) FUNCTION adjcostfn(ld,la)
REAL(8), INTENT(in)	::  ld,la
REAL(8)				::  lx

IF(ld==0.0) THEN
	adjcostfn = 0.0
	RETURN
END IF

lx = ld/max(kappa3,la)

IF(lx == 0.0_8) THEN
	adjcostfn = 0.0
ELSE IF (lx>0.0_8) THEN
	adjcostfn = kappa0_d*abs(lx) + (abs(lx)**(1.0+kappa2_d)) *(kappa1_d**(-kappa2_d)) / (1.0+kappa2_d)
	adjcostfn = adjcostfn*max(kappa3,la)
ELSE IF (lx<0.0_8) THEN
	adjcostfn = kappa0_w*abs(lx) + (abs(lx)**(1.0+kappa2_w)) *(kappa1_w**(-kappa2_w))  / (1.0+kappa2_w)
	adjcostfn = adjcostfn*max(kappa3,la)
END IF

END FUNCTION adjcostfn

!-----------------------------------------------

REAL(8) FUNCTION adjcostfn1(ld,la)
REAL(8), INTENT(in)	::  ld,la
REAL(8)				::  lx

lx = ld/max(kappa3,la)

IF(lx>0.0) adjcostfn1 = kappa0_d + ((lx/kappa1_d)**kappa2_d) 
IF(lx<0.0) adjcostfn1 = -kappa0_w - ((-lx/kappa1_w)**kappa2_w) 

END FUNCTION adjcostfn1

!-----------------------------------------------

REAL(8) FUNCTION adjcostfn1inv(lchi,la)
REAL(8), INTENT(in)	::  lchi,la

IF(lchi>=-kappa0_w .and. lchi<=kappa0_d) adjcostfn1inv = 0.0
IF(lchi>kappa0_d) adjcostfn1inv = kappa1_d * (lchi-kappa0_d)**(1.0/kappa2_d)
IF(lchi<-kappa0_w) adjcostfn1inv = -kappa1_w * (-lchi-kappa0_w)**(1.0/kappa2_w)

adjcostfn1inv = adjcostfn1inv*max(kappa3,la)


IF(adjcostfn1inv>dmax)THEN
	IF(Display>=2)write(*,*),'high optimal d: ',adjcostfn1inv
	adjcostfn1inv = dmax
END IF	
END FUNCTION adjcostfn1inv



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
	
! 	yi(i) = yL(i) + (xi(i)-xL(i))*((yH(i)-yL(i))/(xH(i)-xL(i)))
	yi(i) = yL(i) + ((xi(i)-xL(i))/(xH(i)-xL(i)))*(yH(i)-yL(i))

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

IF(abs(xL-xH)<1.0e-12) THEN
	yi = 0.5*(yL+yH)
ELSE
	yi  = yL  + ((xi -xL )/(xH -xL ))*(yH -yL )
END IF

END SUBROUTINE LinInterp1

!-----------------------------------------------

SUBROUTINE LinInterpCumSum1 (n,x,y,xi,yi)
!this does linear interpolation of (x,cumsum(y)) at points only point,xi
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
yL  = SUM(y(1:locL))
yH  = yL + y(locL +1)


IF(abs(xL-xH)<1.0e-12) THEN
	yi = 0.5*(yL+yH)
ELSE
	yi  = yL  + ((xi -xL )/(xH -xL ))*(yH -yL )
END IF

END SUBROUTINE LinInterpCumSum1

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
! lstring = '('//trim(lstring) // 'E16.8)'
lstring = '('//trim(lstring) // 'E17.8E3)'
DO i1=1,n1
	WRITE(f,lstring) (mat(i1,:))
END DO

CLOSE(f)

END SUBROUTINE WriteMatrixExpon

!----------------------------------------------
SUBROUTINE ConditionalExpectation (n,x,f,g,dx,xL,xH,y)
!finds y = E[g(x) | xL<x<=xH] using a trapezoidal rule at boundaries

IMPLICIT NONE
INTEGER, INTENT(IN)		:: n
REAL(8), INTENT(IN)	    :: x(:),f(:),g(:),dx(:),xL,xH
REAL(8), INTENT(OUT)	:: y
REAL(8), DIMENSION(n)		:: gfdx,fdx
REAL(8)		:: cumgfdxL,cumfdxL,cumgfdxH,cumfdxH


gfdx = g*f*dx
fdx = f*dx

IF(xL>x(1)) THEN
	CALL LinInterpCumSum1(n,x,gfdx,xL,cumgfdxL)
	CALL LinInterpCumSum1(n,x,fdx,xL,cumfdxL)
ELSE
	cumgfdxL = 0.0
	cumfdxL = 0.0
END IF

IF(xH<x(n)) THEN	
	CALL LinInterpCumSum1(n,x,gfdx,xH,cumgfdxH)
	CALL LinInterpCumSum1(n,x,fdx,xH,cumfdxH)
	
ELSE
	cumgfdxH = SUM(gfdx)
	cumfdxH = SUM(fdx)
END IF

y = (cumgfdxH-cumgfdxL)/(cumfdxH-cumfdxL)

END SUBROUTINE ConditionalExpectation


!----------------------------------------------
SUBROUTINE AllocateSolutionType(soln)

IMPLICIT NONE
TYPE(SolutionType), INTENT(INOUT) :: soln

ALLOCATE(soln%V(ngpa,ngpb,ngpy))
ALLOCATE(soln%c(ngpa,ngpb,ngpy))
ALLOCATE(soln%s(ngpa,ngpb,ngpy))
ALLOCATE(soln%h(ngpa,ngpb,ngpy))
ALLOCATE(soln%d(ngpa,ngpb,ngpy))
ALLOCATE(soln%u(ngpa,ngpb,ngpy))
ALLOCATE(soln%gjoint(ngpa,ngpb,ngpy))
ALLOCATE(soln%bdot(ngpa,ngpb,ngpy))
ALLOCATE(soln%gvec(naby))
ALLOCATE(soln%gamarg(ngpa,ngpy))
ALLOCATE(soln%gbmarg(ngpb,ngpy))
ALLOCATE(soln%gmat(nab,ngpy))
ALLOCATE(soln%A(ngpy))
ALLOCATE(soln%B(ngpy))
ALLOCATE(soln%AU(ngpy))
IF (ComputeDiscountedMPC==1) THEN
	ALLOCATE(soln%mpc(ngpa,ngpb,ngpy))
	ALLOCATE(soln%subeff1ass(ngpa,ngpb,ngpy))
	ALLOCATE(soln%subeff2ass(ngpa,ngpb,ngpy))
	ALLOCATE(soln%wealtheff1ass(ngpa,ngpb,ngpy))
	ALLOCATE(soln%wealtheff2ass(ngpa,ngpb,ngpy))
END IF

END SUBROUTINE

!----------------------------------------------
SUBROUTINE AllocateCumulativePolicyType(cum)

IMPLICIT NONE
TYPE(CumulativePolicyType), INTENT(INOUT) :: cum

ALLOCATE(cum%ccum1(ngpa,ngpb,ngpy))
ALLOCATE(cum%ccum2(ngpa,ngpb,ngpy))
ALLOCATE(cum%ccum4(ngpa,ngpb,ngpy))
ALLOCATE(cum%dcum1(ngpa,ngpb,ngpy))
ALLOCATE(cum%dcum2(ngpa,ngpb,ngpy))
ALLOCATE(cum%dcum4(ngpa,ngpb,ngpy))

END SUBROUTINE

!----------------------------------------------

SUBROUTINE LinearUpdate(n,xn,fn,xn1,fn1,xn2)

IMPLICIT NONE

INTEGER, INTENT(IN)					:: n
REAL(8), INTENT(IN), DIMENSION(n)	:: xn,fn,xn1,fn1
REAL(8), INTENT(OUT), DIMENSION(n)	:: xn2
REAL(8), DIMENSION(n) 				:: lm
INTEGER								:: i

DO i = 1,n
	IF(abs(xn1(i)-xn(i))<1.0e-8) THEN
		xn2(i) = 0.5*(xn1(i)+xn(i))
	ELSE	
		lm(i) = (fn1(i)-fn(i)) / (xn1(i)-xn(i))
		IF(lm(i)>0.0) THEN
			write (*,*) 'warning, positive slope'
		END IF
		xn2(i) = (fn(i)-lm(i)*xn(i)) / (1.0-lm(i))		
	END IF
END DO


END SUBROUTINE


!----------------------------------------------

SUBROUTINE LogLinearGrowthUpdate(n,lstep,xn,fn,xn1,fn1,xn2)

IMPLICIT NONE

INTEGER, INTENT(IN)					:: n
REAL(8), INTENT(IN) 				:: lstep
REAL(8), INTENT(IN), DIMENSION(n)	:: xn,fn,xn1,fn1
REAL(8), INTENT(OUT), DIMENSION(n)	:: xn2
REAL(8), DIMENSION(n) 				:: lm,dxn2
REAL(8), DIMENSION(n)				:: dxn,dfn,dxn1,dfn1
INTEGER								:: i

dxn(1) = log(xn(1))
dxn(2:n) = log(xn(2:n)/xn(1:n-1))

dfn(1) = log(fn(1))
dfn(2:n) = log(fn(2:n)/fn(1:n-1))

dxn1(1) = log(xn1(1))
dxn1(2:n) = log(xn1(2:n)/xn1(1:n-1))

dfn1(1) = log(fn1(1))
dfn1(2:n) = log(fn1(2:n)/fn1(1:n-1))


DO i = 1,n
	IF(abs(dxn1(i)-dxn(i))<1.0e-8) THEN
		dxn2(i) = 0.5*(dxn1(i)+dxn(i))
	ELSE	
		lm(i) = (dfn1(i)-dfn(i)) / (dxn1(i)-dxn(i))
! 		IF(lm(i)>0.0) THEN
! 			write (*,*) 'warning, positive slope'
! 		END IF
		dxn2(i) = (dfn(i)-lm(i)*dxn(i)) / (1.0-lm(i))		
	END IF
END DO

xn2(1) = exp(dxn2(1))
DO i = 2,n
	xn2(i) = xn2(i-1)*exp(dxn2(i))
END DO

xn2 = lstep*xn2 + (1.0-lstep)*xn1

END SUBROUTINE

!----------------------------------------------

SUBROUTINE PartialUpdate(n,lstep,xn,fn,xn1)

IMPLICIT NONE

INTEGER, INTENT(IN)					:: n
REAL(8), INTENT(IN) 				:: lstep
REAL(8), INTENT(IN), DIMENSION(n)	:: xn,fn
REAL(8), INTENT(OUT), DIMENSION(n)	:: xn1

xn1 = lstep*fn + (1.0-lstep)*xn

END SUBROUTINE

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

SUBROUTINE WorldBondFunction(r,b,bss,rss,e)

!b = f(r)
IMPLICIT NONE

REAL(8), INTENT(IN)					:: r,bss,rss,e
REAL(8), INTENT(OUT)				:: b

b = bss*exp(e*(r-rss)) 

END SUBROUTINE

!----------------------------------------------
SUBROUTINE WorldBondInverse(b,r,bss,rss,e)

!r = f^-1(b)
IMPLICIT NONE

REAL(8), INTENT(IN)					:: b,bss,rss,e
REAL(8), INTENT(OUT)				:: r

r = rss + (1.0/e)*log(b/bss)

END SUBROUTINE

!----------------------------------------------

SUBROUTINE WorldBondFunction2(r,b,bss,rss,e)

!b = f(r)
IMPLICIT NONE

REAL(8), INTENT(IN)					:: r,bss,rss,e
REAL(8), INTENT(OUT)				:: b

b = bss + e*(r-rss)

END SUBROUTINE

!----------------------------------------------
SUBROUTINE WorldBondInverse2(b,r,bss,rss,e)

!r = f^-1(b)
IMPLICIT NONE

REAL(8), INTENT(IN)		:: b,bss,rss,e
REAL(8), INTENT(OUT)	:: r

r = rss + (1.0/e)*(b - bss)

END SUBROUTINE

!----------------------------------------------
SUBROUTINE CUMSUM(x,y)

REAL(8), INTENT(in)		::  x(:)
REAL(8), DIMENSION(size(x)), INTENT(OUT)	::  y
INTEGER 			:: i,n

n = size(x)
y(1) = x(1) 
DO i =2,n
	y(i) = y(i-1)+ x(i)
END DO

END SUBROUTINE CUMSUM


!----------------------------------------------
SUBROUTINE AdjustDistProportionately(lx,ldelta,lf,ladj,lg)

!requires lower bound of x to be zero.

IMPLICIT NONE

REAL(8), INTENT(in)		::  lx(:),ldelta(:),lf(:),ladj
REAL(8), DIMENSION(size(lx)), INTENT(OUT)	::  lg
INTEGER 	:: n,ix
REAL(8), DIMENSION(size(lx))	:: lf1,lgcum,lg1
REAL(8)		:: lg1mean,lf1mean

n = SIZE(lf)

!make it a valid distribution
lf1 = lf/SUM(lf*ldelta)
! lf1 = lf

!interpolate cumulative distribution
DO ix = 1,n
	CALL LinInterpCumSum1 (n,lx,lf1*ldelta,lx(ix)/ladj,lgcum(ix))
END DO	

!check upper bound
lgcum(n)=1.0
lgcum = MERGE(1.0_8,lgcum,lgcum>1.0_8)

!get marginal
lg1(1) = lgcum(1)/ldelta(1)
DO ix = 2,n
	lg1(ix) = (lgcum(ix)-lgcum(ix-1)) / ldelta(ix)
END DO	
lg1 = MERGE(0.0_8,lg1,abs(lg1)<1.0e-12_8)

!compute means
lg1mean = SUM(lx*ldelta*lg1)
lf1mean = SUM(lx*ldelta*lf1)

!scale to make mean correct
lg = ladj*lf1mean*lg1/lg1mean

!adjust mass at zero so its a valid dist
lg(1) = (1.0-SUM(lg(2:n)*ldelta(2:n))) /ldelta(1)

!adjust it back to have the same sum as on input
lg = lg*SUM(lf*ldelta)

END SUBROUTINE

!----------------------------------------------
SUBROUTINE InvertMatrix(matrix, inverse, n, errorflag)
	IMPLICIT NONE
	!Declarations
	INTEGER, INTENT(IN)  :: n
	INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	REAL(8), INTENT(IN)	:: matrix(n,n)  !Input matrix
	REAL(8), INTENT(OUT)	:: inverse(n,n) !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k
	REAL(8) :: m
	REAL(8), DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
END SUBROUTINE InvertMatrix

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
subroutine amux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
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
!    Input, real X(*), and array of length equal to the column dimension 
!    of A.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(*)
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(n)

  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
    do k = ia(i), ia(i+1)-1
      t = t + a(k) * x(ja(k))
    end do

    y(i) = t

  end do

  return
end subroutine amux

END MODULE Procedures