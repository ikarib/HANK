SUBROUTINE rtflsp(func,lx1,lx2,xacc,ftol,iflag,maxit)
!modified by greg to change initial region to allow for it to be outside range
!lx1 should evaluate to the low value of func

USE Globals

IMPLICIT NONE

REAL(8), INTENT(IN) :: lx1,lx2,xacc,ftol
INTEGER, INTENT(IN) :: maxit
INTEGER, INTENT(OUT) :: iflag
REAL(8) :: lrtflsp,x1,x2
INTEGER :: ieval,j
REAL(8) :: del,dx,f,fh,fl,xh,xl,ltemp
REAL(8), EXTERNAL	:: func

iflag = 0
ieval = 0

x1 = lx1
x2 = lx2

fl=func(x1)
ieval = ieval+1
if (ieval>maxit) then
	iflag = 1
	return
endif
if (abs(fl)<ftol) return
if (fl<0.0) goto 10
if (fl>0.0) goto 20


10 fh=func(x2)
ieval = ieval+1
if (ieval>maxit) then
	iflag = 1
	return
endif
if (abs(fh)<ftol) return
if (fh>0.0) goto 40
if (fh<0.0) then
   goto 30
endif

20 fh = fl
x2  = x1
x1 = 0.8*x1
fl = func(x1)
ieval = ieval+1
if (ieval>maxit) then
	iflag = 1
	return
endif
if (abs(fl)<ftol) return
if (fl<0.0) goto 40
if (fl>0.0) goto 20

30 fl = fh
x1 = x2
x2  = x2*1.2
fh = func(x2)
ieval = ieval+1
if (ieval>maxit) then
	iflag = 1
	return
endif
if (abs(fh)<ftol) return
if (fh>0.0) go to 40
if (fh<0.0) go to 30

if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) then
	write(*,*)'rtflsp: root must be bracketed between arguments'
	return
end if

40 if (fl < 0.0) then
	xl=x1
	xh=x2
else
	xl=x2
	xh=x1
	ltemp = fl
	fl = fh
	fh = ltemp
end if

dx=xh-xl
do j=1,maxit-ieval
	
	lrtflsp=xl+dx*fl/(fl-fh)
	f=func(lrtflsp)
	ieval = ieval+1
	
	if (f < 0.0) then
		del=xl-lrtflsp
		xl=lrtflsp
		fl=f
	else
		del=xh-lrtflsp
		xh=lrtflsp
		fh=f
	end if
	dx=xh-xl
	if (abs(del) < xacc .or. abs(f)<ftol) return
end do

! write(*,*)('rtflsp exceed maximum iterations')
iflag = 1

END SUBROUTINE rtflsp
