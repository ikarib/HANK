SUBROUTINE rtbis(func,x1,x2,xacc,ftol,xout)
IMPLICIT NONE
REAL(8), INTENT(IN) :: x1,x2,xacc,ftol
REAL(8), INTENT(OUT) :: xout
REAL(8) :: lrtbis
INTEGER, PARAMETER :: MAXIT=80
INTEGER :: j
REAL(8) :: dx,f,fmid,xmid
REAL(8), EXTERNAL	:: func
fmid=func(x2)
f=func(x1)
if (f*fmid >= 0.0) then
	if (abs(f)<ftol .or. abs(fmid)<ftol) then
		if (abs(f)<ftol .and. abs(fmid)<ftol) write(*,*) 'rtbis: both initial points have abs <ftol'
		if(abs(f)<=abs(fmid)) xout = x1
		if(abs(fmid)<abs(f)) xout = x2
! 		f = func(xout)	!so that last one called is correct to set globals.
		return
	end if
	write(*,*)'rtbis: root must be bracketed'
	return
end if
if (f < 0.0) then
	lrtbis=x1
	dx=x2-x1
else
	lrtbis=x2
	dx=x1-x2
end if
do j=1,MAXIT
	dx=dx*0.5_8
	xmid=lrtbis+dx
	fmid=func(xmid)
	if (fmid <= 0.0) lrtbis=xmid
	if (abs(dx) < xacc .or. abs(fmid) < ftol) then
		xout = lrtbis
		RETURN
	end if
end do
write(*,*)'rtbis: too many bisections'

END SUBROUTINE rtbis
