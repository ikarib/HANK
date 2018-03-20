SUBROUTINE rtsec(func,lx1,lx2,lfacc,lxsol,iflag)

IMPLICIT NONE

REAL(8), INTENT(IN) :: lx1,lx2,lfacc
REAL(8), INTENT(OUT) :: lxsol
INTEGER, INTENT(OUT) :: iflag
REAL(8) :: lrtsec,x1,x2
INTEGER, PARAMETER :: MAXIT=20
INTEGER :: j
REAL(8) :: dx,f,fl,xl,ltemp
REAL(8), EXTERNAL	:: func

iflag = 0

x1 = lx1
x2 = lx2
! write(*,*)'-----------'
fl=func(x1)
! write(*,*)x1,fl

if (abs(fl)<lfacc) then
	lxsol = x1
	return
end if

f=func(x2)
! write(*,*)x2,f

if (abs(f)<lfacc) then
	lxsol = x2
	return
end if


if (abs(fl) < abs(f)) then
	lrtsec=x1
	xl=x2
	ltemp = fl
	fl = f
	f = ltemp
else
	xl=x1
	lrtsec=x2
end if

do j=1,MAXIT
	dx=(xl-lrtsec)*f/(f-fl)
	xl=lrtsec
	fl=f
	lrtsec=lrtsec+dx
	f=func(lrtsec)
! 	write(*,*)lrtsec,f
    if (abs(f)<lfacc) then
		lxsol = lrtsec
		return
	end if
end do
	
! write(*,*)('rtsec: exceed maximum iterations')
lxsol = lrtsec
iflag = 1

END SUBROUTINE rtsec