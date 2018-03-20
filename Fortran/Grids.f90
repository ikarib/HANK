SUBROUTINE Grids

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER			:: ia,ib,iy,iaby,iab
REAL(8)         :: lmean
REAL(8)         :: leye(ngpy,ngpy)

REAL(8), EXTERNAL :: golden

!identity matrix
leye = 0.0; DO iy = 1,ngpy
	leye(iy,iy) = 1.0
END DO

!productivity process
IF(ngpy==1) THEN
	logygrid = 0.0
	ygrid = 1.0
	ymarkov = 0.0
	ytrans = 1.0
	ydist = 1.0

ELSE IF (ReadEarningsProcess==1) THEN
	OPEN(1, FILE = trim(EarningsProcessDir) // '/ygrid_combined.txt');READ(1,*) logygrid;CLOSE(1)
	OPEN(1, FILE = trim(EarningsProcessDir) // '/ydist_combined.txt');READ(1,*) ydist;CLOSE(1)
	OPEN(1, FILE = trim(EarningsProcessDir) // '/ymarkov_combined.txt');READ(1,*) ymarkov;CLOSE(1)
	IF (AdjustProdGridFrisch==1) logygrid = logygrid/ (1.0+adjfricshgridfrac*frisch)
	ygrid = exp(logygrid)
	ymarkov = TRANSPOSE(ymarkov) !since fortran reads in column major order	
	!fix up rounding in markov matrix
	DO iy = 1,ngpy
		ymarkov(iy,iy) = ymarkov(iy,iy) - SUM(ymarkov(iy,:))
	END DO
	
ELSE IF (TwoPointWageProcess==1) THEN
	ygrid(1) = 0.8
	ygrid(2) = 1.2
	logygrid = log(ygrid)
	
	ytrans(1,1) = 1.0-0.06667
	ytrans(1,2) = 0.06667
	ytrans(2,1) = 0.06667
	ytrans(2,2) = 1.0-0.06667
	ydist(1) = 0.5
	ydist(2) = 0.5
	ymarkov = (ytrans-leye)/1.0 !assumes ytrans is quarterly	
END IF

ymarkovdiag = 0.0
DO iy = 1,ngpy
	ymarkovdiag(iy,iy) = ymarkov(iy,iy)
END DO
ymarkovoff = ymarkov-ymarkovdiag

!adjust mean productivity
lmean = DOT_PRODUCT(ydist,ygrid)
ygrid = meanlabeff*ygrid/lmean

CALL PowerSpacedGrid (ngpa,agridparam,0.0_8,amax,agrid)

!with low gridparam points get bunched close to zero, so evenly space first 8 points;
IF(ngpa>10) THEN
	DO ia = 1,10-1
		agrid(ia) = (ia-1)*agrid(10)/(10.0-1.0)
	END DO
END IF
dagrid = agrid(2:ngpa)-agrid(1:ngpa-1)

!liquid grid
IF(Borrowing==0) THEN
	CALL PowerSpacedGrid (ngpb,bgridparam,0.0_8,bmax,bgrid)
ELSE IF(Borrowing==1) THEN
	CALL PowerSpacedGrid (ngpbPOS,bgridparam,0.0_8,bmax,bgrid(ngpbNEG+1:ngpb))
	nbl = -lumptransfer/(rborr+PerfectAnnuityMarkets*deathrate)
 	abl = max(nbl + cmin,blim)
	IF (Display>=2) write(*,*) 'natural borrowing limit = ',nbl
	IF (Display>=2) write(*,*) 'actual borrowing limit = ',abl
	CALL PowerSpacedGrid (ngpbNEG/2+1,bgridparamNEG,abl,(abl+bgrid(ngpbNEG+1))/2.0,bgrid(1:ngpbNEG/2+1))
	DO ib = ngpbNEG/2+2,ngpbNEG
		bgrid(ib) = bgrid(ngpbNEG+1) -(bgrid(ngpbNEG+2-ib)-bgrid(1))
	END DO
END IF
dbgrid = bgrid(2:ngpb)-bgrid(1:ngpb-1)


adelta(1) = 0.5*dagrid(1)
adelta(2:ngpa-1) = 0.5*(dagrid(1:ngpa-2)+dagrid(2:ngpa-1))
adelta(ngpa) = 0.5*dagrid(ngpa-1)

bdelta(1) = 0.5*dbgrid(1)
bdelta(2:ngpb-1) = 0.5*(dbgrid(1:ngpb-2)+dbgrid(2:ngpb-1))
bdelta(ngpb) = 0.5*dbgrid(ngpb-1)


!combined grids for openmp loops
iab = 1
DO ia = 1,ngpa
DO ib = 1,ngpb
	afromab(iab) = ia
	bfromab(iab) = ib
	abfromab(ia,ib) = iab
	abdelta(iab) = adelta(ia)*bdelta(ib)
	iab = iab+1
END DO
END DO

iaby = 1
DO ia = 1,ngpa
DO ib = 1,ngpb
DO iy = 1,ngpy
	afromaby(iaby) = ia
	bfromaby(iaby) = ib
	yfromaby(iaby) = iy
	abyfromaby(ia,ib,iy) = iaby
	abydelta(iaby) = adelta(ia)*bdelta(ib)
	abfromaby(iaby) = abfromab(ia,ib)
	iaby = iaby+1
END DO
END DO
END DO


END SUBROUTINE Grids

