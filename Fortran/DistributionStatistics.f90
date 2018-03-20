SUBROUTINE DistributionStatistics

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER			:: ia,ib,iy,iaby,ip,iab
REAL(8), DIMENSION(naby) :: la,lb,lc,ld,lwage,lnetlabinc,lgrosslabinc,ladelta,lbdelta,lhours,ladjcost,lgrossinc,lnw,llabor,lgrossprofinc,lnetprofinc
REAL(8)	:: lbmargdist(ngpb),lamargdist(ngpa),lbmargcum(ngpb),lamargcum(ngpa),lpvec(11),linterp1,linterp2
INTEGER, DIMENSION(nab)		:: ordernw
REAL(8), DIMENSION(nab)		:: lnwmargdist,lnwmargcum,lnwgrid,lnwdelta,lnw_a,lnw_b,lnw_c,lnw_inc,lnw_h
REAL(8), DIMENSION(naby)	:: lcmargdist,lcmargcum,lcgrid,lcdelta,lmargdist
REAL(8), DIMENSION(naby)	:: lincgrid,lincmargdist,lincdelta,lincmargcum,linc_a,linc_b,linc_c,linc_nw,linc_h
INTEGER, DIMENSION(naby)	:: orderinc


!vectorize everything
!$OMP PARALLEL DO PRIVATE(ia,ib,iy,iab)
DO iaby = 1,naby
	ia = afromaby(iaby)
	ib = bfromaby(iaby)
	iy = yfromaby(iaby)	
	iab = abfromab(ia,ib)
	orderinc(iaby) = iaby
	
	la(iaby) = agrid(ia)
	lb(iaby) = bgrid(ib)
	lnw(iaby) = agrid(ia)+ bgrid(ib)
	lc(iaby) = c(ia,ib,iy)
	ld(iaby) = d(ia,ib,iy)
	ladjcost(iaby) = adjcostfn(ld(iaby),la(iaby))
	lwage(iaby) = ygrid(iy)*wage
	lhours(iaby) = h(ia,ib,iy)
	llabor(iaby) = h(ia,ib,iy)*ygrid(iy)
	lgrosslabinc(iaby) = lhours(iaby)*ygrid(iy)*wage
	lnetlabinc(iaby) = lhours(iaby)*ygrid(iy)*netwage + lumptransfer
	IF(DistributeProfitsInProportion==0) lgrossprofinc(iaby) = 0.0 
	IF(DistributeProfitsInProportion==1) lgrossprofinc(iaby) = (1.0-profdistfrac)*profit*ygrid(iy)/meanlabeff
	IF(TaxHHProfitIncome==0) lnetprofinc(iaby) = lgrossprofinc(iaby)
	IF(TaxHHProfitIncome==1) lnetprofinc(iaby) = (1.0-labtax)*lgrossprofinc(iaby)
	
	lgrossinc(iaby) = lgrosslabinc(iaby) + lgrossprofinc(iaby) + (rb+PerfectAnnuityMarkets*deathrate)*max(bgrid(ib),0.0) + (rborr+PerfectAnnuityMarkets*deathrate)*min(bgrid(ib),0.0) + (ra+PerfectAnnuityMarkets*deathrate)*agrid(ia)
	ladelta(iaby) = adelta(ia)
	lbdelta(iaby) = bdelta(ib)
	lmargdist(iaby) = gmat(iab,iy)*abydelta(iaby)
	
	IF(iy==1) lnwgrid(iab) = agrid(ia) + bgrid(ib)
	IF(iy==1) ordernw(iab) = iab
	
END DO
!$OMP END PARALLEL DO

!liquid wealth marginal dist
lbmargdist = SUM(gbmarg,DIM=2)*bdelta
CALL CUMSUM(lbmargdist,lbmargcum)

!illiquid wealth marginal dist
lamargdist = SUM(gamarg,DIM=2)*adelta
CALL CUMSUM(lamargdist,lamargcum)

!other marginal dist
IF ((iteratingtransition==.false.) .and. (calibrating==.false.) ) THEN

	!networth grid and marginal dist
	CALL isort2(nab,lnwgrid,ordernw)	
	lnwdelta = 0.0
	lnwdelta(1:nab-1) = 0.5*(lnwgrid(2:nab)-lnwgrid(1:nab-1))
	lnwdelta(2:nab) = lnwdelta(2:nab) + 0.5*(lnwgrid(2:nab)-lnwgrid(1:nab-1))
	
	!$OMP PARALLEL DO 	
	DO iab = 1,nab
		lnwmargdist(iab) = SUM(gmat(ordernw(iab),:))*abdelta(ordernw(iab)) / lnwdelta(iab)
		IF( SUM(lmargdist, MASK = abfromaby == ordernw(iab)) > 1.0e-12) THEN
			lnw_a(iab)  = SUM(la*lmargdist, MASK = abfromaby == ordernw(iab)) /  SUM(lmargdist, MASK = abfromaby == ordernw(iab))
			lnw_b(iab)  = SUM(lb*lmargdist, MASK = abfromaby == ordernw(iab)) /  SUM(lmargdist, MASK = abfromaby == ordernw(iab))
			lnw_c(iab)  = SUM(lc*lmargdist, MASK = abfromaby == ordernw(iab)) /  SUM(lmargdist, MASK = abfromaby == ordernw(iab))
			lnw_h(iab)  = SUM(lhours*lmargdist, MASK = abfromaby == ordernw(iab)) /  SUM(lmargdist, MASK = abfromaby == ordernw(iab))
			lnw_inc(iab)  = SUM(lgrossinc*lmargdist, MASK = abfromaby == ordernw(iab)) /  SUM(lmargdist, MASK = abfromaby == ordernw(iab))
		ELSE
			lnw_a(iab)  = 0.0
			lnw_b(iab)  = 0.0
			lnw_c(iab)  = 0.0
			lnw_h(iab)  = 0.0
			lnw_inc(iab)  = 0.0
		END IF
	END DO
	!$OMP END PARALLEL DO
	CALL CUMSUM(lnwmargdist*lnwdelta,lnwmargcum)

	!consumption
	lcgrid = lc
	lcmargdist = lmargdist
	CALL sort2(naby,lcgrid,lcmargdist)
	lcdelta = 0.0
	lcdelta(1:naby-1) = 0.5*(lcgrid(2:naby)-lcgrid(1:naby-1))
	lcdelta(2:naby) = lcdelta(2:naby) + 0.5*(lcgrid(2:naby)-lcgrid(1:naby-1))
	CALL CUMSUM(lcmargdist,	lcmargcum)

	!gross total income
	lincgrid = lgrossinc
	CALL isort2(naby,lincgrid,orderinc)
	lincdelta = 0.0
	lincdelta(1:naby-1) = 0.5*(lincgrid(2:naby)-lincgrid(1:naby-1))
	lincdelta(2:naby) = lincdelta(2:naby) + 0.5*(lincgrid(2:naby)-lincgrid(1:naby-1))

	lincmargdist = lmargdist(orderinc)/lincdelta
	linc_a = la(orderinc)
	linc_b = lb(orderinc)
	linc_c = lc(orderinc)
	linc_h = lhours(orderinc)
	linc_nw = lnw(orderinc)
	CALL CUMSUM(lincmargdist*lincdelta,lincmargcum)
	
END IF

!means
Ehours = SUM(lhours*gvec*ladelta*lbdelta)
Elabor = SUM(llabor*gvec*ladelta*lbdelta)
Ewage = SUM(lwage*gvec*ladelta*lbdelta)
Enetlabinc = SUM(lnetlabinc*gvec*ladelta*lbdelta)
Egrosslabinc = SUM(lgrosslabinc*gvec*ladelta*lbdelta)
Enetprofinc = SUM(lnetprofinc*gvec*ladelta*lbdelta)
Egrossprofinc = SUM(lgrossprofinc*gvec*ladelta*lbdelta)
Einc = SUM(lincgrid*lincmargdist*lincdelta)
Ea = SUM(la*gvec*ladelta*lbdelta)
Eb = SUM(lb*gvec*ladelta*lbdelta)
Ec = SUM(lc*gvec*ladelta*lbdelta)
Ed = SUM(ld*gvec*ladelta*lbdelta)
Eadjcost = SUM(ladjcost*gvec*ladelta*lbdelta)
EbP = SUM(lb*gvec*ladelta*lbdelta,MASK= lb>0)
EbN = SUM(lb*gvec*ladelta*lbdelta,MASK= lb<0)
Enw = SUM(lnwgrid*lnwmargdist*lnwdelta)

!frac at exactly zero
FRACa0 = lamargdist(1)
IF(Borrowing==1) FRACb0 = lbmargdist(ngpbNEG+1)
IF(Borrowing==0) FRACb0 = lbmargdist(1)
FRACnw0 = SUM(MERGE(1.0_8,0.0_8,abs(lnwgrid)<1.0e-8)*lnwmargdist*lnwdelta)
FRACb0a0 = SUM(MERGE(1.0_8,0.0_8,abs(la)<1.0e-8)*MERGE(1.0_8,0.0_8,abs(lb)<1.0e-8)*gvec*ladelta*lbdelta)
FRACb0aP = SUM(MERGE(1.0_8,0.0_8,la>1.0e-8)*MERGE(1.0_8,0.0_8,abs(lb)<1.0e-8)*gvec*ladelta*lbdelta)

!frac zero or close but greater than zero
CALL LinInterp1 (ngpb,bgrid,lbmargcum,defnbclose*Egrosslabinc,linterp2)
IF(Borrowing==1) FRACb0close = linterp2-lbmargcum(ngpbNEG)
IF(Borrowing==0) FRACb0close = linterp2
CALL LinInterp1 (ngpa,agrid,lamargcum,defnaclose*Egrosslabinc,FRACa0close)
IF ((iteratingtransition==.false.) .and. (calibrating==.false.)) THEN
	CALL LinInterp1 (nab,lnwgrid,lnwmargcum,(defnbclose+defnaclose)*Egrosslabinc,linterp2)
	IF(Borrowing==0) FRACnw0close = linterp2
	IF(Borrowing==1) CALL LinInterp1 (nab,lnwgrid,lnwmargcum,-1.0e-8_8,linterp1)
	IF(Borrowing==1) FRACnw0close = linterp2-linterp1
END IF
IF (iteratingtransition==.false.) THEN !since we don't store gabcum
	CALL BiLinInterp1 (ngpa,agrid,ngpb,bgrid,gabcum,defnaclose*Egrosslabinc,defnbclose*Egrosslabinc,linterp2)
	IF(Borrowing==0) FRACb0a0close = linterp2
	IF(Borrowing==1) THEN
		CALL LinInterp1 (ngpa,agrid,gabcum(:,ngpbNEG),defnaclose*Egrosslabinc,linterp1)
		FRACb0a0close = linterp2-linterp1
	END IF
END IF

FRACbN = SUM(MERGE(1.0_8,0.0_8,lb<-1.0e-12)*gvec*ladelta*lbdelta)

Ec_bN = SUM(lc*gvec*ladelta*lbdelta, MASK = lb<-1.0e-8)/FRACbN
Ec_b0close = SUM(lc*gvec*ladelta*lbdelta, MASK = lb>-1.0e-8 .and. lb<defnbclose*Egrosslabinc)/FRACb0close
Ec_b0far = SUM(lc*gvec*ladelta*lbdelta, MASK = lb>=defnbclose*Egrosslabinc)/(1.0-FRACb0close-FRACbN)

!percentiles: use cumulative marginal distributions
lpvec = (/0.01,0.02,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.98,0.99 /)

IF ((iteratingtransition==.false.) .and. (calibrating==.false. .or. MatchMedianLiq==1) .and. (ngpy>1 .or. deathrate>0.0)) THEN
	DO ip = 1,11
		
		!liquid wealth	
		IF (lbmargcum(1)>=lpvec(ip)) THEN
			PERCb(ip) = bgrid(1)
		ELSE	
			CALL LinInterp1 (ngpb,lbmargcum,bgrid,lpvec(ip),PERCb(ip))
		END IF
	END DO
END IF	

IF ((iteratingtransition==.false.) .and. (calibrating==.false. .or. MatchMedianIll==1 .or. MatchP75Ill==1) .and. (ngpy>1 .or. deathrate>0.0)) THEN
	DO ip = 1,11
	
		!iliquid wealth
		IF (lamargcum(1)>=lpvec(ip)) THEN
			PERCa(ip) = agrid(1)
		ELSE	
			CALL LinInterp1 (ngpa,lamargcum,agrid,lpvec(ip),PERCa(ip))
		END IF
	END DO
END IF	
	
IF ((iteratingtransition==.false.) .and. (calibrating==.false.) .and. (ngpy>1 .or. deathrate>0.0)) THEN
	DO ip = 1,11

		!net worth
		IF (lnwmargcum(1)>=lpvec(ip)) THEN
			PERCnw(ip) = lnwgrid(1)
		ELSE	
			CALL LinInterp1 (nab,lnwmargcum,lnwgrid,lpvec(ip),PERCnw(ip))
		END IF

		!consumption
		IF (lcmargcum(1)>=lpvec(ip)) THEN
			PERCc(ip) = lcgrid(1)
		ELSE	
			CALL LinInterp1 (naby,lcmargcum,lcgrid,lpvec(ip),PERCc(ip))
		END IF

		!gross income
		IF (lincmargcum(1)>=lpvec(ip)) THEN
			PERCinc(ip) = lincgrid(1)
		ELSE	
			CALL LinInterp1 (naby,lincmargcum,lincgrid,lpvec(ip),PERCinc(ip))
		END IF

	END DO

	!gini coefficient
	IF(Ea>0.0) GINIa = SUM(lamargcum*(1.0-lamargcum)*adelta) / Ea
	IF(Ea==0.0) GINIa = 0.0
	GINIb = SUM(lbmargcum*(1.0-lbmargcum)*bdelta) / Eb
	GINInw = SUM(lnwmargcum*(1.0-lnwmargcum)*lnwdelta) / Enw
	GINIc = SUM(lcmargcum*(1.0-lcmargcum)*lcdelta) / Ec
	GINIinc = SUM(lincmargcum*(1.0-lincmargcum)*lincdelta) / Einc

	!in transition use groupings based on steady sate	
	
	!statistics conditional on quartile of net worth distributions	
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_a,lnwdelta,lnwgrid(1),PERCnw(5),Ea_nwQ(1))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_a,lnwdelta,PERCnw(5),PERCnw(6),Ea_nwQ(2))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_a,lnwdelta,PERCnw(6),PERCnw(7),Ea_nwQ(3))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_a,lnwdelta,PERCnw(7),lnwgrid(nab),Ea_nwQ(4))

	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_b,lnwdelta,lnwgrid(1),PERCnw(5),Eb_nwQ(1))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_b,lnwdelta,PERCnw(5),PERCnw(6),Eb_nwQ(2))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_b,lnwdelta,PERCnw(6),PERCnw(7),Eb_nwQ(3))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_b,lnwdelta,PERCnw(7),lnwgrid(nab),Eb_nwQ(4))

	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_c,lnwdelta,lnwgrid(1),PERCnw(5),Ec_nwQ(1))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_c,lnwdelta,PERCnw(5),PERCnw(6),Ec_nwQ(2))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_c,lnwdelta,PERCnw(6),PERCnw(7),Ec_nwQ(3))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnw_c,lnwdelta,PERCnw(7),lnwgrid(nab),Ec_nwQ(4))

	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnwgrid,lnwdelta,lnwgrid(1),PERCnw(5),Einc_nwQ(1))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnwgrid,lnwdelta,PERCnw(5),PERCnw(6),Einc_nwQ(2))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnwgrid,lnwdelta,PERCnw(6),PERCnw(7),Einc_nwQ(3))
	CALL ConditionalExpectation (nab,lnwgrid,lnwmargdist,lnwgrid,lnwdelta,PERCnw(7),lnwgrid(nab),Einc_nwQ(4))
	
	!statistics conditional on quartile of total gross income distributions
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_a,lincdelta,lincgrid(1),PERCinc(5),Ea_incQ(1))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_a,lincdelta,PERCinc(5),PERCinc(6),Ea_incQ(2))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_a,lincdelta,PERCinc(6),PERCinc(7),Ea_incQ(3))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_a,lincdelta,PERCinc(7),lincgrid(naby),Ea_incQ(4))

	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_b,lincdelta,lincgrid(1),PERCinc(5),Eb_incQ(1))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_b,lincdelta,PERCinc(5),PERCinc(6),Eb_incQ(2))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_b,lincdelta,PERCinc(6),PERCinc(7),Eb_incQ(3))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_b,lincdelta,PERCinc(7),lincgrid(naby),Eb_incQ(4))

	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_c,lincdelta,lincgrid(1),PERCinc(5),Ec_incQ(1))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_c,lincdelta,PERCinc(5),PERCinc(6),Ec_incQ(2))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_c,lincdelta,PERCinc(6),PERCinc(7),Ec_incQ(3))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,linc_c,lincdelta,PERCinc(7),lincgrid(naby),Ec_incQ(4))

	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,lincgrid,lincdelta,lincgrid(1),PERCinc(5),Einc_incQ(1))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,lincgrid,lincdelta,PERCinc(5),PERCinc(6),Einc_incQ(2))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,lincgrid,lincdelta,PERCinc(6),PERCinc(7),Einc_incQ(3))
	CALL ConditionalExpectation (naby,lincgrid,lincmargdist,lincgrid,lincdelta,PERCinc(7),lincgrid(naby),Einc_incQ(4))

END IF

END SUBROUTINE DistributionStatistics