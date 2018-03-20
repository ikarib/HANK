MODULE Globals
USE Parameters

IMPLICIT NONE

character(len=100)	OutputDir
integer				:: Display,SaveSimulations,RevertToMedianWithSkew,UseDoubleParetoDist,UseNormalDist,Include2ndProcess,AdditiveDrift
integer 			:: MatchVarLogY,MatchVarD1LogY,MatchSkewD1LogY,MatchKurtD1LogY,MatchVarD5LogY,MatchSkewD5LogY,MatchKurtD5LogY,MatchFracD1Less5,MatchFracD1Less10,MatchFracD1Less20,MatchFracD1Less50
real(8) 			:: TargetVarLogY,TargetVarD1LogY,TargetSkewD1LogY,TargetKurtD1LogY,TargetVarD5LogY,TargetSkewD5LogY,TargetKurtD5LogY,TargetFracD1Less5,TargetFracD1Less10,TargetFracD1Less20,TargetFracD1Less50
integer 			:: EstimateLambda1,EstimateZeta1P,EstimateZeta1N,EstimateLambda2,EstimateZeta2P,EstimateZeta2N,EstimateRho1,EstimateRho2
integer 			:: EstimateSigma1,EstimateSigma2,EstimateDelta1,EstimateDelta2


real(8), dimension(nsim,Tsim)	:: y1jumprand,y2jumprand,y1rand,y2rand,ysim,ylevsim,y1sim,y2sim
integer, dimension(nsim,Tsim)	:: y1jumpI,y2jumpI
real(8), dimension(nsim,5)		:: yannsim,yannlevsim

real(8)		:: lambda1,zeta1P,zeta1N,lambda2,zeta2P,zeta2N,rho1,rho2,sigma1,sigma2,delta1,delta2
real(8)		:: lambda1guess,zeta1guess,lambda2guess,zeta2guess,rho1guess,rho2guess,sigma1guess,sigma2guess,delta1guess,delta2guess
real(8)		:: lambdamax,sigmamax,zetamax,zetamin,deltamax
real(8), dimension(:),allocatable	:: paramguess,paramout,paramscale

integer		:: nparam,nmoments,ndfls,iprint,maxfun,objeval
real(8) 	:: rhobeg,rhoend

real(8)		:: muy,mu2y,mu3y,mu4y,gam3y,gam4y
real(8)		:: muylev,mu2ylev,mu3ylev,mu4ylev,gam3ylev,gam4ylev
real(8)		:: mudy1,mu2dy1,mu3dy1,mu4dy1,gam3dy1,gam4dy1
real(8)		:: mudy5,mu2dy5,mu3dy5,mu4dy5,gam3dy5,gam4dy5
real(8)		:: fracdy1less5,fracdy1less10,fracdy1less20,fracdy1less50


END MODULE Globals