MODULE Globals
USE Parameters

IMPLICIT NONE

character(len=100)	OutputDir
integer				:: Display,RestrictJumpDomain,SaveSimulations,DriftPointsVersion,DefaultY1GridPar1,DefaultY2GridPar1,AddPointsCloseZeroY2
integer 			:: MatchVarLogY,MatchVarD1LogY,MatchSkewD1LogY,MatchKurtD1LogY,MatchVarD5LogY,MatchSkewD5LogY,MatchKurtD5LogY,MatchFracD1Less5,MatchFracD1Less10,MatchFracD1Less20,MatchFracD1Less50
real(8) 			:: TargetVarLogY,TargetVarD1LogY,TargetSkewD1LogY,TargetKurtD1LogY,TargetVarD5LogY,TargetSkewD5LogY,TargetKurtD5LogY,TargetFracD1Less5,TargetFracD1Less10,TargetFracD1Less20,TargetFracD1Less50
integer 			:: EstimateY1width,EstimateY2width,EstimateBeta,EstimateGamma,EstimateSigma,EstimateLambda1,EstimateZeta1P,EstimateZeta1N,EstimateLambda2,EstimateZeta2P,EstimateZeta2N
integer 			:: EstimateY1GridPar,EstimateY2GridPar,EstimateRho1,EstimateRho2,EstimateSigma1,EstimateSigma2

real(8), dimension(ngpy1) 		:: y1grid,y1dist
real(8), dimension(ngpy1,ngpy1)  	:: y1trans,y1markov,y1trans_qu
real(8), dimension(ngpy2) 		:: y2grid,y2dist
real(8), dimension(ngpy2,ngpy2)  	:: y2trans,y2markov,y2trans_qu

real(8), dimension(ngpy1*ngpy2) 		:: ygrid_combined,ydist_combined
real(8), dimension(ngpy1*ngpy2,ngpy1*ngpy2)  	:: ymarkov_combined,ytrans_qu_combined


real(8), dimension(nsim,Tsim)	:: y1rand,y2rand,ysim,ylevsim
integer, dimension(nsim,Tsim)	:: ysimI,y1simI,y2simI
real(8), dimension(nsim,5)		:: yannsim,yannlevsim

real(8)		:: y1width,y2width,y1gridpar,y2gridpar,beta1,beta2,sigma1,sigma2,lambda1,lambda2,rho1,rho2
real(8)		:: y1widthguess,y2widthguess,y1gridparguess,y2gridparguess
real(8)		:: y1widthmin,y1widthmax,y2widthmin,y2widthmax,y1gridmin,y2gridmin
real(8) 	:: deltaforapprox
real(8), dimension(:),allocatable	:: paramguess,paramout,paramscale

integer		:: nparam,nmoments,ndfls,iprint,maxfun,objeval
real(8) 	:: rhobeg,rhoend

real(8)		:: muy,mu2y,mu3y,mu4y,gam3y,gam4y
real(8)		:: muylev,mu2ylev,mu3ylev,mu4ylev,gam3ylev,gam4ylev
real(8)		:: mudy1,mu2dy1,mu3dy1,mu4dy1,gam3dy1,gam4dy1
real(8)		:: mudy5,mu2dy5,mu3dy5,mu4dy5,gam3dy5,gam4dy5
real(8)		:: fracdy1less5,fracdy1less10,fracdy1less20,fracdy1less50


END MODULE Globals