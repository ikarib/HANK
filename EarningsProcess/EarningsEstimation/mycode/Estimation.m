% Globals:
% Integers:
% nparam,nmoments,maxfun,objeval
% Real:
% lambda1,zeta1P,zeta1N,lambda2,zeta2P,zeta2N,rho1,rho2,sigma1,sigma2,delta1,delta2
% muy,mu2y,mu3y,mu4y,gam3y,gam4y
% muylev,mu2ylev,mu3ylev,mu4ylev,gam3ylev,gam4ylev
% mudy1,mu2dy1,mu3dy1,mu4dy1,gam3dy1,gam4dy1
% mudy5,mu2dy5,mu3dy5,mu4dy5,gam3dy5,gam4dy5
% fracdy1less5,fracdy1less10,fracdy1less20,fracdy1less50

% Arrays:
% yjumprand,yrand,ysim = nan(2,nsim,Tsim)
% ysim,ylevsim = nan(nsim,Tsim)
% yjumpI = false(nsim,Tsim,2)
% yannsim,yannlevsim = nan(nsim,5)

% local integer :: ip,j,npt
% arrays :: x(:),w(:)


% set up parameters and assign guesses
% nparam = numel(paramguess);

% set up moments
% nmoments = sum(Match);

% dfls estimation
% diary([OutputDir 'iterations.txt'])
% xmin = bounds(Estimate,[1 1]);
xmax = bounds(Estimate,[2 2]);
% x = guess(Estimate,:);
x = -log(xmax./guess(Estimate,:)-1);
x = lsqnonlin(@dfovec,x); % ,xmin,xmax
guess(Estimate,:) = xmax./(1+exp(-x));
% npt = 2*nparam+1;
% rhobeg = 5.0;
% rhoend = 1.0e-4;
% iprint = 4;
% maxfun = 500*(nparam+1);
% ALLOCATE(w((npt+11)*(npt+nparam) +nparam*(3*nparam+11)/2) )
% NEWUOA_H(nparam,npt,x,rhobeg,rhoend,iprint,maxfun,w,nmoments)
% options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
%     'MaxFunctionEvaluations',1500)
