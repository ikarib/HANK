% INTEGER            :: in,it,it1,itN
% REAL(8)            :: lssvar1,lssvar2

lssvar = sigma.^2 ./ ((1-rho.^2)+2*delta./lambda);

if rho(1) == 1
    lssvar(1) = sigma(1)^2 / (1-0.99^2);
end

% load('../earnings_estimation_output/ysim.txt'); ysim_=permute(reshape(ysim,nsim,Tsim,2),[1 3 2]); clear ysim

% draw initial from normal distribution with same mean and variance
if UseNormalDist
    ysim = sqrt(lssvar).*yrand(:,:,1);
else
    ysim = yrand(:,:,1)./zetaP;
end

% simulate income path in dt increments
ylevsim = nan(nsim,Tsim);
ylevsim(:,1) = exp(sum(ysim,2));
for it = 1:Tsim-1
    yjumpI = yjumprand(:,:,it) > 1-dt*lambda;
    if UseNormalDist
        ysim1 = rho.*ysim + sigma.*yrand(:,:,it+1);
    else
        ysim1 = rho.*ysim + yrand(:,:,it+1)./zetaP;
    end
    ysim(yjumpI) = ysim1(yjumpI);
    if AdditiveDrift
        ysim1 = ysim - dt*delta.*(2*(ysim>0)-1);
    else
        ysim1 = (1-dt*delta).*ysim;
    end
    ysim(~yjumpI) = ysim1(~yjumpI);
    ylevsim(:,it+1) = exp(sum(ysim,2));
end

% aggregate to annual income
Yper = floor(4/dt);
yannlevsim = squeeze(sum(reshape(ylevsim(:,Tburn+1:Tburn+Yper*Tann),nsim,Yper,Tann),2));
% max(max(abs(log(yannlevsim)-load('../earnings_estimation_output/yannsim.txt'))))