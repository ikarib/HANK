yannsim = log(yannlevsim(1:nsim,:));

% central moments: logs
muy = sum(yannsim(:,1))/nsim;
mu2y = sum((yannsim(:,1)-muy).^2)/nsim;
mu3y = sum((yannsim(:,1)-muy).^3)/nsim;
mu4y = sum((yannsim(:,1)-muy).^4)/nsim;

% standardised moments: logs
if mu2y>0
    gam3y = mu3y/mu2y^1.5;
    gam4y = mu4y/mu2y^2;
else
    gam3y = 0;
    gam4y = 0;
end

% central moments: levels
muylev = sum(yannlevsim(:,1))/nsim;
mu2ylev = sum((yannlevsim(:,1)-muylev).^2)/nsim;
mu3ylev = sum((yannlevsim(:,1)-muylev).^3)/nsim;
mu4ylev = sum((yannlevsim(:,1)-muylev).^4)/nsim;

% standardised moments: levels
if mu2ylev>0
    gam3ylev = mu3ylev/mu2ylev^1.5;
    gam4ylev = mu4ylev/mu2ylev^2;
else
    gam3ylev = 0;
    gam4ylev = 0;
end

% central moments: 1 year log changes
mudy1 = sum(yannsim(:,2)-yannsim(:,1))/nsim;
mu2dy1 = sum((yannsim(:,2)-yannsim(:,1)-mudy1).^2)/nsim;
mu3dy1 = sum((yannsim(:,2)-yannsim(:,1)-mudy1).^3)/nsim;
mu4dy1 = sum((yannsim(:,2)-yannsim(:,1)-mudy1).^4)/nsim;

% standardised moments: 1 year log changes
if mu2dy1>0
    gam3dy1 = mu3dy1/mu2dy1^1.5;
    gam4dy1 = mu4dy1/mu2dy1^2;
else
    gam3dy1 = 0;
    gam4dy1 = 0;
end

% central moments: 5 year log changes
mudy5 = sum(yannsim(:,5)-yannsim(:,1))/nsim;
mu2dy5 = sum((yannsim(:,5)-yannsim(:,1)-mudy5).^2)/nsim;
mu3dy5 = sum((yannsim(:,5)-yannsim(:,1)-mudy5).^3)/nsim;
mu4dy5 = sum((yannsim(:,5)-yannsim(:,1)-mudy5).^4)/nsim;

% standardised moments: 5 year log changes
if mu2dy5>0
    gam3dy5 = mu3dy5/mu2dy5^1.5;
    gam4dy5 = mu4dy5/mu2dy5^2;
else
    gam3dy5 = 0;
    gam4dy5 = 0;
end

% fraction 1 year log changes in ranges
fracdy1less5 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.05)/nsim;
fracdy1less10 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.1)/nsim;
fracdy1less20 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.2)/nsim;
fracdy1less50 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.5)/nsim;
