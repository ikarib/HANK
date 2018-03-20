function lg=AdjustDistProportionately(lx,ldelta,lf,ladj)
% requires lower bound of x to be zero.
[n,m,k]=size(lf);
lf = lf(:,:);

% make it a valid distribution
lf1 = lf./(ldelta'*lf);

% interpolate cumulative distribution
lgcum = interp1(lx,cumsum(lf1.*ldelta),lx/ladj);

% check upper bound
lgcum(n,:) = 1;
lgcum(lgcum>1) = 1;

% get marginal
lg1 = [lgcum(1,:); diff(lgcum)]./ldelta;
lg1(abs(lg1)<1e-12) = 0;

% compute means
lxdelta = lx.*ldelta;
lg1mean = lxdelta'*lg1;
lf1mean = lxdelta'*lf1;

% scale to make mean correct
lg = ladj*lf1mean.*lg1./lg1mean;

% adjust mass at zero so its a valid dist
lg(1,:) = (1-ldelta(2:n)'*lg(2:n,:))/ldelta(1);

% adjust it back to have the same sum as on input
lg = lg.*(ldelta'*lf);
lg = reshape(lg,n,m,k);