function f = adjcostfn(ld,la)
global kappa0_w kappa1_w kappa2_w kappa0_d kappa1_d kappa2_d kappa3
la = max(la,kappa3);
lx = abs(ld)./la;
f = zeros(size(ld));
id = ld>0; iw = ld<0;
f(id) = kappa0_d*lx(id) + ((kappa1_d^-kappa2_d)/(1+kappa2_d))*lx(id).^(1+kappa2_d);
f(iw) = kappa0_w*lx(iw) + ((kappa1_w^-kappa2_w)/(1+kappa2_w))*lx(iw).^(1+kappa2_w);
f = f.*la;