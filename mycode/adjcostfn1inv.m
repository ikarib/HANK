function f = adjcostfn1inv(lchi,la)
global kappa0_w kappa1_w kappa2_w kappa0_d kappa1_d kappa2_d kappa3 dmax Display
la = max(la,kappa3);
id = lchi>kappa0_d; iw = lchi<-kappa0_w;
f = zeros(size(lchi));
f(id) = kappa1_d * (lchi(id)-kappa0_d).^(1/kappa2_d);
f(iw) = -kappa1_w * (-lchi(iw)-kappa0_w).^(1/kappa2_w);
f = f.*la;
imax = f>dmax;
if any(any(any(imax(:,1:end-1,:),3)))
    if Display>1
        fprintf('high optimal d: %e\n',f(imax));
    end
    f(imax) = dmax;
end
