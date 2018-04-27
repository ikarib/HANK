function irf = DiscountedMPCTransition(irf,equmINITSS,solnFINALSS)
global agrid bgrid dagrid dbgrid utilfn1 Display Ttransition deltatransvec

if Display; fprintf(' Solving discounted MPCs backward\n'); end

for it = Ttransition:-1:1
    if Display>1; fprintf('   Solving discounted MPCs backward at time: %g\n',it); end
    delta = deltatransvec(it);

%     rho = irf.equmSTICKY.rho(it);
    rb = irf.equmSTICKY.rb(it);
    ldrb = irf.equmSTICKY.rb(it) - equmINITSS.rb;
    c = irf.solnSTICKY.c{it};
    d = irf.solnSTICKY.d{it};
    AU = irf.solnSTICKY.AU{it};
    
    if it==Ttransition
        subeff1ass = solnFINALSS.subeff1ass;
        subeff2ass = solnFINALSS.subeff2ass;
        wealtheff1ass = solnFINALSS.wealtheff1ass;
        wealtheff2ass = solnFINALSS.wealtheff2ass;
    else
        subeff1ass = irf.solnSTICKY.subeff1ass{it+1};
        subeff2ass = irf.solnSTICKY.subeff2ass{it+1};
        wealtheff1ass = irf.solnSTICKY.wealtheff1ass{it+1};
        wealtheff2ass = irf.solnSTICKY.wealtheff2ass{it+1};
    end

    % compute MPCs and MPDs
    mpc = diff(c,1,2)./dbgrid; mpc = [mpc mpc(:,end,:)]; irf.solnSTICKY.mpc{it} = mpc;
    mpd_a = diff(d,1,1)./dagrid; mpd_a = [mpd_a; mpd_a(end,:,:)];
    mpd_b = diff(d,1,2)./dbgrid; mpd_b = [mpd_b  mpd_b(:,end,:)];
    effdisc = mpc+(1+adjcostfn1(d,agrid)).*mpd_b - mpd_a;
    muc = utilfn1(c);
    badj = bgrid.*(mpc./c).*muc;

    % SUBSTITUTION EFFECT: ONE ASSET
    irf.solnSTICKY.subeff1ass{it} = effect(mpc,muc,subeff1ass);

    % SUBSTITUTION EFFECT: TWO ASSET
    irf.solnSTICKY.subeff2ass{it} = effect(effdisc,muc,subeff2ass);

    % WEALTH EFFECT: ONE ASSET
    irf.solnSTICKY.wealtheff1ass{it} = effect(mpc,badj,wealtheff1ass);
    
    % WEALTH EFFECT: TWO ASSET: effective discounting of MPC and decay leta
    irf.solnSTICKY.wealtheff2ass{it} = effect(effdisc,badj,wealtheff2ass);
end

function lvec = effect(disc,rhs1,lvec)
    global nab ngpa ngpb ngpy rho ymarkovoff ymarkovdiag

    % have to iterate because need to solve for all y points simultaneously
    disc = reshape(disc,nab,1,ngpy);
    rhs1 = reshape(rhs1,nab,ngpy);
    lvec = reshape(lvec,nab,ngpy);
    % initialize iterations without y transitions
    B = [1+delta*(sum(AU,2)+rho-rb+disc-ymarkovdiag) -delta*AU];
    lvec = delta*rhs1*ldrb + lvec*(eye(ngpy) + delta*ymarkovoff');
    for iy = 1:ngpy %par
        lvec(:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\lvec(:,iy);
    end
    lvec = reshape(lvec,ngpa,ngpb,ngpy);
end

end