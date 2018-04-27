function INITSS=DiscountedMPC(c,d,AU)
global agrid bgrid dagrid dbgrid utilfn1 Display

if Display; fprintf('Computing discounted MPC and wealth in steady state\n'); end

% compute MPCs and MPDs
mpc = diff(c,1,2)./dbgrid; mpc = [mpc mpc(:,end,:)]; INITSS.mpc = mpc;
mpd_a = diff(d,1,1)./dagrid; mpd_a = [mpd_a; mpd_a(end,:,:)];
mpd_b = diff(d,1,2)./dbgrid; mpd_b = [mpd_b  mpd_b(:,end,:)];
effdisc = mpc+(1+adjcostfn1(d,agrid)).*mpd_b - mpd_a;
muc = utilfn1(c);
badj = bgrid.*(mpc./c).*muc;

% SUBSTITUTION EFFECT: ONE ASSET: effective discounting of xi (see notes) and decay eta
if Display>1; fprintf('Sub effect 1 asset:\n'); end
INITSS.subeff1ass = effect(mpc,muc);

% SUBSTITUTION EFFECT: TWO ASSET: effective discounting of (MPC + (1+ chi'(d)*MPD_B-MPD_A) and decay eta
if Display>1; fprintf('Sub effect 2 asset:\n'); end
INITSS.subeff2ass = effect(effdisc,muc);

% WEALTH EFFECT: ONE ASSET: effective discounting of xi and decay eta
if Display>1; fprintf('Wealth effect 1 asset:\n'); end
INITSS.wealtheff1ass = effect(mpc,badj);

% WEALTH EFFECT: TWO ASSET: effective discounting of (MPC + (1+ chi'(d)*MPD_B-MPD_A) and decay eta
if Display>1; fprintf('Wealth effect 2 asset:\n'); end
INITSS.wealtheff2ass = effect(effdisc,badj);

function lvec = effect(disc,rhs)
    global nab ngpa ngpb ngpy MonetaryShockPers rho rb ymarkovoff ymarkovdiag
    maxiter_discMPC = 50;
    discMPCtol = 1e-10;
    eta = -log(MonetaryShockPers);

    % have to iterate because need to solve for all y points simultaneously
    disc = reshape(disc,nab,1,ngpy);
    rhs = reshape(rhs,nab,ngpy);
    % initialize iterations without y transitions
    B = [sum(AU,2)+eta+rho-rb+disc -AU];
    lvec = rhs;
    for iy=1:ngpy %par
        lvec(:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\lvec(:,iy);
    end
    B = [sum(AU,2)+eta+rho-rb+disc-ymarkovdiag -AU];
    % iterate to convergence
    for it = 1:maxiter_discMPC
        f1 = lvec;
        lvec = rhs + lvec*ymarkovoff';
        for iy = 1:ngpy %par
            lvec(:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\lvec(:,iy);
        end
        ldiff = max(max(abs(lvec-f1)));
        if Display>1; fprintf('iteration %d max change: %g\n',it,ldiff); end
        if ldiff<discMPCtol; break; end
    end
    lvec = reshape(lvec,ngpa,ngpb,ngpy);
end

end