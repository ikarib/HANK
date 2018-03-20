global utilfn utilfninv utilfn1 utilfn1inv
if gam==1
    utilfn = @(lx) prefshock*log(lx);
    utilfninv = @(lu) exp(lu/prefshock);
    utilfn1 = @(lx) prefshock./lx;
    utilfn1inv = @(lu) prefshock./lu;
else
    utilfn = @(lx) prefshock*(lx.^(1-gam)-1)/(1-gam);
    utilfninv = @(lu) ((1-gam)*(lu/prefshock)+1).^(1/(1-gam));
    utilfn1 = @(lx) prefshock*lx.^-gam;
    utilfn1inv = @(lu) (lu/prefshock).^(-1/gam);
end
WorldBondFunction2 = @(r,bss,rss,e) bss + e*(r-rss); % b = f(r)
WorldBondInverse2 = @(b,bss,rss,e) rss + (b-bss)/e; % r = f^-1(b)
PartialUpdate = @(n,lstep,xn,fn) lstep*fn + (1-lstep)*xn;