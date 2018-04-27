function [Vnew,c,h,s,d,dadj] = HJBUpdate(V,adrift,bdrift,netwage,lumptransfer,labtax,profit,meanlabeff,delta)
global ngpa ngpb ngpy nab dagrid dbgrid dVamin dVbmin agrid ygrid 
global DistributeProfitsInProportion TaxHHProfitIncome profdistfrac corptax
global LaborSupply chi frisch ScaleDisutilityIdio utilfn utilfn1 ymarkovdiag ymarkovoff rho deathrate ReportNonMonotonicity
% Hamilton-Jacobi-Bellman equation

c = zeros(ngpa,ngpb,ngpy); h = c; s = c; d = c;

% derivatives wrt a: forward and backward
lVaF = max(diff(V,1,1)./dagrid,dVamin);
lVaB(2:ngpa,:,:) = lVaF;
lVaF(ngpa,:,:) = dVamin;

% derivatives wrt b: forward and backward
lVbF = max(diff(V,1,2)./dbgrid,dVbmin);
lVbB(:,2:ngpb,:) = lVbF;
lVbF(:,ngpb,:) = dVbmin;

gidioprod = ygrid;
gnetwage = netwage*gidioprod;
gbdrift = bdrift + lumptransfer;
if DistributeProfitsInProportion
    if TaxHHProfitIncome
        gbdrift = gbdrift + (1-labtax)*(1-profdistfrac)*(1-corptax)*profit/meanlabeff*ygrid;
    else
        gbdrift = gbdrift + (1-profdistfrac)*(1-corptax)*profit/meanlabeff*ygrid;
    end
end

% consumption decision
[lcF,lhF,lsF,lHcF] = OptimalConsumption(lVbF,gbdrift,gnetwage,gidioprod);
lsF(:,ngpb,:) = 0;
lHcF(:,ngpb,:) = -1e12;
validF = lsF>0;

lVbB(:,1,:) = -999.9;
[lcB,lhB,lsB,lHcB] = OptimalConsumption(lVbB,gbdrift,gnetwage,gidioprod);
validB = lsB<0;

[lc0,lh0,ls0,lHc0] = OptimalConsumption(repmat(-999.9,ngpa,ngpb,ngpy),gbdrift,gnetwage,gidioprod);

iF = validF & (~validB | lHcF>=lHcB) & lHcF>=lHc0; % forward
c(iF) = lcF(iF);
h(iF) = lhF(iF);
s(iF) = lsF(iF);
iB = validB & (~validF | lHcB>=lHcF) & lHcB>=lHc0; % backward
c(iB) = lcB(iB);
h(iB) = lhB(iB);
s(iB) = lsB(iB);
i0 = ~iF & ~iB;
c(i0) = lc0(i0);
h(i0) = lh0(i0);
s(i0) = ls0(i0);

% deposit decision
ldFB = adjcostfn1inv(lVaF./lVbB-1,agrid); % a forward, b backward
lHdFB = lVaF.*ldFB - lVbB.*(ldFB + adjcostfn(ldFB,agrid));
lHdFB(ngpa,:,:) = -1e12;
lHdFB(:,1,:) = -1e12;
validFB = ldFB>0 & lHdFB>0;

ldBF = adjcostfn1inv(lVaB./lVbF-1,agrid); % a backward, b forward
ldBFadj = ldBF + adjcostfn(ldBF,agrid);
lHdBF = lVaB.*ldBF - lVbF.*ldBFadj;
lHdBF(1,:,:) = -1e12;
lHdBF(:,ngpb,:) = -1e12;
validBF = ldBFadj<=0 & lHdBF>0;

ldBB = adjcostfn1inv(lVaB./lVbB-1,agrid); % a backward, b backward
ldBBadj = ldBB + adjcostfn(ldBB,agrid);
lHdBB = lVaB.*ldBB - lVbB.*ldBBadj;
lHdBB(1,:,:) = -1e12;
validBB = ldBBadj>0 & ldBB<=0 & lHdBB>0;

llabdisutil = lhB(:,1,:).^(1+1/frisch)/(1+1/frisch);
if ScaleDisutilityIdio; llabdisutil = ygrid.*llabdisutil; end
if LaborSupply==0 || LaborSupply==1; lVbB(:,1,:) = utilfn1(lcB(:,1,:)); end
if LaborSupply==2; lVbB(:,1,:) = utilfn1(lcB(:,1,:)-chi*llabdisutil); end

iFB = validFB & (~validBF | lHdFB>=lHdBF) & (~validBB | lHdFB>=lHdBB); % forward, backward
iBF = validBF & (~validFB | lHdBF>=lHdFB) & (~validBB | lHdBF>=lHdBB); % backward,forward
iBB = validBB & (~validFB | lHdBB>=lHdFB) & (~validBF | lHdBB>=lHdBF); % backward,backward
d(iFB) = ldFB(iFB);
d(iBF) = ldBF(iBF);
d(iBB) = ldBB(iBB);
d(~validFB & ~validBF & ~validBB) = 0; % none
if any(any(any(~iFB & ~iBF & ~iBB & (validFB | validBF | validBB)))) % more than 1
    error('should never be here!')
end

llabdisutil = h.^(1+1/frisch)/(1+1/frisch);
if ScaleDisutilityIdio; llabdisutil = ygrid.*llabdisutil; end
switch LaborSupply
    case 0; u = utilfn(c);
    case 1; u = utilfn(c) - chi*llabdisutil;
    case 2; u = utilfn(c - chi*llabdisutil);
end
dadj = d + adjcostfn(d,agrid);

% vector of constants
lbvec = reshape(delta*u + V + delta*shiftdim(reshape(ymarkovoff*reshape(shiftdim(V,2),ngpy,nab),ngpy,ngpa,ngpb),1),nab,ngpy);

% a drifts
ladriftB = min(d,0) + min(adrift,0);
ladriftF = max(d,0) + max(adrift,0);

% b drift
lbdriftB = min(-dadj,0) + min(s,0);
lbdriftF = max(-dadj,0) + max(s,0);

% construct A matrix in sparse form and vectors by filling in coefficients on V
% construct B matrix = I +delta*(rho*I -A): assumes that all diagonal terms in A matrix are non-zero
% do not impose that diagonal elements are non-zero, keep them even if zero since B matrix below will modify them
A = reshape([  ladriftF./[dagrid;1] ... % a+1
              -ladriftB./[1;dagrid] ... % a-1
               lbdriftF./[dbgrid 1] ... % b+1
              -lbdriftB./[1 dbgrid] ... % b-1
            ], nab,4,ngpy);
B = [1+delta*(sum(A,2)+rho+deathrate-ymarkovdiag) -delta*A]; % diagonals
for iy = 1:ngpy %par
    lbvec(:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\lbvec(:,iy);
end
Vnew = reshape(lbvec,ngpa,ngpb,ngpy);

% check for non-monotonicity
if ReportNonMonotonicity
    if any(any(any(diff(Vnew,1,1)<-1e-8)))
        error('non-monotonicity in a dimension')
    end
    if any(any(any(diff(Vnew,1,2)<-1e-8)))
        error('non-monotonicity in b dimension')
    end
end
