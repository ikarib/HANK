function [lc,lh,ls,lHc] = OptimalConsumption(lVb,gbdrift,gnetwage,gidioprod)
% set lVb = -999.9 to use budget constraint rather than FOC for consumption
global LaborSupply chi frisch ScaleDisutilityIdio ImposeMaxHours ngpb ngpa utilfn utilfn1 utilfn1inv prefshock

if LaborSupply==0; lh = 1; else; lh = (gnetwage/chi).^frisch; end
if ScaleDisutilityIdio; lh = lh./(gidioprod.^frisch); end
if LaborSupply==1; lh = lh.*lVb.^frisch; end
if ImposeMaxHours; lh = min(lh, 1); end
lc = utilfn1inv(lVb);
ls = gbdrift + lh.*gnetwage - lc;
i = any(lVb <= -999); % stationary point or limit
if any(i(:))
    gbdrift = gbdrift(i);
    gnetwage = repmat(gnetwage,1,ngpb,1);
    gnetwage = gnetwage(i);
    if ScaleDisutilityIdio; gidioprod = squeeze(gidioprod(any(i))); end
    switch LaborSupply
        case 0
            lH = 1;
        case 1
            if ImposeMaxHours; lhmax = 1; else; lhmax = 100; end
            if frisch==1 && prefshock==1
                lhmin = gbdrift./gnetwage/2;
                if ScaleDisutilityIdio
                    lH = sqrt(lhmin.^2+1/chi./gidioprod)-lhmin;
                else
                    lH = sqrt(lhmin.^2+1/chi)-lhmin;
                end
            else
                lhmin = max(0,-gbdrift./gnetwage+1e-5);
                if ScaleDisutilityIdio
                    FnHoursBC = @(lhmin,gbdrift,gnetwage,gidioprod) fzero(@(lh) chi*(lh.^(1/frisch)) - utilfn1(gbdrift + lh.*gnetwage) .* (gnetwage./gidioprod),[lhmin lhmax]);
                    [lH,~,iflag] = arrayfun(FnHoursBC,lhmin,gbdrift,gnetwage,gidioprod);
                else
                    FnHoursBC = @(lhmin,gbdrift,gnetwage) fzero(@(lh) chi*(lh.^(1/frisch)) - utilfn1(gbdrift + lh.*gnetwage) .* gnetwage,[lhmin lhmax]);
                    [lH,~,iflag] = arrayfun(FnHoursBC,lhmin,gbdrift,gnetwage);
                end
                if any(iflag~=1); error('fzero'); end
            end
        case 2        
            lH = (gnetwage/chi).^frisch;
            if ScaleDisutilityIdio; lH = lH./gidioprod.^frisch; end
            if ImposeMaxHours; lH = min(lH, 1); end
    end
    lc(:,i) = repmat(gbdrift + lH.*gnetwage,1,ngpa)';
    lh(:,i) = repmat(lH,1,ngpa)';
    ls(:,i) = 0;
end
llabdisutil = chi*lh.^(1+1/frisch)/(1+1/frisch);
if ScaleDisutilityIdio; llabdisutil = gidioprod.*llabdisutil; end
if LaborSupply==2; ls = ls-llabdisutil; end
lHc = utilfn(lc) + lVb.*ls;
if LaborSupply==1; lHc = lHc-llabdisutil; end
lHc(lc<=0) = -1e12;
if LaborSupply==2; lc = lc+llabdisutil; end