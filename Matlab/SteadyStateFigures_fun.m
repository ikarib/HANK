function out = SteadyStateFigures_fun(options)

out = 1;

%% Directories
BaseDir         = options.BASEDIR;
BaseSaveDir     = [BaseDir,'/','Output'];
BaseOutputDir   = BaseDir;

%% Options
Experiment              = options.Experiment;
Save                    = options.Save;
datagrosslabinc         = options.datagrosslabinc;

%% Grids
agrid       = load([BaseOutputDir '/agrid.txt']);
dagrid      = diff(agrid);
ngpa        = size(agrid,1);
bgrid       = load([BaseOutputDir '/bgrid.txt']);
dbgrid      = diff(bgrid);
ngpb        = size(bgrid,1);
b0point     = find(bgrid==0);
ygrid       = load([BaseOutputDir '/ygrid.txt']);
ngpy        = size(ygrid,1);
adelta      = load([BaseOutputDir '/adelta.txt']);
bdelta      = load([BaseOutputDir '/bdelta.txt']);
abdelta     = adelta*bdelta';
abydelta    = repmat(abdelta,[1,1,ngpy]);


%% Initial steady state
temp = importdata([BaseOutputDir '/InitialSteadyStateParameters.txt']);
for i = 1:size(temp.data,1)
    initss.(temp.textdata{i}) = temp.data(i,1);
end
initss.priceadjust = 0;

annlabinc = initss.Egrosslabinc.*4;
annoutput = initss.output.*4;

V           = zeros(ngpa,ngpb,ngpy);
dep         = zeros(ngpa,ngpb,ngpy);
con         = zeros(ngpa,ngpb,ngpy);
hour        = zeros(ngpa,ngpb,ngpy);
ccum1       = zeros(ngpa,ngpb,ngpy);
ccum2       = zeros(ngpa,ngpb,ngpy);
ccum4       = zeros(ngpa,ngpb,ngpy);
bdot        = zeros(ngpa,ngpb,ngpy);
gjoint      = zeros(ngpa,ngpb,ngpy);
for iy = 1:ngpy
    V(:,:,iy)       = load([BaseOutputDir '/INITSS/V_INITSS_y' int2str(iy) '.txt']);
    dep(:,:,iy)     = load([BaseOutputDir '/INITSS/dep_INITSS_y' int2str(iy) '.txt']);
    con(:,:,iy)     = load([BaseOutputDir '/INITSS/con_INITSS_y' int2str(iy) '.txt']);
    hour(:,:,iy)     = load([BaseOutputDir '/INITSS/hour_INITSS_y' int2str(iy) '.txt']);
    bdot(:,:,iy)    = load([BaseOutputDir '/INITSS/bdot_INITSS_y' int2str(iy) '.txt']);
    ccum1(:,:,iy)   = load([BaseOutputDir '/INITSS/ccum1_INITSS_y' int2str(iy) '.txt']);
    ccum2(:,:,iy)   = load([BaseOutputDir '/INITSS/ccum2_INITSS_y' int2str(iy) '.txt']);
    ccum4(:,:,iy)   = load([BaseOutputDir '/INITSS/ccum4_INITSS_y' int2str(iy) '.txt']);
    dcum1(:,:,iy)   = load([BaseOutputDir '/INITSS/dcum1_INITSS_y' int2str(iy) '.txt']);
    dcum2(:,:,iy)   = load([BaseOutputDir '/INITSS/dcum2_INITSS_y' int2str(iy) '.txt']);
    dcum4(:,:,iy)   = load([BaseOutputDir '/INITSS/dcum4_INITSS_y' int2str(iy) '.txt']);
    gjoint(:,:,iy)  = load([BaseOutputDir '/INITSS/gjoint_INITSS_y' int2str(iy) '.txt']);    
end    
gamarg          = load([BaseOutputDir '/INITSS/gamarg_INITSS.txt']);    
gbmarg          = load([BaseOutputDir '/INITSS/gbmarg_INITSS.txt']);
gamargallinc    = sum(gamarg,2);
gbmargallinc    = sum(gbmarg,2);
gjointallinc    = sum(gjoint,3);

%%
initss.PERCa = load([BaseOutputDir '/INITSS/PERCa.txt']);    
initss.PERCb = load([BaseOutputDir '/INITSS/PERCb.txt']);    
initss.PERCc = load([BaseOutputDir '/INITSS/PERCc.txt']);    
initss.PERCinc = load([BaseOutputDir '/INITSS/PERCinc.txt']);    
initss.PERCnw  = load([BaseOutputDir '/INITSS/PERCnw.txt']);    

initss.Ea_incQ = load([BaseOutputDir '/INITSS/Ea_incQ.txt']);    
initss.Ea_nwQ = load([BaseOutputDir '/INITSS/Ea_nwQ.txt']);    
initss.Eb_incQ = load([BaseOutputDir '/INITSS/Eb_incQ.txt']);    
initss.Eb_nwQ = load([BaseOutputDir '/INITSS/Eb_nwQ.txt']);    
initss.Ec_incQ = load([BaseOutputDir '/INITSS/Ec_incQ.txt']);    
initss.Ec_nwQ = load([BaseOutputDir '/INITSS/Ec_nwQ.txt']);    
initss.Einc_incQ = load([BaseOutputDir '/INITSS/Einc_incQ.txt']);    
initss.Einc_nwQ = load([BaseOutputDir '/INITSS/Einc_nwQ.txt']);    

%%
ydist       = sum(gbmarg.*(bdelta*ones(1,ngpy)))';
Eb          = sum(gbmarg.*(bgrid.*bdelta*ones(1,ngpy)))'./ydist;
Ea          = sum(gamarg.*(agrid.*adelta*ones(1,ngpy)))'./ydist;
Ed0         = sum(sum((dep==0).*gjoint.*repmat(abdelta,[1,1,ngpy])));
Ed0         = squeeze(Ed0)./ydist;

FRACdNEG    = squeeze(sum(sum((dep<0).*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;
FRACdPOS    = squeeze(sum(sum((dep>0).*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;
FRACd0      = squeeze(sum(sum((dep==0).*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

%%
FRACbNEG    = sum(gbmarg(1:b0point-1,:).*(bdelta(1:b0point-1)*ones(1,ngpy)))'./ydist;
FRACb0      = (gbmarg(b0point,:).*(bdelta(b0point)*ones(1,ngpy)))'./ydist;
FRACbPOS    = sum(gbmarg(b0point+1:ngpb,:).*(bdelta(b0point+1:ngpb)*ones(1,ngpy)))'./ydist;

FRACa0      = adelta(1).*gamarg(1,:) ./ ydist';

% use 5% of average quarterly labor income (approx $750)
b0closepoints   = find(and(bgrid>=0,bgrid<=0.05*initss.Egrosslabinc));
b0farpoints     = find(bgrid>0.05*initss.Egrosslabinc);
FRACb0close     = sum(gbmarg(b0closepoints,:).*(bdelta(b0closepoints)*ones(1,ngpy)))'./ydist; 
FRACb0far       = sum(gbmarg(b0farpoints,:).*(bdelta(b0farpoints)*ones(1,ngpy)))'./ydist; 

bNEG0closepoints    = find(and(bgrid<0,bgrid>=-0.05*initss.Egrosslabinc));
FRACbNEG0close      = sum(gbmarg(bNEG0closepoints,:).*(bdelta(bNEG0closepoints)*ones(1,ngpy)))'./ydist; 

% FRACb0a0 = sum(sum((dep==0).*gjoint.*repmat(abdelta,1,1,ngpy)));

%%
dephistpoints   = [min(min(min(dep))):0.05:max(max(max(dep)))]';
dep_cumdist     = zeros(size(dephistpoints));
for i = 1:length(dephistpoints)
    dep_cumdist(i) = sum(sum(sum(gjoint(dep<=dephistpoints(i)).* abydelta(dep<=dephistpoints(i)) )));
end
dephist         = diff([0; dep_cumdist]);



%%
mpc         = (con(:,2:ngpb,:) - con(:,1:ngpb-1,:))./ repmat(dbgrid',[ngpa,1,ngpy]);
mpc         = [mpc mpc(:,ngpb-1,:)];
Empc        = sum(sum(sum(mpc.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empc_by_y   = squeeze(sum(sum(mpc.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;
Empc_b0close    = sum(sum(sum(mpc(:,b0closepoints,:).*gjoint(:,b0closepoints,:).*repmat(abdelta(:,b0closepoints),[1,1,ngpy])))) ...
                    ./ sum(sum(sum(gjoint(:,b0closepoints,:).*repmat(abdelta(:,b0closepoints),[1,1,ngpy]))));
Empc_b0         = sum(sum(sum(mpc(:,b0point,:).*gjoint(:,b0point,:).*repmat(abdelta(:,b0point),[1,1,ngpy])))) ...
                    ./ sum(sum(sum(gjoint(:,b0point,:).*repmat(abdelta(:,b0point),[1,1,ngpy]))));
Empc_blim       = sum(sum(sum(mpc(:,1,:).*gjoint(:,1,:).*repmat(abdelta(:,1),[1,1,ngpy])))) ...
                    ./ sum(sum(sum(gjoint(:,1,:).*repmat(abdelta(:,1),[1,1,ngpy]))));
Empc_b0far      = sum(sum(sum(mpc(:,b0farpoints,:).*gjoint(:,b0farpoints,:).*repmat(abdelta(:,b0farpoints),[1,1,ngpy])))) ...
                    ./ sum(sum(sum(gjoint(:,b0farpoints,:).*repmat(abdelta(:,b0farpoints),[1,1,ngpy]))));
Empc_bnegint    = sum(sum(sum(mpc(:,1:b0point-1,:).*gjoint(:,1:b0point-1,:).*repmat(abdelta(:,1:b0point-1),[1,1,ngpy])))) ...
                    ./ sum(sum(sum(gjoint(:,1:b0point-1,:).*repmat(abdelta(:,1:b0point-1),[1,1,ngpy]))));

mpchistpoints   = [[0.05:0.05:1]'; max(max(max(mpc)))];
mpc_cumdist     = zeros(size(mpchistpoints));
for i = 1:length(mpchistpoints)
    mpc_cumdist(i) = sum(sum(sum(gjoint(mpc<=mpchistpoints(i)).* abydelta(mpc<=mpchistpoints(i)) )));
end
mpchist         = diff([0; mpc_cumdist]);
mpchistpoints(length(mpchistpoints)) = 1.05;


%%
mpcum1          = (ccum1(:,2:ngpb,:) - ccum1(:,1:ngpb-1,:))./ repmat(dbgrid',[ngpa,1,ngpy]);
mpcum1          = [mpcum1 mpcum1(:,ngpb-1,:)];
Empcum1         = sum(sum(sum(mpcum1.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empcum1_by_y    = squeeze(sum(sum(mpcum1.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

mpcum1histpoints    = [0.05:0.05:1]';
mpcum1_cumdist      = zeros(size(mpcum1histpoints));
for i = 1:length(mpcum1histpoints)
    mpcum1_cumdist(i) = sum(sum(sum(gjoint(mpcum1<=mpcum1histpoints(i)).* abydelta(mpcum1<=mpcum1histpoints(i)) )));
end
mpcum1hist          = diff([0; mpcum1_cumdist]);


%%
mpcum4          = (ccum4(:,2:ngpb,:) - ccum4(:,1:ngpb-1,:))./ repmat(dbgrid',[ngpa,1,ngpy]);
mpcum4          = [mpcum4 mpcum4(:,ngpb-1,:)];
Empcum4         = sum(sum(sum(mpcum4.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empcum4_by_y    = squeeze(sum(sum(mpcum4.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

mpcum4histpoints    = [0.05:0.05:1]';
mpcum4_cumdist      = zeros(size(mpcum4histpoints));
for i = 1:length(mpcum4histpoints)
    mpcum4_cumdist(i) = sum(sum(sum(gjoint(mpcum4<=mpcum4histpoints(i)).* abydelta(mpcum4<=mpcum4histpoints(i)) )));
end
mpcum4hist          = diff([0; mpcum4_cumdist]);

%%
mpd             = (dep(:,2:ngpb,:) - dep(:,1:ngpb-1,:))./ repmat(dbgrid',[ngpa,1,ngpy]);
mpd             = [mpd mpd(:,ngpb-1,:)];
Empd            = sum(sum(sum(mpd.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empd_by_y       = squeeze(sum(sum(mpd.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

%%
rebamount       = (500/datagrosslabinc).* (initss.Egrosslabinc*4);

mpreb1          = zeros(ngpa,ngpb,ngpy);
mpreb4          = zeros(ngpa,ngpb,ngpy);
for ia = 1:ngpa
    for iy = 1:ngpy
        mpreb1(ia,:,iy) = interp1(bgrid',ccum1(ia,:,iy),bgrid'+rebamount,'linear','extrap');
        mpreb1(ia,:,iy) = (mpreb1(ia,:,iy) - ccum1(ia,:,iy))./rebamount;
        mpreb4(ia,:,iy) = interp1(bgrid',ccum4(ia,:,iy),bgrid'+rebamount,'linear','extrap');
        mpreb4(ia,:,iy) = (mpreb4(ia,:,iy) - ccum4(ia,:,iy))./rebamount;
    end
end   

Empreb1         = sum(sum(sum(mpreb1.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empreb1_by_y    = squeeze(sum(sum(mpreb1.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;
Empreb4         = sum(sum(sum(mpreb4.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empreb4_by_y    = squeeze(sum(sum(mpreb4.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

mpreb1histpoints = [0.05:0.05:1]';
mpreb1_cumdist = zeros(size(mpreb1histpoints));
for i = 1:length(mpreb1histpoints)
    mpreb1_cumdist(i) = sum(sum(sum(gjoint(mpreb1<=mpreb1histpoints(i)).* abydelta(mpreb1<=mpreb1histpoints(i)) )));
end
mpreb1hist      = diff([0; mpreb1_cumdist]);

%%
reblargeamount  = (2500/datagrosslabinc).* (initss.Egrosslabinc*4);

mpreblarge1     = zeros(ngpa,ngpb,ngpy);
mpreblarge4     = zeros(ngpa,ngpb,ngpy);
for ia = 1:ngpa
    for iy = 1:ngpy
        mpreblarge1(ia,:,iy) = interp1(bgrid',ccum1(ia,:,iy),bgrid'+reblargeamount,'linear','extrap');
        mpreblarge1(ia,:,iy) = (mpreblarge1(ia,:,iy) - ccum1(ia,:,iy))./reblargeamount;
        mpreblarge4(ia,:,iy) = interp1(bgrid',ccum4(ia,:,iy),bgrid'+reblargeamount,'linear','extrap');
        mpreblarge4(ia,:,iy) = (mpreblarge4(ia,:,iy) - ccum4(ia,:,iy))./reblargeamount;
    end
end   

Empreblarge1        = sum(sum(sum(mpreblarge1.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empreblarge1_by_y   = squeeze(sum(sum(mpreblarge1.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;
Empreblarge4        = sum(sum(sum(mpreblarge4.*gjoint.*repmat(abdelta,[1,1,ngpy]))));
Empreblarge4_by_y   = squeeze(sum(sum(mpreblarge4.*gjoint.*repmat(abdelta,[1,1,ngpy]))))./ydist;

mpreblarge1histpoints = [0.05:0.05:1]';
mpreblarge1_cumdist = zeros(size(mpreblarge1histpoints));
for i = 1:length(mpreblarge1histpoints)
    mpreblarge1_cumdist(i) = sum(sum(sum(gjoint(mpreblarge1<=mpreblarge1histpoints(i)).* abydelta(mpreblarge1<=mpreblarge1histpoints(i)) )));
end
mpreblarge1hist     = diff([0; mpreblarge1_cumdist]);





%% Save the workspace
save(['Steadystate_workspace.mat']);



