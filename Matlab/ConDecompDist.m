clc;
clear;
close all;

%% PARAMETERS

BaseInputDir = '~/FortranOutputDir/BaselineOutputSubdir/'; 
SaveFigsDir = '~/FiguresDir'; 


%% OPTIONS
STIM                = {'NOFS','PE1','PE2','PE3','PE4','PE5','PE6','PE7','PE8','PE9','PE10','PE11','PE12','PE13','PE14','PE15'};
explist = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
InputDirAlt = BaseInputDir;
lSave = 1;
UseSteadyStateDecomp = 0;
Smoothing = 1;
tstep       = load([BaseInputDir '/deltatransvec.txt']);

%% Grids
agrid       = load([BaseInputDir '/agrid.txt']);
dagrid      = diff(agrid);
ngpa        = size(agrid,1);
bgrid       = load([BaseInputDir '/bgrid.txt']);
dbgrid      = diff(bgrid);
ngpb        = size(bgrid,1);
b0point     = find(bgrid==0);
ygrid       = load([BaseInputDir '/ygrid.txt']);
ngpy        = size(ygrid,1);
adelta      = load([BaseInputDir '/adelta.txt']);
bdelta      = load([BaseInputDir '/bdelta.txt']);
abdelta     = adelta*bdelta';
abydelta    = repmat(abdelta,[1,1,ngpy]);

% datagrosslabinc = 69100;
dataannoutput = 115000;

% Load initial steady state policy functions
INITSS.V           = zeros(ngpa,ngpb,ngpy);
INITSS.dep         = zeros(ngpa,ngpb,ngpy);
INITSS.con         = zeros(ngpa,ngpb,ngpy);
INITSS.bdot        = zeros(ngpa,ngpb,ngpy);
INITSS.sav         = zeros(ngpa,ngpb,ngpy);
INITSS.gjoint      = zeros(ngpa,ngpb,ngpy);
INITSS.mpc         = zeros(ngpa,ngpb,ngpy);
INITSS.subeff1ass     = zeros(ngpa,ngpb,ngpy);
INITSS.subeff2ass     = zeros(ngpa,ngpb,ngpy);
INITSS.wealtheff1ass     = zeros(ngpa,ngpb,ngpy);
INITSS.wealtheff2ass     = zeros(ngpa,ngpb,ngpy);

for iy = 1:ngpy
    INITSS.V(:,:,iy)       = load([BaseInputDir '/INITSS/V_INITSS_y' int2str(iy) '.txt']);
    INITSS.dep(:,:,iy)     = load([BaseInputDir '/INITSS/dep_INITSS_y' int2str(iy) '.txt']);
    INITSS.con(:,:,iy)     = load([BaseInputDir '/INITSS/con_INITSS_y' int2str(iy) '.txt']);
    INITSS.bdot(:,:,iy)    = load([BaseInputDir '/INITSS/bdot_INITSS_y' int2str(iy) '.txt']);
    INITSS.sav(:,:,iy)    = load([BaseInputDir '/INITSS/sav_INITSS_y' int2str(iy) '.txt']);
    INITSS.gjoint(:,:,iy)  = load([BaseInputDir '/INITSS/gjoint_INITSS_y' int2str(iy) '.txt']);    
    INITSS.mpc(:,:,iy)    = load([BaseInputDir '/INITSS/mpc_INITSS_y' int2str(iy) '.txt']);
    if UseSteadyStateDecomp==0
        INITSS.subeff1ass(:,:,iy)    = load([BaseInputDir '/INITSS/subeff1ass_INITSS_y' int2str(iy) '.txt']);
        INITSS.subeff2ass(:,:,iy)    = load([BaseInputDir '/INITSS/subeff2ass_INITSS_y' int2str(iy) '.txt']);
        INITSS.wealtheff1ass(:,:,iy)    = load([BaseInputDir '/INITSS/wealtheff1ass_INITSS_y' int2str(iy) '.txt']);
        INITSS.wealtheff2ass(:,:,iy)    = load([BaseInputDir '/INITSS/wealtheff2ass_INITSS_y' int2str(iy) '.txt']);
    else    
        INITSS.subeff1ass(:,:,iy)    = load([InputDirAlt '/INITSS/subeff1ass_INITSS_y' int2str(iy) '.txt']);
        INITSS.subeff2ass(:,:,iy)    = load([InputDirAlt '/INITSS/subeff2ass_INITSS_y' int2str(iy) '.txt']);
        INITSS.wealtheff1ass(:,:,iy)    = load([InputDirAlt '/INITSS/wealtheff1ass_INITSS_y' int2str(iy) '.txt']);
        INITSS.wealtheff2ass(:,:,iy)    = load([InputDirAlt '/INITSS/wealtheff2ass_INITSS_y' int2str(iy) '.txt']);
    end

    
end    


%%
load([BaseInputDir '/Steadystate_workspace.mat']);

NOFS = load([BaseInputDir '/IRF_Monetary_NOFS_workspace.mat']);
for e =explist
    PE{e} = load([BaseInputDir '/IRF_Monetary_PE' int2str(e) '_workspace.mat']);
end
%% Load time 1 policy functions
for i = 1:numel(STIM)
    InputDir = [BaseInputDir '/IRF_Monetary/' STIM{i} '/STICKY'];

    eval([STIM{i} '.V           = zeros(ngpa,ngpb,ngpy);']);
    eval([STIM{i} '.dep           = zeros(ngpa,ngpb,ngpy);']);
    eval([STIM{i} '.con           = zeros(ngpa,ngpb,ngpy);']);
    eval([STIM{i} '.bdot           = zeros(ngpa,ngpb,ngpy);']);
    eval([STIM{i} '.sav           = zeros(ngpa,ngpb,ngpy);']);
    eval([STIM{i} '.gjoint         = zeros(ngpa,ngpb,ngpy);']);
    eval([STIM{i} '.mpc          = zeros(ngpa,ngpb,ngpy);']);
    eval([STIM{i} '.subeff1ass          = zeros(ngpa,ngpb,ngpy);']);
    eval([STIM{i} '.subeff2ass          = zeros(ngpa,ngpb,ngpy);']);
    eval([STIM{i} '.wealtheff1ass          = zeros(ngpa,ngpb,ngpy);']);
    eval([STIM{i} '.wealtheff2ass          = zeros(ngpa,ngpb,ngpy);']);
    
end    
for iy = 1:ngpy
    InputDir = [BaseInputDir '/IRF_Monetary/NOFS/STICKY'];
    NOFS.V(:,:,iy)       = load([InputDir '/V_T1_y' int2str(iy) '.txt']);
    NOFS.dep(:,:,iy)     = load([InputDir '/dep_T1_y' int2str(iy) '.txt']);
    NOFS.con(:,:,iy)     = load([InputDir '/con_T1_y' int2str(iy) '.txt']);
    NOFS.bdot(:,:,iy)    = load([InputDir '/bdot_T1_y' int2str(iy) '.txt']);
    NOFS.sav(:,:,iy)    = load([InputDir '/sav_T1_y' int2str(iy) '.txt']);
    NOFS.gjoint(:,:,iy)  = load([InputDir '/gjoint_T1_y' int2str(iy) '.txt']);    
    NOFS.mpc(:,:,iy)    = load([InputDir '/mpc_T1_y' int2str(iy) '.txt']);
    if UseSteadyStateDecomp==0    
        NOFS.subeff1ass(:,:,iy)    = load([InputDir '/subeff1ass_T1_y' int2str(iy) '.txt']);
        NOFS.subeff2ass(:,:,iy)    = load([InputDir '/subeff2ass_T1_y' int2str(iy) '.txt']);
        NOFS.wealtheff1ass(:,:,iy)    = load([InputDir '/wealtheff1ass_T1_y' int2str(iy) '.txt']);
        NOFS.wealtheff2ass(:,:,iy)    = load([InputDir '/wealtheff2ass_T1_y' int2str(iy) '.txt']);
    else 
        NOFS.subeff1ass(:,:,iy)    = load([InputDirAlt '/INITSS/subeff1ass_INITSS_y' int2str(iy) '.txt']);
        NOFS.subeff2ass(:,:,iy)    = load([InputDirAlt '/INITSS/subeff2ass_INITSS_y' int2str(iy) '.txt']);
        NOFS.wealtheff1ass(:,:,iy)    = load([InputDirAlt '/INITSS/wealtheff1ass_INITSS_y' int2str(iy) '.txt']);
        NOFS.wealtheff2ass(:,:,iy)    = load([InputDirAlt '/INITSS/wealtheff2ass_INITSS_y' int2str(iy) '.txt']);
    end 
    
    for e =explist

        InputDir = [BaseInputDir '/IRF_Monetary/PE' int2str(e) '/STICKY'];
        PE{e}.V(:,:,iy)       = load([InputDir '/V_T1_y' int2str(iy) '.txt']);
        PE{e}.dep(:,:,iy)     = load([InputDir '/dep_T1_y' int2str(iy) '.txt']);
        PE{e}.con(:,:,iy)     = load([InputDir '/con_T1_y' int2str(iy) '.txt']);
        PE{e}.bdot(:,:,iy)    = load([InputDir '/bdot_T1_y' int2str(iy) '.txt']);
        PE{e}.sav(:,:,iy)    = load([InputDir '/sav_T1_y' int2str(iy) '.txt']);
        PE{e}.gjoint(:,:,iy)  = load([InputDir '/gjoint_T1_y' int2str(iy) '.txt']);    
        PE{e}.mpc(:,:,iy)     = load([InputDir '/mpc_T1_y' int2str(iy) '.txt']);
       if UseSteadyStateDecomp==0    
            PE{e}.subeff1ass(:,:,iy)    = load([InputDir '/subeff1ass_T1_y' int2str(iy) '.txt']);
            PE{e}.subeff2ass(:,:,iy)    = load([InputDir '/subeff2ass_T1_y' int2str(iy) '.txt']);
            PE{e}.wealtheff1ass(:,:,iy)    = load([InputDir '/wealtheff1ass_T1_y' int2str(iy) '.txt']);
            PE{e}.wealtheff2ass(:,:,iy)    = load([InputDir '/wealtheff2ass_T1_y' int2str(iy) '.txt']);
        else
            PE{e}.subeff1ass(:,:,iy)    = load([InputDirAlt '/INITSS/subeff1ass_INITSS_y' int2str(iy) '.txt']);
            PE{e}.subeff2ass(:,:,iy)    = load([InputDirAlt '/INITSS/subeff2ass_INITSS_y' int2str(iy) '.txt']);
            PE{e}.wealtheff1ass(:,:,iy)    = load([InputDirAlt '/INITSS/wealtheff1ass_INITSS_y' int2str(iy) '.txt']);
            PE{e}.wealtheff2ass(:,:,iy)    = load([InputDirAlt '/INITSS/wealtheff2ass_INITSS_y' int2str(iy) '.txt']);
        end
    end

end    



%%
Save = lSave;

%% change in prices
drb = (NOFS.sticky.rb(2)-initss.rb);
dra = (NOFS.sticky.ra(2)-initss.ra);

%% smooth effects of upwinding direction change on consumption and deposits
for e = explist
    PE{e}.consmooth = PE{e}.con;
    PE{e}.initconsmooth = INITSS.con;
    
for ia = 1:ngpa     %17: 17% 1:ngpa    
for iy = 1:ngpy %21:21 %1:ngpy
    
    %find points where bdot=0 at adjacent points 
    clear ib1 ib2 ib3;
    ib1 = find((INITSS.sav(ia,:,iy)==0).*(PE{e}.sav(ia,:,iy)~=0)==1);
    ib2 = find((INITSS.sav(ia,:,iy)~=0).*(PE{e}.sav(ia,:,iy)==0)==1);
    ib3 = find((INITSS.sav(ia,:,iy)==0).*(PE{e}.sav(ia,:,iy)==0)==1);
    
    ibvec = [ib1 ib2 ib3];
    ibL = min(ibvec);
    ibH = max(ibvec);
    
    %positive points
    if isempty(ibvec)==0 && ibL>11 && ibH<50
     
        if ibH - ibL <=2
            for ib = ibL: ibH 
                PE{e}.consmooth(ia,ib,iy) = interp1(bgrid([ibL-1 ibH+1]),PE{e}.con(ia,[ibL-1 ibH+1],iy),bgrid(ib));
                PE{e}.initconsmooth(ia,ib,iy) = interp1(bgrid([ibL-1 ibH+1]),INITSS.con(ia,[ibL-1 ibH+1],iy),bgrid(ib));
            end    
        end
    end
    
    %negative points
    if isempty(ibvec)==0 && ibL>1 && ibH<9
        if ibH - ibL <=2
            for ib = ibL: ibH 
                PE{e}.consmooth(ia,ib,iy) = interp1(bgrid([ibL-1 ibH+1]),PE{e}.con(ia,[ibL-1 ibH+1],iy),bgrid(ib));
                PE{e}.initconsmooth(ia,ib,iy) = interp1(bgrid([ibL-1 ibH+1]),INITSS.con(ia,[ibL-1 ibH+1],iy),bgrid(ib));
            end    
        end
    end
    
end
end
end


%% consumption elasticities

% consumption elasticities and components at every point in statespace

for e =explist
    elastc{e} = (log(PE{e}.consmooth) - log(PE{e}.initconsmooth) ) ./(-drb);   
    direct_sub{e} = (1./initss.gam) .* (PE{e}.initconsmooth.^initss.gam) .* PE{e}.subeff2ass;
    direct_wealth{e} = - PE{e}.wealtheff2ass .* (PE{e}.initconsmooth.^initss.gam);
    direct_tot{e} = direct_sub{e} + direct_wealth{e};
    direct_port{e} = direct_tot{e} ...
                    -(1./initss.gam) .* (PE{e}.initconsmooth.^initss.gam) .* PE{e}.subeff1ass...
                    + PE{e}.wealtheff1ass .* (PE{e}.initconsmooth.^initss.gam);
end


% consumption weights for aggregation
ctot_init = sum(sum(sum(INITSS.con.*abydelta.*gjoint)));
for e =explist
%     ctot_pe{e} = sum(sum(sum(PE{e}.con.*abydelta.*gjoint)));
    ctot_pe{e} = sum(sum(sum(PE{e}.consmooth.*abydelta.*gjoint)));
end
cshare_aby = (INITSS.con.*abydelta.* gjoint)./ctot_init;
cshare_b = sum(sum(cshare_aby,1),3)';
cshare_a = sum(sum(cshare_aby,2),3);
cshare_ab = sum(cshare_aby,3);
cshare_a_in_b = cshare_ab./repmat(cshare_b',[ngpa,1]);



cshare_aby_in_b = cshare_aby./repmat(cshare_b',[ngpa,1,ngpy]);
cshare_aby_in_ab = cshare_aby./repmat(cshare_ab,[1,1,ngpy]);

% consumption elasticities and components by liquid wealth

for e =explist
    elastc_b{e} = sum(sum(elastc{e}.*cshare_aby_in_b,1),3)';
    direct_sub_b{e} = sum(sum(direct_sub{e}.*cshare_aby_in_b,1),3)';
    direct_port_b{e} = sum(sum(direct_port{e}.*cshare_aby_in_b,1),3)';
    direct_wealth_b{e} = sum(sum(direct_wealth{e}.*cshare_aby_in_b,1),3)';
    direct_tot_b{e} = sum(sum(direct_tot{e}.*cshare_aby_in_b,1),3)';

end

% consumption elasticities and components by liquid wealth, conditonal on illiquid wealth
for e =explist
    elastc_ab{e} = sum(elastc{e}.*cshare_aby_in_ab,3);
    direct_sub_ab{e} = sum(direct_sub{e}.*cshare_aby_in_ab,3);
    direct_port_ab{e} = sum(direct_port{e}.*cshare_aby_in_ab,3);
    direct_wealth_ab{e} = sum(direct_wealth{e}.*cshare_aby_in_ab,3);
    direct_tot_ab{e} = sum(direct_tot{e}.*cshare_aby_in_ab,3);
end

% consumption elasticities and components overall
for e =explist
    elastc_all{e} = sum(elastc{e}(:).*cshare_aby(:));
    direct_sub_all{e} = sum(direct_sub{e}(:).*cshare_aby(:));
    direct_port_all{e} = sum(direct_port{e}(:).*cshare_aby(:));
    direct_wealth_all{e} = sum(direct_wealth{e}(:).*cshare_aby(:));
    direct_tot_all{e} = sum(direct_tot{e}(:).*cshare_aby(:));
end
%% liquid wealth consumption share histogratm

liqspace = 1; %thousands;
liqcdf = cumsum(cshare_b);
liqbinmin = floor(bgrid(1).*dataannoutput/(4*initss.output*1000));
liqbinmax = ceil(bgrid(ngpb).*dataannoutput/(4*initss.output*1000));
liqbin = [liqbinmin:liqspace:liqbinmax]';
nliqbin = size(liqbin,1);
liqbincdf = interp1(bgrid.*dataannoutput/(4*initss.output*1000),liqcdf,liqbin);
liqbincdf(1) = 0;
liqbincdf(nliqbin) = 1;
liqbinpdf = [liqbincdf(1:nliqbin)] - [0;liqbincdf(1:nliqbin-1)];

maxliqhist = 150;
% maxliqhist = 750;
liqbintrunc = liqbin(liqbin<=maxliqhist);
liqbinpdftrunc = liqbinpdf(liqbin<=maxliqhist);
liqbinpdftrunc(end) = liqbinpdftrunc(end) + sum(liqbinpdf(liqbin>maxliqhist));


%% moving average smoothing because aggregation weights
for e =explist
    
    if Smoothing ==0
        sm_elastc_b{e} = elastc_b{e};
        sm_direct_sub_b{e} = direct_sub_b{e};
        sm_direct_wealth_b{e} = direct_wealth_b{e};
        sm_direct_tot_b{e} = direct_tot_b{e};
        sm_direct_port_b{e} = direct_port_b{e};
    end

    if Smoothing ==1
        sm_elastc_b{e}(1:10) = smooth(bgrid(1:10),elastc_b{e}(1:10),0.3,'loess');
        sm_elastc_b{e}(11) = elastc_b{e}(11);
        sm_elastc_b{e}(12:ngpb) = smooth(bgrid(12:ngpb),elastc_b{e}(12:ngpb),0.4,'loess');

        sm_direct_sub_b{e}(1:10) = smooth(bgrid(1:10),direct_sub_b{e}(1:10),0.3,'loess');
        sm_direct_sub_b{e}(11) = direct_sub_b{e}(11);
        sm_direct_sub_b{e}(12:ngpb) = smooth(bgrid(12:ngpb),direct_sub_b{e}(12:ngpb),0.4,'loess');

        sm_direct_wealth_b{e}(1:10) = smooth(bgrid(1:10),direct_wealth_b{e}(1:10),0.3,'loess');
        sm_direct_wealth_b{e}(11) = direct_wealth_b{e}(11);
        sm_direct_wealth_b{e}(12:ngpb) = smooth(bgrid(12:ngpb),direct_wealth_b{e}(12:ngpb),0.4,'loess');

        sm_direct_tot_b{e}(1:10) = smooth(bgrid(1:10),direct_tot_b{e}(1:10),0.3,'loess');
        sm_direct_tot_b{e}(11) = direct_tot_b{e}(11);
        sm_direct_tot_b{e}(12:ngpb) = smooth(bgrid(12:ngpb),direct_tot_b{e}(12:ngpb),0.4,'loess');

        sm_direct_port_b{e}(1:10) = smooth(bgrid(1:10),direct_port_b{e}(1:10),0.3,'loess');
        sm_direct_port_b{e}(11:11) = direct_port_b{e}(11:11);
        sm_direct_port_b{e}(12:ngpb) = smooth(bgrid(12:ngpb),direct_port_b{e}(12:ngpb),0.4,'loess');
    end
end


%% FIGURES FOR PAPER
%% Elasticity C figure: Total elasticity
%8 is all at once
figure;

[fbar]   = bar(liqbintrunc, liqbinpdftrunc,'histc');
axbar = gca;
sh  = findall(gcf,'marker','*'); delete(sh);
set(fbar,'FaceColor','blue','EdgeColor','none');
text(0,0.06,['$\leftarrow (b=0)$ share $=' num2str(round(cshare_b(11),2)) '$'],'FontSize',20,'interpreter','latex','Color','k');
xlim([-30 maxliqhist+2]);
ylim([0 0.1]);
alpha(0.4);

axbar.YAxisLocation = 'right';
axline = axes('Position',axbar.Position,...
    'XAxisLocation','bottom','YAxisLocation','left','Color','none','XLim',xlim);
grid on;
hold on;
fline = plot(bgrid.*dataannoutput/(initss.output*4*1000), sm_elastc_b{8},'Parent',axline, 'Color','k','LineWidth',3.0);
% axline.YLim = ([-0.05 1.15]);
% axline.YLim = ([-1 4]);
ylabel('Elasticity','FontSize',20,'interpreter','latex');
xlabel('\$ Thousands','FontSize',20,'interpreter','latex');
set(axbar,'FontSize',16) ;
set(axline,'FontSize',16) ;
grid on;
ylim([-0.6 7])

plot(bgrid.*datagrosslabinc/(annlabinc*1000),100.*zeros(size(bgrid)),'k-','LineWidth',1.0);
hold off;

if Save==1
    print('-depsc',[SaveFigsDir '/fig_5a']);
end

%% Elasticity C figure: Breakdown into Direct  and Indirect Effects
%4 is rb, so 4 is direct effect
%5 is ra, 3 is W+profit, 7 is tau, 6 is equity, so 3+5+6+7is indirect effect

figure;

[fbar]   = bar(liqbintrunc, liqbinpdftrunc,'histc');
axbar = gca;
sh  = findall(gcf,'marker','*'); delete(sh);
set(fbar,'FaceColor','blue','EdgeColor','none');
text(0,0.06,['$\leftarrow (b=0)$ share $=' num2str(round(cshare_b(11),2)) '$'],'FontSize',20,'interpreter','latex','Color','k');
xlim([-30 maxliqhist+2]);
ylim([0 0.1]);
alpha(0.4);

axbar.YAxisLocation = 'right';
axline = axes('Position',axbar.Position,...
    'XAxisLocation','bottom','YAxisLocation','left','Color','none','XLim',xlim);
grid on;
hold on;
fline = plot(bgrid.*dataannoutput/(initss.output*4*1000), sm_elastc_b{4},'Parent',axline, 'Color','b','LineWidth',3.0);
fline2 = plot(bgrid.*dataannoutput/(initss.output*4*1000), sm_elastc_b{3} + sm_elastc_b{5}+ sm_elastc_b{6}+sm_elastc_b{7},'Parent',axline, 'Color','r','LineWidth',3.0,'LineStyle','--');
% axline.YLim = ([-0.05 1.15]);
% axline.YLim = ([-1 4]);
ylabel('Elasticity','FontSize',20,'interpreter','latex');
xlabel('\$ Thousands','FontSize',20,'interpreter','latex');
set(axbar,'FontSize',16) ;
set(axline,'FontSize',16) ;
grid on;
ylim([-0.6 7])
leg = legend({'Direct Effects' 'Indirect Effects' },'Location','NorthEast','Interpreter','latex');
set(leg,'FontSize',16);

plot(bgrid.*datagrosslabinc/(annlabinc*1000),100.*zeros(size(bgrid)),'k-','LineWidth',1.0);
hold off;

if Save==1
    print('-depsc',[SaveFigsDir '/fig_5b']);
end


%% Elasticity C figure: Breakdown of Direct Effects
%4 is rb, so 4 is direct effect
figure;

[fbar]   = bar(liqbintrunc, liqbinpdftrunc,'histc');
axbar = gca;
sh  = findall(gcf,'marker','*'); delete(sh);
set(fbar,'FaceColor','blue','EdgeColor','none');
text(0,0.06,['$\leftarrow (b=0)$ share $=' num2str(round(cshare_b(11),2)) '$'],'FontSize',20,'interpreter','latex','Color','k');
xlim([-30 maxliqhist+2]);
ylim([0 0.1]);
alpha(0.4);

axbar.YAxisLocation = 'right';
axline = axes('Position',axbar.Position,...
    'XAxisLocation','bottom','YAxisLocation','left','Color','none','XLim',xlim);
grid on;
hold on;
fline = plot(bgrid.*dataannoutput/(initss.output*4*1000), sm_elastc_b{4},'Parent',axline, 'Color','b','LineWidth',3.0);
fline2 = plot(bgrid.*dataannoutput/(initss.output*4*1000), sm_direct_sub_b{4},'Parent',axline, 'Color','r','LineWidth',3.0,'LineStyle','--');
fline4 = plot(bgrid.*dataannoutput/(initss.output*4*1000), sm_direct_wealth_b{4},'Parent',axline, 'Color','g','LineWidth',3.0,'LineStyle','-.');
fline3 = plot(bgrid.*dataannoutput/(initss.output*4*1000), sm_direct_port_b{4},'Parent',axline, 'Color','m','LineWidth',2.5,'LineStyle','-');
ylabel('Elasticity','FontSize',20,'interpreter','latex');
xlabel('\$ Thousands','FontSize',20,'interpreter','latex');
set(axbar,'FontSize',16) ;
set(axline,'FontSize',16) ;
grid on;
ylim([-0.6 7])
leg = legend({'Direct Effects' 'Substitution Effect' 'Income Effect' 'Portfolio Reallocation Contribution'  },'Location','NorthEast','Interpreter','latex');
set(leg,'FontSize',16);

plot(bgrid.*datagrosslabinc/(annlabinc*1000),100.*zeros(size(bgrid)),'k-','LineWidth',1.0);
hold off;

if Save==1
    print('-depsc',[SaveFigsDir '/fig_6a']);
end


%% Elasticity C figure: Breakdown of Indirect Effects
%5 +6 is ra +equity
%3 is W +Pi
%7 is tau
figure;

[fbar]   = bar(liqbintrunc, liqbinpdftrunc,'histc');
axbar = gca;
sh  = findall(gcf,'marker','*'); delete(sh);
set(fbar,'FaceColor','blue','EdgeColor','none');
text(0,0.06,['$\leftarrow (b=0)$ share $=' num2str(round(cshare_b(11),2)) '$'],'FontSize',20,'interpreter','latex','Color','k');
xlim([-30 maxliqhist+2]);
ylim([0 0.1]);
alpha(0.4);

axbar.YAxisLocation = 'right';
axline = axes('Position',axbar.Position,...
    'XAxisLocation','bottom','YAxisLocation','left','Color','none','XLim',xlim);
grid on;
hold on;
fline = plot(bgrid.*dataannoutput/(initss.output*4*1000), sm_elastc_b{5}+sm_elastc_b{6},'Parent',axline, 'Color','m','LineWidth',2.5,'LineStyle','-');
fline2 = plot(bgrid.*dataannoutput/(initss.output*4*1000), sm_elastc_b{3},'Parent',axline, 'Color','c','LineWidth',3.0,'LineStyle','-.');
fline3 = plot(bgrid.*dataannoutput/(initss.output*4*1000), sm_elastc_b{7},'Parent',axline, 'Color','b','LineWidth',3.0,'LineStyle',':');
% axline.YLim = ([-0.05 1.15]);
% axline.YLim = ([-1 4]);
ylabel('Elasticity','FontSize',20,'interpreter','latex');
xlabel('\$ Thousands','FontSize',20,'interpreter','latex');
set(axbar,'FontSize',16) ;
set(axline,'FontSize',16) ;
grid on;
ylim([-0.6 7])
leg = legend({'Indirect Effect: $r^a$' 'Indirect Effect: $\omega$' 'Indirect Effect: $T$' },'Location','NorthEast','Interpreter','latex');
set(leg,'FontSize',16);

plot(bgrid.*datagrosslabinc/(annlabinc*1000),100.*zeros(size(bgrid)),'k-','LineWidth',1.0);
hold off;

if Save==1
    print('-depsc',[SaveFigsDir '/fig_6b']);
end








