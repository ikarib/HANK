clear;
close all;

%% PARAMETERS

BaseInputDir = '~/FortranOutputDir/'; 
lSaveDir =  '~/FiguresDir'; 

modeldir{1} = 'philcurve_T'; %folders should be named modeldir{1}_i for i = 1 to nphilpoints;
modeldir{2} = 'philcurve_B';

nmodel = 2;
nphilpoints = 14;

SaveDir = lSaveDir;
Save = 1;



%% LOAD DATA
%rank with capital;
temp = importdata([SaveDir '/rank_capital_tradeoff.txt']);
for i = 1:size(temp.data,2)
    rank_capital.(temp.textdata{i}) = temp.data(:,i);
end

%initial steady state;
temp = importdata([BaseInputDir modeldir{1} '_' int2str(1) '/InitialSteadyStateParameters.txt']);
for i = 1:size(temp.data,1)
    initss.(temp.textdata{i}) = temp.data(i,1);
end

% timestep
tstep   = load([BaseInputDir modeldir{1} '_' int2str(1) '/deltatransvec.txt']);
% tset = [1:12];
% tsetRb = [2:13];
tset = [1:3];
tsetRb = [2:4];


for im = 1:nmodel
    
    model{im}.output = zeros(nphilpoints,1);
    model{im}.pi = zeros(nphilpoints,1);
    model{im}.rb = zeros(nphilpoints,1);
    model{im}.rnom = zeros(nphilpoints,1);
    model{im}.con = zeros(nphilpoints,1);

    for ip = 1:nphilpoints
            
        temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/output.txt']);
        model{im}.output(ip) = sum(temp(tset).*tstep(tset))./sum(tstep(tset))./ initss.output - 1;
        
        temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/pi.txt']);
        model{im}.pi(ip) = sum(temp(tset).*tstep(tset))./sum(tstep(tset)) -  initss.pi;
        
        temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/rb.txt']);
        model{im}.rb(ip) = sum(temp(tsetRb).*tstep(tsetRb))./sum(tstep(tsetRb)) - initss.rb;
        
        temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/rnom.txt']);
        model{im}.rnom(ip) = sum(temp(tsetRb).*tstep(tsetRb))./sum(tstep(tsetRb)) - initss.rnom;

        temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/mc.txt']);
        model{im}.mc(ip) = sum(temp(tset).*tstep(tset))./sum(tstep(tset))./ initss.mc- 1;

%         temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/Ec.txt']);
%         model{im}.con(ip) = temp(1);
%     
%         temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/investment.txt']);
%         model{im}.inv(ip) = temp(1);
% 
        temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/wage.txt']);
        model{im}.wage(ip) = sum(temp(tset).*tstep(tset))./sum(tstep(tset))./ initss.wage - 1;
        
        temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/labor.txt']);
        model{im}.labor(ip) = sum(temp(tset).*tstep(tset))./sum(tstep(tset))./ initss.labor - 1;
        
        % 
%         temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/mc.txt']);
%         model{im}.mc(ip) = temp(1);
% 
%         temp = load([BaseInputDir modeldir{im} '_' int2str(ip) '/IRF_Monetary/NOFS/STICKY/rcapital.txt']);
%         model{im}.rcapital(ip) = temp(1);

    end
end


%% RANK MODEL
shockvec = [-0.0125; -0.01; -0.0075; -0.005; -0.00375; -0.0025; -0.00125;...
             0.00125; 0.0025; 0.00375; 0.005; 0.0075; 0.01; 0.0125];

% rank.rho = initss.rho;
rank.rho = initss.rb;
rank.eta = 0.5;
rank.phi = 1.25;
rank.gam = 1;
rank.eps = 1/(1-initss.mc);
rank.frisch = 1;
rank.chi = (1/3).^(-(1+rank.frisch)./rank.frisch);
rank.theta = 100;
rank.alpha = 0;
rank.kappa = (rank.eps-1).*(1./rank.frisch + rank.gam)./rank.theta;
rank.kappa = rank.kappa .* (1-rank.alpha)/(1- rank.alpha+ rank.alpha*rank.eps);

rank.denom = (rank.phi-1).*rank.kappa + rank.gam.*rank.eta.*(rank.rho + rank.eta);
rank.outputgap = -(rank.rho+rank.eta).*shockvec./rank.denom;
rank.outputgap = rank.outputgap .* (1-rank.alpha)/(1- rank.alpha+ rank.alpha*rank.eps);

rank.inflation = -rank.kappa.*shockvec./rank.denom;
rank.rb = rank.rho + rank.gam.*rank.eta*(rank.rho+rank.eta).*shockvec./rank.denom;
rank.rnom = rank.rb + rank.inflation;
rank.mc = rank.outputgap.*(rank.gam + 1/rank.frisch).*(rank.eps-1)./rank.eps;
rank.mc = rank.mc .* (1-rank.alpha)/(1- rank.alpha+ rank.alpha*rank.eps);


%% PLOT: INFLATION vs OUTPUT GAP, SHORTER RANGE
figure;
hold on;
plot(model{1}.output(4:11).*100,model{1}.pi(4:11)*4*100,'-ok','LineWidth',2.5,'MarkerSize',8);
plot(model{3}.output(4:11).*100,model{3}.pi(4:11)*4*100,'->b','LineWidth',2.5,'MarkerSize',8);
plot(rank_capital.outputgap(4:11),rank_capital.inflation(4:11),'c--s','LineWidth',2.5,'MarkerSize',8,'MarkerFace','c');

ylim([-2.5 2.5]);
xlim([-1.5 1.5]);

hold off;
grid on;
title('Taylor Rule Shocks $\in [-2\%,+2\%]$, \% p.a.','FontSize',20,'interpreter','latex');
ylabel('Inflation, \% p.a.','FontSize',20,'interpreter','latex');
xlabel('Output Gap, \%','FontSize',20,'interpreter','latex');
set(gca,'FontSize',16) ;
leg = legend({'$T$ adjusts' '$B^g$ adjusts'  'RANK model' },'Location','NorthWest','Interpreter','latex');
set(leg,'FontSize',16);

if Save==1
    print('-depsc',[SaveDir '/fig_9a']);
end

%% PLOT: INFLATION vs MARGINAL COSTS, SHORTER RANGE
figure;
hold on;
plot(model{1}.mc(4:11).*100,model{1}.pi(4:11)*4*100,'-ok','LineWidth',2.5,'MarkerSize',8);
plot(model{2}.mc(4:11).*100,model{2}.pi(4:11)*4*100,'->b','LineWidth',2.5,'MarkerSize',8);
plot(rank_capital.mc(4:11),rank_capital.inflation(4:11),'c--s','LineWidth',2.5,'MarkerSize',8,'MarkerFace','c');

ylim([-2.5 2.5]);
xlim([-4 4]);
hold off;
grid on;
title('Taylor Rule Shocks $\in [-2\%,+2\%]$, \% p.a.','FontSize',20,'interpreter','latex');
ylabel('Inflation, \% p.a.','FontSize',20,'interpreter','latex');
xlabel('Marginal Costs, \% dev','FontSize',20,'interpreter','latex');
set(gca,'FontSize',16) ;
leg = legend({'$T$ adjusts' '$B^g$ adjusts'  'RANK model' },'Location','SouthEast','Interpreter','latex');
set(leg,'FontSize',16);

if Save==1
    print('-depsc',[SaveDir '/fig_9b']);
end

%% PLOT: MARGINAL COST vs OUTPUT GAP, SHORTER RANGE
figure;
hold on;
plot(model{1}.output(4:11).*100,model{1}.mc(4:11)*100,'-ok','LineWidth',2.5,'MarkerSize',8);
plot(model{2}.output(4:11).*100,model{2}.mc(4:11)*100,'->b','LineWidth',2.5,'MarkerSize',8);
plot(rank_capital.outputgap(4:11),rank_capital.mc(4:11),'c--s','LineWidth',2.5,'MarkerSize',8,'MarkerFace','c');

ylim([-4 4]);
xlim([-1.5 1.5]);

hold off;
grid on;
title('Taylor Rule Shocks $\in [-2\%,+2\%]$, \% p.a.','FontSize',20,'interpreter','latex');
ylabel('Marginal Cost, \% dev.','FontSize',20,'interpreter','latex');
xlabel('Output Gap, \%','FontSize',20,'interpreter','latex');
set(gca,'FontSize',16) ;
leg = legend({'$T$ adjusts' '$B^g$ adjusts' 'RANK model' },'Location','NorthWest','Interpreter','latex');
set(leg,'FontSize',16);

if Save==1
    print('-depsc',[SaveDir '/fig_9c']);
end

