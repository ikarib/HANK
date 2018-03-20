clc;
clear;
close all;

%% DIRECTORY PATHS
lSaveDir =  '~/FiguresDir'; 

% T adjusts
InputDirT{1} = '~/FortranOutputDir/T_pers05';
InputDirT{2} = '~/FortranOutputDir/T_pers10';
InputDirT{3} = '~/FortranOutputDir/T_pers15';
InputDirT{4} = '~/FortranOutputDir/T_pers20';
InputDirT{5} = '~/FortranOutputDir/T_pers25';
InputDirT{6} = '~/FortranOutputDir/T_pers50';
InputDirT{7} = '~/FortranOutputDir/T_pers75';
InputDirT{8} = '~/FortranOutputDir/T_pers80';
InputDirT{9} = '~/FortranOutputDir/T_pers85';
InputDirT{10} = '~/FortranOutputDir/T_pers90';
InputDirT{11} = '~/FortranOutputDir/T_pers95';


% B adjusts
InputDirB{1} = '~/FortranOutputDir/B_pers05';
InputDirB{2} = '~/FortranOutputDir/B_pers10';
InputDirB{3} = '~/FortranOutputDir/B_pers15';
InputDirB{4} = '~/FortranOutputDir/B_pers20';
InputDirB{5} = '~/FortranOutputDir/B_pers25';
InputDirB{6} = '~/FortranOutputDir/B_pers50';
InputDirB{7} = '~/FortranOutputDir/B_pers75';
InputDirB{8} = '~/FortranOutputDir/B_pers80';
InputDirB{9} = '~/FortranOutputDir/B_pers85';
InputDirB{10} = '~/FortranOutputDir/B_pers90';
InputDirB{11} = '~/FortranOutputDir/B_pers95';


%% LOAD WORKSPACES

ndirs = 11;
Save = 1;

for j = 1:ndirs
    IRF{j,1} = load([InputDirT{j} '/IRF_Monetary_NOFS_workspace.mat']); %T adjusts
    IRF{j,2} = load([InputDirB{j} '/IRF_Monetary_NOFS_workspace.mat']); %B adjusts
end

SaveDir = lSaveDir;

tstep   = load([InputDirT{1} '/deltatransvec.txt']);
tpoints = IRF{1,1}.tpoints-tstep(1);

Lambda = 0.3; %0.17;
Lambda_T = 0.3;
gamma = 1;
BY_ratio = 1; %0.26;


tlim    = [0 16];
return
%% FIGURES


%% integral under:

for j = 1:11
    for gov = 1:2
       tot_rb_dev{gov}(j) = sum((IRF{j,gov}.sticky.rb-IRF{j,gov}.initss.rb).*tstep);    
       cumsum_rb_dev{gov}(:,j)= cumsum((IRF{j,gov}.sticky.rb-IRF{j,gov}.initss.rb).*tstep,'reverse');
       fourq_rb_dev{gov}(j) = sum(cumsum_rb_dev{gov}(1:12,j).*tstep(1:12));
       init_c_dev{gov}(j) = log(IRF{j,gov}.sticky.Ec(1)./IRF{j,gov}.initss.Ec);
       fourq_c_dev{gov}(j) = sum(log(IRF{j,gov}.sticky.Ec(1:12)./IRF{j,gov}.initss.Ec).*tstep(1:12));
    end
end

eta = zeros(11,1);
eta(1) = 0.05;
eta(2) = 0.10;
eta(3) = 0.15;
eta(4) = 0.20;
eta(5) = 0.25;
eta(6) = 0.50;
eta(7) = 0.75;
eta(8) = 0.80;
eta(9) = 0.85;
eta(10) = 0.90;
eta(11) = 0.95;

%% persistence plots for paper with TANK, separately for T adjusts and B adjusts
etamaxplot = 11;
etaminplot = 2;
neta = etamaxplot - etaminplot+1;

figure;
hold on;
plot(exp(-eta(etaminplot:etamaxplot)),ones(neta,1)./gamma,'k','Linewidth',3.5);
plot(exp(-eta(etaminplot:etamaxplot)),Lambda_T/(1-Lambda).*BY_ratio.*eta(etaminplot:etamaxplot) + 1/gamma,'b--','Linewidth',3.5);
plot(exp(-eta(etaminplot:etamaxplot)),-init_c_dev{1}(etaminplot:etamaxplot)'./tot_rb_dev{1}(etaminplot:etamaxplot)','r-o','Linewidth',3.5,'MarkerSize',12);
hold off;
grid on;
xlim([0.35,1]);
ylim([0.1,1.5]);
xlabel('Quarterly Autocorrelation of Monetary Shock', 'FontSize', 30, 'Interpreter','Latex');
ylabel('(Abs) Cum Elasticity of $C_0$', 'FontSize', 30, 'Interpreter','Latex');
l = legend('RANK','TANK','HANK: $T$ adjusts','Location','Best');
set(gca,'FontSize',24);
set(l,'FontSize',16,'Interpreter','Latex');
if Save==1
    print('-depsc',[SaveDir '/fig_8a']);
end

figure;
hold on;
plot(exp(-eta(etaminplot:etamaxplot)),ones(neta,1)./gamma,'k','Linewidth',3.5);
plot(exp(-eta(etaminplot:etamaxplot)),ones(neta,1)/gamma,'b--','Linewidth',3.5);
plot(exp(-eta(etaminplot:etamaxplot)),-init_c_dev{2}(etaminplot:etamaxplot)'./tot_rb_dev{2}(etaminplot:etamaxplot)','m-s','Linewidth',3.5,'MarkerSize',12);
hold off;
grid on;
xlim([0.35,1]);
ylim([0.1,1.5]);
xlabel('Quarterly Autocorrelation of Monetary Shock', 'FontSize', 30, 'Interpreter','Latex');
ylabel('(Abs) Cum Elasticity of $C_0$', 'FontSize', 30, 'Interpreter','Latex');
l = legend('RANK','TANK','HANK: $B^g$ adjusts','Location','Best');
set(gca,'FontSize',24);
set(l,'FontSize',16,'Interpreter','Latex');
if Save==1
    print('-depsc',[SaveDir '/fig_8b']);
end

