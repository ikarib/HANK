clear;
close all;

%% PARAMETERS

BaseInputDir = '~/FortranOutputDir/';
lSaveDir =  '~/FiguresDir'; 

%put the paths to the subdirectories for each vale of rho here
modeldir{1} = 'rho_value_1'; 
modeldir{2} = 'rho_value_2';
modeldir{3} = 'etc...';
modeldir{4} = 'etc...';
modeldir{5} = 'etc...';
modeldir{6} = 'etc...';
modeldir{7} = 'etc...';
modeldir{8} = 'etc...';
modeldir{9} = 'etc...';
modeldir{10} = 'etc...';

nmodel = 10;

Save = 1;

%% LOAD DATA
rhovec = zeros(nmodel,1);
BYvec = zeros(nmodel,1);
MPCvec = zeros(nmodel,1);
Celast = zeros(nmodel,1);
Celast_partial = zeros(nmodel,1);

for im = 1:nmodel
    InputDir = [BaseInputDir modeldir{im}];
    NOFS{im} = load([InputDir '/IRF_Monetary_NOFS_workspace.mat']);
    PE4{im} = load([InputDir '/IRF_Monetary_PE4_workspace.mat']);
    PE9{im} = load([InputDir '/IRF_Monetary_PE9_workspace.mat']);
    SS{im} = load([InputDir '/Steadystate_workspace.mat']);
    rhovec(im,1) = SS{im}.initss.rho;
    BYvec(im,1) = SS{im}.initss.Eb ./ (4.*SS{im}.initss.output);
    MPCvec(im,1) = SS{im}.Empreb1;
    
    
 
    tstep = NOFS{im}.tstep;
     tset = [1:12];
     tsetRb = [2:13];
    
    elastdenom = sum(NOFS{im}.sticky.rb(tsetRb).*tstep(tsetRb))./ sum(tstep(tsetRb))- SS{im}.initss.rb;
    %total consumption elasticity;
    dC =  NOFS{im}.sticky.Ec(tset) - SS{im}.initss.Ec;
    Celast(im,1) = - ( sum(dC.*tstep(tset)./SS{im}.initss.Ec)./ sum(tstep(tset))) ./ elastdenom;

    %partial consumption elasticity;
    dC4 =  PE4{im}.sticky.Ec(tset) - SS{im}.initss.Ec;
    Celast_partial(im,1) = - ( sum(dC4.*tstep(tset)./SS{im}.initss.Ec)./ sum(tstep(tset))) ./ elastdenom;

end

SaveDir = lSaveDir;

%% PLOT: MPC, B/Y
figure;
[ax h1 h2]=plotyy(rhovec.*400, BYvec,rhovec.*400, MPCvec); 
h1.LineWidth = 2.5;
h2.LineWidth = 2.5;
h1.MarkerSize= 10;
h2.MarkerSize= 10;
h1.Marker= 'o';
h2.Marker= '^';

xlim(ax(1),[2 7.5]);
xlim(ax(2),[2 7.5]);
set(ax,'FontSize',16) ;
grid on;
xlabel(ax(1),'Annual Discount Rate (\%)','FontSize',20,'interpreter','latex');
ylabel(ax(1),'Mean Wealth (relative to annual GDP) ','FontSize',20,'interpreter','latex');
ylabel(ax(2),'Quarterly MPC \$500','FontSize',20,'interpreter','latex');
if Save==1
    print('-depsc',[SaveDir '/fig_7a']);
end

%% PLOT: C elast
figure;
plot(rhovec.*400, Celast,'b-o',rhovec.*400, Celast_partial,'r-.s','MarkerSize',10,'LineWidth',2.5); 
set(gca,'FontSize',16) ;
ylim([-1 4.5])
xlim([2 7.5]);

grid on;
xlabel('Annual Discount Rate (\%)','FontSize',20,'interpreter','latex');
ylabel('Elasticity ','FontSize',20,'interpreter','latex');
legend({'Total Elasticity' 'Direct Elasticity'},'Location','Best','Interpreter','latex','FontSize',16);

if Save==1
    print('-depsc',[SaveDir '/fig_7b']);
end

