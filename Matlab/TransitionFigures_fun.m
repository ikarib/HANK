function out = TransitionFigures_fun(options)

out = 1;

%% Directories
BaseDir         = options.BASEDIR;
BaseSaveDir     = [BaseDir,'/','Output'];
ShockType       = options.IRF;                          
StimulusType    = options.STIM;                         

%% Options
Experiment              = options.Experiment;
Save                    = options.Save;
prices                  = options.prices;

%% Create directories for saving and copy base .tex file over to that directory

BaseOutputDir   = BaseDir;
% SaveDir         = [BaseSaveDir,'/',ShockType,'/',StimulusType];
% eval(sprintf('mkdir %s',SaveDir));
% copyfile(sprintf('%s/FIGURES.tex',BaseDir),SaveDir); 

%% Grids
tstep       = load([BaseOutputDir '/deltatransvec.txt']);

%% Initial steady state
temp = importdata([BaseOutputDir '/InitialSteadyStateParameters.txt']);
for i = 1:size(temp.data,1)
    initss.(temp.textdata{i}) = temp.data(i,1);
end

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
%% Other transitions
OutputDir           = [BaseOutputDir '/IRF_' ShockType '/' StimulusType];

%% Load all transitions
% QQ      = 1;
% Qlist   = {'a'};
for j = 1:numel(prices);
    eval(sprintf('Z = options.%s;',prices{j}));
    if (Z==1);
        s                   = dir([OutputDir,'/',upper(prices{j}),'/*.txt']);
        s                   = {s.name};
        for i = 1:numel(s);
            var     = s{i}(1:end-4);    % Takes off .txt which is last 4 characters
            ZZ      = load(sprintf('%s',[OutputDir '/' upper(prices{j}) '/' var '.txt']));
            eval(sprintf('%s.%s = ZZ;',prices{j},var));
%             if size(ZZ,2)==4
%                Qlist{QQ}    = var;
%                QQ           = QQ+1;
%             end
        end 
    end
    if (j>1)&&(Z==0)
         eval(sprintf('%s = %s;',prices{j},'flex'));
%        eval(sprintf('%s = %s;',prices{j},'sticky'));
    end
end

% Qlist = unique(Qlist);

%% Transition path: aggregate variables
% cd(SaveDir);
initss.housefrac = 0;

T                       = length(tstep);
tpoints                 = cumsum(tstep);
% tlim                    = [0 min(tmaxplot,max(tpoints))];
ylimits                 = [0.95 1.05];

initss.mpshock          = 0;
initss.markup           = 1./initss.mc-1;
initss.hhwealth         = initss.Ea+initss.Eb;
initss.hhwliqpfshare    = initss.Eb./(initss.Ea+initss.Eb);
initss.SHAREc_incQ      = initss.Ec_incQ./(sum(initss.Ec_incQ,2)*ones(1,4));
initss.SHAREc_nwQ       = initss.Ec_nwQ./(sum(initss.Ec_nwQ,2)*ones(1,4));
initss.SHAREinc_incQ    = initss.Einc_incQ./(sum(initss.Einc_incQ,2)*ones(1,4));
initss.SHAREinc_nwQ     = initss.Einc_nwQ./(sum(initss.Einc_nwQ,2)*ones(1,4));

initss.inv_hh           = ((1-initss.housefrac)./(1-initss.fundlev)).* initss.deprec.*initss.Ea;
initss.netexports       = initss.rb .* initss.worldbond;

for j = 1:numel(prices);
    eval(sprintf('%s.markup = 1./%s.mc-1;',prices{j},prices{j}));
    eval(sprintf('%s.hhwealth = %s.Ea+%s.Eb;',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.hhwliqpfshare = %s.Eb./(%s.Ea+%s.Eb);',prices{j},prices{j},prices{j},prices{j}));
    eval(sprintf('%s.SHAREc_incQ = %s.Ec_incQ./(sum(%s.Ec_incQ,2)*ones(1,4));',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.SHAREc_nwQ = %s.Ec_nwQ./(sum(%s.Ec_nwQ,2)*ones(1,4));',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.SHAREinc_incQ = %s.Einc_incQ./(sum(%s.Einc_incQ,2)*ones(1,4));',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.SHAREinc_nwQ = %s.Einc_nwQ./(sum(%s.Einc_nwQ,2)*ones(1,4));',prices{j},prices{j},prices{j}));
    eval(sprintf('%s.inv_hh = ([%s.Ea(2:T); initss.Ea] - %s.Ea)./tstep  + initss.deprec.*%s.Ea;',prices{j},prices{j},prices{j},prices{j}));
    eval(sprintf('%s.inv_hh = ((1-initss.housefrac)./(1-initss.fundlev)).* %s.inv_hh;',prices{j},prices{j}));
    eval(sprintf('%s.netexports = -([%s.worldbond(2:T); initss.worldbond] - %s.worldbond)./tstep  + %s.rb.*%s.worldbond;',prices{j},prices{j},prices{j},prices{j},prices{j}));
end





%% Save the workspace
cd(BaseDir);
save(['IRF_' ShockType '_' StimulusType '_workspace.mat']);


