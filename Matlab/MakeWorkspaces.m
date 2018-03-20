clear all;clc;

%% Parameters

%path to folder containing output directories
BaseOutputDir = '../FortranOutputDir/';

%list of output directories to process
explist{1} = 'BaselineOutputSubdir';

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NO INPUT REQUIRED AFTER THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IRF                 = {'Monetary'};
STIM                = {'NOFS','PE1','PE2','PE3','PE4','PE5','PE6','PE7','PE8','PE9','PE10','PE11','PE12','PE13','PE14','PE15'};
PRICES              = {'sticky'};  

options.OnlyWorkspace           = 1;  %for the PE experiments does not make figures
options.ipoint                  = 17; %17;
options.Quartiles               = 'Y';
options.Save                    = 1;
options.LoadDistributions       = 1;
options.datagrosslabinc         = 69100;
options.dataannoutput           = 115000;

HOMEDIR = pwd;

for ie = 1:numel(explist)
    options.Experiment              =  explist{ie}; 
    

    %% Set up directories
    copyfile([HOMEDIR '/SteadyStateFigures_fun.m'],[BaseOutputDir,'/',options.Experiment],'f');
    copyfile([HOMEDIR '/TransitionFigures_fun.m'],[BaseOutputDir,'/',options.Experiment],'f');
    
    cd([BaseOutputDir '/' options.Experiment]); 
    BASEDIR                         = pwd;
    options.BASEDIR                 = BASEDIR;

    %% Steady state
    disp('Computing steady state');
                
    SteadyStateFigures_fun(options);

    %% Transitions 

    options.prices = PRICES;
    for i = 1:1:numel(IRF)
        for j = 1:1:numel(STIM)
            irfdir  = [BASEDIR,'/','IRF_',IRF{i},'/',STIM{j}];

            if (exist(irfdir,'dir')==7)
                disp(['Computing - ' IRF{i} '/' STIM{j}]);
    
                options.IRF                     = IRF{i};
                options.STIM                    = STIM{j};

                % Set options for types of price transitions
                for k = 1:numel(options.prices);
                    Z = (exist([irfdir,'/',upper(options.prices{k})],'dir')==7);
                    eval(sprintf('options.%s = Z;',options.prices{k}));
                end

                % Compute transition figures
                TransitionFigures_fun(options); 
            end
        end
    end

    disp(' ');

    %% Back to base directory
    cd ..
end            
         
