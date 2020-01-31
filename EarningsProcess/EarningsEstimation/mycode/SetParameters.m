global guess Estimate nsim dt Tburn Tsim Tann UseNormalDist AdditiveDrift Match Target xmax yjumprand yrand yannlevsim

nsim = 4992; %5000; %50000
dt = 0.25; %0.005  %time step in quarters 
Tburn = floor(100/dt)+1;
Tsim = Tburn + floor(20/dt)+1; % only need 20 quarters
Tann = floor((Tsim-Tburn)*dt/4);

Display                = false;
SaveSimulations        = true;
RevertToMedianWithSkew = false;
UseNormalDist          = true;
AdditiveDrift          = false;


Estimate = [ true  % Lambda
             false % ZetaP
             false % ZetaN
             false % Rho
             true  % Sigma
             true  % Delta
           ];

Match = [
    true  % VarLogY
    true  % VarD1LogY
    false % SkewD1LogY
    true  % KurtD1LogY
    true  % VarD5LogY
    false % SkewD5LogY
    true  % KurtD5LogY
    false % FracD1Less5
    true  % FracD1Less10
    true  % FracD1Less20
    true  % FracD1Less50
  ];

Target = [ % USA
    0.7  % VarLogY
    0.23 % VarD1LogY
   -1.35 % SkewD1LogY
   17.8  % KurtD1LogY
    0.46 % VarD5LogY
   -1.01 % SkewD5LogY
   11.55 % KurtD5LogY
    0.35 % FracD1Less5
    0.54 % FracD1Less10
    0.71 % FracD1Less20
    0.86 % FracD1Less50
  ];
Target = [ % Canada
    0.76  % VarLogY
    0.217 % VarD1LogY
   nan % SkewD1LogY
   13.377  % KurtD1LogY
    0.437 % VarD5LogY
   nan % SkewD5LogY
   8.782 % KurtD5LogY
    nan % FracD1Less5
    0.51 % FracD1Less10
    0.68 % FracD1Less20
    0.85 % FracD1Less50
  ];

bounds = [
    0 2 % lambda, jump process cant arrive more frequently than quarterly on average
    1 1000 % zetaP
    1 1000 % zetaN
    0 1 % rho
    0 2 % sigma
    0 1 % deltaUseNormalDist
  ];

OutputDir = 'earnings_estimation_output/';
if ~exist(OutputDir,'dir'); mkdir(OutputDir); end

if 0
    % set random seed
    rng default
    yjumprand = rand(nsim,2,Tsim);
    if UseNormalDist
        yrand = randn(nsim,2,Tsim);
    else
        p = 2*rand(Tsim,2,nsim); % DoubleParetoInverseCDF
        yrand = log(p);
        yrand(p>1) = -log(2-p(p>1));
    end
else
    if exist('fortran.mat','file')
        load fortran
    else
        load('../earnings_estimation_output/yjumprand.txt'); yjumprand=permute(reshape(yjumprand,nsim,Tsim,2),[1 3 2]);
        load('../earnings_estimation_output/yrand.txt'); yrand=permute(reshape(yrand,nsim,Tsim,2),[1 3 2]);
        save fortran yjumprand yrand
    end
end