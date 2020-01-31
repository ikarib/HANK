function f = dfovec(x)
% tic
global guess Estimate nsim dt Tburn Tsim Tann UseNormalDist AdditiveDrift Match Target xmax yjumprand yrand yannlevsim
% objective function for DFLS minimization
% output f vector of least squares objective

% n = numel(x);
% lp(12),lweight(m),lf
% im,ip

% fprintf('*********************************');
% fprintf('EVALUATION NUMBER: %d',objeval);

% extract parameter vector
guess(Estimate,:) = xmax./(1+exp(-x));
lambda = guess(1,:);
zetaP =  guess(2,:);
zetaN = -guess(3,:);
rho = guess(4,:);
sigma = guess(5,:);
delta = guess(6,:);

% fprintf(' lambda: %f  %f\n',lambda);
% fprintf(' zetaP:  %f  %f\n',zetaP);
% fprintf(' zetaN:  %f  %f\n',zetaN);
% fprintf(' rho:    %f  %f\n',rho);
% fprintf(' sigma:  %f  %f\n',sigma);
% fprintf(' delta:  %f  %f\n',delta);

Simulate
ComputeMoments

Model = [mu2y mu2dy1 gam3dy1 gam4dy1 mu2dy5 gam3dy5 gam4dy5 fracdy1less5 fracdy1less10 fracdy1less20 fracdy1less50]';
weight = [1 1 1 .5 1 .5 .5 1 1 1 1]';
f = (Model(Match)./Target(Match)-1).*weight(Match);
% moments = {'VarLogY','VarD1LogY','SkewD1LogY','KurtD1LogY','VarD5LogY','SkewD5LogY','KurtD5LogY','FracD1Less5','FracD1Less10','FracD1Less20','FracD1Less50'};
% fprintf(' Moment:\tTarget:\tModel:\n')
% for i=1:numel(moments)
%     if Match(i)
%         fprintf(' %s\t%g\t%g\n',moments{i},Target(i),Model(i))
%     end
% end
obj = sqrt(sum(f.^2)/sum(weight(Match)));
fprintf(' objective fun: %e\n',obj)
% [muy mu2y gam3y gam4y-3; mudy1 mu2dy1 gam3dy1 gam4dy1-3; mudy5 mu2dy5 gam3dy5 gam4dy5-3]
% toc
% [obj 7*(obj)^2/2]