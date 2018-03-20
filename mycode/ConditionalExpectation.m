function varargout = ConditionalExpectation(n,x,f,dx,p,varargin)
% finds y = E[g(x) | xL<x<=xH] using a trapezoidal rule at boundaries
% inputs: x,f,dx,g - vectors of size n
%         xL,xH - vectors of size m
% outputs: y - vector of same size as g in varargin

p(p<x(1)) = x(1);
p(p>x(n)) = x(n);
xp = [ x(1) p    % xL
       p x(n) ]; % xH
fdx = f.*dx;
df = diff(interp1(x,cumsum(fdx),xp));
for i=1:nargin-5
    varargout{i} = diff(interp1(x,cumsum(varargin{i}.*fdx),xp))./df; %#ok<AGROW>
end