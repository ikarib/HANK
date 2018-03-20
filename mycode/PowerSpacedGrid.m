function y=PowerSpacedGrid(n,k,low,high)
% gives a grid of n points spaced between low and high based on the unit interval with a function x^(1/k)
% k = 1 is linear, k = 0 is L-shaped

if n<2
    error('n must be at least 2 to make grids')
end
if n==2
    y = [low high];
    return
end

y = low + (high-low)*linspace(0,1,n).^(1/k);
