function y = nakarushton(x,n)
%NAKARUSHTON computes nakarushton with parameters chosen to match sigmoid(0) and sigmoid'(0)
%
% Y = NAKARUSHTON(X,N)
%
% 2018, Alexander Heimel

if nargin<2 || isempty(n)
    n = 2;
end
    
a = -2;
b = 2;
y = (x-a).^n./(b^2+(x-a).^n);