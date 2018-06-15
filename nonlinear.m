function y = nonlinear(x,n)
%NONLINEAR is a threshold-nonlinear function
%
%  Y = NONLINEAR(X,N)
%
%   parameters chosen to match sigmoid(0)=1/2 and sigmoid'(0)
%
%
% 2018, Alexander Heimel

if nargin<2
    n = 2;
    a = 1/32;
    t = -4;
else
    a = 1/(2*(2*n)^n);
    t = -2*n;
end

x = x-t;
x = x.*(1+sign(x))/2;% equiv to x(x<0) = 0 but faster
y = a*(x.^n);

