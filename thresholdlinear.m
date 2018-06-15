function y = thresholdlinear(x)
%THRESHOLDLINEAR returns x if x>0, otherwise 0
%
% shifted to match sigmoid(0) and sigmoid'(0)
%
% 2018, Alexander Heimel

y = x/4+1/2;
y(y<0) = 0;
