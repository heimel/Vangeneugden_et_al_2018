%COMPARE_NONLINEARITIES
%
% script to compare the various nonlinear transfer functions
%
% 2018, Alexander Heimel

c = linspace(-2,2,100);

figure
hold on

plot(c,sigmoid(c));
plot(c,nonlinear(c));
plot(c,thresholdlinear(c));

plot(c, 1./(-c+2))
n=1;
a=-2;
b=2;
plot(c, (c-a).^n./(b^2+(c-a).^n))
%plot(c, c.^2./(1+c.^2))

legend('sig','non','thr','hyp','nk','location','northwest');
legend boxoff
ylim([0 2]);