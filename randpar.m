function par = randpar( par,scale)
%RANDPAR randomly perturbs the values in a struct
%
% PAR = RANDPAR( PAR, SCALE)
%
% 2018, Alexander Heimel

if nargin<2
    scale = 0.1;
end

flds = fieldnames(par);
for i=1:length(flds)
    par.(flds{i}) = (1+scale*(rand(1)-0.5))*par.(flds{i});
end
