function disppar( p )
%DISPPAR displays a struct
%
% 2018, Alexander Heimel

flds = fieldnames(p);
for i=1:length(flds)
    if isnumeric(p.(flds{i})) ||  islogical(p.(flds{i}))
        disp(['par.' flds{i} ' = ' mat2str(p.(flds{i}),3) ';']);
    elseif isa(p.(flds{i}),'function_handle')
        disp(['par.' flds{i} ' = ' func2str(p.(flds{i})) ';']);
    else
        disp(['par.' flds{i} ' = ' eval(p.(flds{i})) ';']);
    end
end
