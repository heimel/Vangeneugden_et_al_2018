%TRACK_SOLUTIONS 
%
%  tracks solutions from local-integrating, wide-projecting to
%  wide-integrating local-projecting inhibition
%
% 2018, Alexander Heimel

if ~exist('fullpar','var') || isempty(fullpar)
    w_is = 0:0.05:1;
    
    prev_fullpar = model_parsets(0);
    prev_par = model_parsets(10);
    err = [];
    par = prev_par([]);
    fullpar = prev_fullpar([]);
    for i=1:length(w_is)
        logmsg(['Evaluating w_is = ' num2str( w_is(i))]);
        prev_fullpar.w_is = w_is(i);
        [par(i),fullpar(i),err(i)] = feedback_model(prev_par,prev_fullpar,true,false); %#ok<SAGROW>
    end
    
    save('parameters_errors.mat','par','fullpar','err');
end

load('parameters_errors.mat','par','fullpar','err');

rows = 3;
cols = 2;

figure;
subplot(rows,cols,1)
[~,ind] = sort([fullpar.w_is]);
plot([fullpar(ind).w_is],[fullpar(ind).w_et])
axis square 
box off
xlabel('W_{ISE}');
ylabel('W_{ESI}');

subplot(rows,cols,2)
[~,ind] = sort([fullpar.w_is]);
plot([fullpar(ind).w_is],[fullpar(ind).w_es])
axis square 
box off
xlabel('W_{ISE}');
ylabel('W_{ESE}');

subplot(rows,cols,3)
[~,ind] = sort([fullpar.w_is]);
plot([fullpar(ind).w_is],[fullpar(ind).w_ie])
axis square 
box off
xlabel('W_{ISE}');
ylabel('W_{IE}');

subplot(rows,cols,4)
[~,ind] = sort([fullpar.w_is]);
plot([fullpar(ind).w_is],[fullpar(ind).w_ei])
axis square 
box off
xlabel('W_{ISE}');
ylabel('W_{EI}');

subplot(rows,cols,5)
[~,ind] = sort([fullpar.w_is]);
plot([fullpar(ind).w_is],[fullpar(ind).ge])
axis square 
box off
xlabel('W_{ISE}');
ylabel('g_{E}');

subplot(rows,cols,6)
[~,ind] = sort([fullpar.w_is]);
plot([fullpar(ind).w_is],[fullpar(ind).gi])
axis square 
box off
xlabel('W_{ISE}');
ylabel('g_{I}');

