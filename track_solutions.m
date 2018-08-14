%TRACK_SOLUTIONS 
%
%  tracks solutions from local-integrating, wide-projecting to
%  wide-integrating local-projecting inhibition
%
% 2018, Alexander Heimel

if ~exist('fullpar','var') || isempty(fullpar)
    w_ise = 0:0.025:1;
    
    prev_fullpar = model_parsets(0);
    prev_par = model_parsets(10);
     err = [];
    par = prev_par([]);
    fullpar = prev_fullpar([]);
    for i=1:length(w_ise)
        logmsg(['Evaluating w_ise = ' num2str( w_ise(i))]);
        prev_fullpar.w_ise = w_ise(i);
        prev_par.w_esi = 0.62 - 0.32 * w_ise(i);
        
        
        [par(i),fullpar(i),err(i)] = feedback_model(prev_par,prev_fullpar,true,false); %#ok<SAGROW>
        
        prev_par = par(i);
        prev_full = fullpar(i);
    end
    
    save('parameters_errors.mat','par','fullpar','err');
end

load('parameters_errors.mat','par','fullpar','err');

% only take converged results
ind = find(err<0.0001); 
 fullpar = fullpar(ind);


rows = 3;
cols = 2;

yl = [0 2];

figure;
subplot(rows,cols,1)
[~,ind] = sort([fullpar.w_ise]);
plot([fullpar(ind).w_ise],[fullpar(ind).w_esi])
axis square 
box off
xlabel('W_{ISE}');
ylabel('W_{ESI}');
ylim(yl);
set(gca,'xtick',0:0.2:1);
set(gca,'ytick',0:0.4:2);

subplot(rows,cols,2)
[~,ind] = sort([fullpar.w_ise]);
plot([fullpar(ind).w_ise],[fullpar(ind).w_ese])
axis square 
box off
xlabel('W_{ISE}');
ylabel('W_{ESE}');
ylim(yl);
set(gca,'xtick',0:0.2:1);
set(gca,'ytick',0:0.4:2);

subplot(rows,cols,3)
[~,ind] = sort([fullpar.w_ise]);
yyaxis left
plot([fullpar(ind).w_ise],[fullpar(ind).w_ie])
axis square 
box off
xlabel('W_{ISE}');
ylabel('W_{IE}');
hold on 
ylim(yl);
set(gca,'xtick',0:0.2:1);
set(gca,'ytick',0:0.4:2);

yyaxis right
ylim(yl);
set(gca,'ytick',0:0.4:2);
ylabel('W_{EI}');
plot([fullpar(ind).w_ise],[fullpar(ind).w_ei])


subplot(rows,cols,4)
[~,ind] = sort([fullpar.w_ise]);
plot([fullpar(ind).w_ise],[fullpar(ind).w_ei])
axis square 
box off
xlabel('W_{ISE}');
ylabel('W_{EI}');
ylim(yl);
set(gca,'xtick',0:0.2:1);
set(gca,'ytick',0:0.4:2);

subplot(rows,cols,5)
[~,ind] = sort([fullpar.w_ise]);
plot([fullpar(ind).w_ise],[fullpar(ind).ge])
axis square 
box off
xlabel('W_{ISE}');
ylabel('g_{E}');
set(gca,'xtick',0:0.2:1);

subplot(rows,cols,6)
[~,ind] = sort([fullpar.w_ise]);
plot([fullpar(ind).w_ise],[fullpar(ind).gi])
axis square 
box off
xlabel('W_{ISE}');
ylabel('g_{I}');
set(gca,'xtick',0:0.2:1);

