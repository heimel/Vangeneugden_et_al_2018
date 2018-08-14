function [par,fullpar,err] = feedback_model( par,fullpar,optim,verbose,simmode)
% FEEDBACK_MODEL made to explain Vangeneugden et al data
%
%   [PAR,FULLPAR,ERR] = FEEDBACK_MODEL(PAR=1,FULLPAR=0,OPTIM=FALSE,VERBOSE=TRUE,SIMMODE='figures')
%    model just a single column, with surround input and feedback
%
%       PAR is parameter struct or parameter set number from model_parsets
%             when OPTIM is true, only parameters in PAR will be optimized
%             PAR can also be a struct with vectors for multiple parameter
%             values. PAR will then be expanded to a struct array
%             containing all combinations of parameters
%       FULLPAR is the full parameter struct or parameter struct set number
%             from which to take the parameters that are not set in PAR
%       OPTIM when true, a FMINSEARCH search will be done to find minimum
%       VERBOSE when true, makes figures and gives extra output
%       SIMMODE when 'figures' has simulation parameters for producing
%             figures for presentation, when 'optim' simulates a shorter
%             time for faster parameter optimization
%
% 2018, Alexander Heimel

if exist ('OCTAVE_VERSION', 'builtin') 
    more off
else
    rng('default');
end

if nargin<1 || isempty(par)
    % model parameters
    par = 1;
end
if ~isstruct(par)
    par = model_parsets(par);
end
if nargin<2 || isempty(fullpar)
    fullpar = 0; % minimal empty model
end
if ~isstruct(fullpar)
    fullpar = model_parsets(fullpar);
end
fullpar =  mergestruct(model_parsets(0),fullpar);

if nargin<3 || isempty(optim)
    optim = false;
end
if nargin<4 || isempty(verbose)
    verbose = true;
end
if nargin<5 || isempty(simmode)
    simmode = 'optim';
end

warning('Off','MATLAB:legend:IgnoringExtraEntries'); 

% Simulation parameters
simpar.dt = 0.001; % s, simulation time step
simpar.r_onsetwidth = 0.1; %s
switch simmode
    case 'figures'
        simpar.max_t = 1; % s, simulation time
        simpar.stim_onset = 0.5; %s
        simpar.stim_offset = 1; %s
        simpar.noise = 0; 
    case 'optim'
        simpar.max_t = 0.5; % s, simulation time
        simpar.stim_onset = 0; %s
        simpar.stim_offset = 0.5; %s
        simpar.noise = 0; % 0.2 for more robust optimization 
end
simpar.nsamples = ceil(simpar.max_t/simpar.dt);

onset_i = round(simpar.stim_onset/simpar.dt)+1;
offset_i =  round(simpar.stim_offset/simpar.dt)+1;

simpar.R = zeros(simpar.nsamples,1);
contrast = 1;
simpar.R(onset_i:min(simpar.nsamples,offset_i)) = contrast;
simpar.R(onset_i:onset_i+round(simpar.r_onsetwidth/simpar.dt)) = linspace(0,contrast,round(simpar.r_onsetwidth/simpar.dt)+1);

disp(simpar);

% Measured si(aw,fb)
mea.si(2,2) = 0.44; % awake, with feedback
mea.si(2,1) = 0.26; % awake, without feedback
mea.si(1,2) = 0.27; % anesthetized, with feedback, estimate from Fig. 3f
mea.si(1,1) = 0.31; % anesthetized, without feedback, estimate from Fig. 3f

% Only test awake pars
mea.si(1,:) = mea.si(2,:); 

if isfield(par,'fi')
    simpar.fi = par.fi;
    par = rmfield(par,'fi');
elseif isfield(fullpar,'fi')
    simpar.fi = fullpar.fi;
else
    simpar.fi = @sigmoid;
end

par = find_bestpar(par,fullpar,simpar,mea,false);

if optim
    [err,sim] = compute_statstate(par,fullpar,simpar,mea,false);
    disp('Starting model');
    disp(['Measured SI = ' mat2str(mea.si,2)]);
    disp(['Model SI    = ' mat2str(sim.si,2)]);
    disp(['Error = ' num2str(err,4) ]);
    
    flds = fieldnames(par);
    fun = @(p) compute_statstate(cell2struct(num2cell(p)',flds),fullpar,simpar,mea,false);
    parv = struct2cell(par);
    parv = [parv{:}];
    options = optimset;
    options.MaxFunEvals = 10000;
    best_par_vec = fminsearch(fun,parv,options);
    par = cell2struct(num2cell(best_par_vec)',flds);
end

[err,sim] = compute_statstate(par,fullpar,simpar,mea,verbose);
disp('Final model');
disppar( par )
disp(['Measured SI = ' mat2str(mea.si,2)]);
disp(['Model SI    = ' mat2str(sim.si,2)]);
disp(['Error = ' num2str(err,4) ]);

if verbose
    compute_contrast_curve(par,fullpar,simpar)
end

if verbose % extra paper figures
    figure
    
    % bar plot measured vs model
    subplot(2,2,1);
    hold on
    bar([0.5 1.8 ],[mea.si(1,2)  sim.si(1,2) ],0.35,'k');
    bar([1 2.3 ],[mea.si(1,1)  sim.si(1,1) ],0.35,'b');
    ylabel('SI')
    ylim([0 0.6]);
    set(gca,'ytick',0:0.1:0.6);
    set(gca,'xtick',[0.75 2.05]);
    set(gca,'xticklabel',{'Exp.','Model'});
    xlim([0. 5]);
    
    % nonlinear fig
    if isfield(par,'fi_par')
        fi = @(x) simpar.fi(x,par.fi_par * par.fi_para);
    else
        fi = simpar.fi;
    end
    subplot(2,2,2);
    x = (-2:0.01:1);
    y = fi(x-4);
    y = y/max(y);
    plot(x,y,'k');
    box off
    axis square
    xlim([-0.5 1]);
    xlabel('Input');
    ylabel('Output');
end

fullpar = mergestruct(fullpar,par);
if isfield(fullpar,'matchlva2v1') && fullpar.matchlva2v1
fullpar = make_matchlva2v1(fullpar);
end


function best_par = find_bestpar(par,fullpar,simpar,mea,verbose)

% expand parameters
flds = fieldnames(par);
for f=1:length(flds)
    field = flds{f};
    i = 1;
    while i<=length(par)
        p = par(i);
        v = par(i).(field);
        if length(v)>1
            par(i).(field) = v(1);
            for j = 2:length(v)
                par(end+1) = p; %#ok<AGROW>
                par(end).(field) = v(j);
            end
        else
            i = i+1;
        end
    end
end

best_par = par(1);
min_error = compute_statstate(best_par,fullpar,simpar,mea,false);
par = par(randperm(length(par)));
disp(['Error = ' num2str(min_error,4)]);
disp(['Testing ' num2str(length(par)) ' conditions']);

for i=1:length(par)
    if mod(i,10)==0
        disp(['Par # = ' num2str(i) ' of ' num2str(length(par)) ...
            ', besterror=' num2str(min_error,4)]);
    end
    
    err = compute_statstate(par(i),fullpar,simpar,mea,verbose);
    if err<min_error
        min_error = err;
        best_par = par(i);
        disppar(best_par);
        disp(['Error = ' num2str(err,4)]);
        disp(['Par # = ' num2str(i) ' of ' num2str(length(par))]);
        compute_statstate(best_par,fullpar,simpar,mea,false);
    end
end



function compute_contrast_curve(par,fullpar,simpar)
%COMPUTE_CONTRAST_CURVE varies retinal contrast from 0 to 1 to calculate a contrast tuning curve 

par = mergestruct(fullpar,par);

contrast = linspace(0,2,40);
E = zeros(size(contrast));
I = zeros(size(contrast));
F = zeros(size(contrast));

for aw = 1:2
    figure('Name',['Awake=' num2str(aw)]);
    for fb = 1:2
        par.fb = fb-1;
        
        hold on
        for sr = 1:2
            par.sr = sr-1;
            
            subplot(2,2,(fb-1)*2+sr);
            title(['fb=' num2str(par.fb)  ',' 'sr=' num2str(par.sr) ]);
            ylim([0 1.2]);
            hold on
            for i = 1:length(contrast)
                par.aw = aw-1;
                
                csimpar = simpar;
                csimpar.R = contrast(i)*csimpar.R;
                obs = compute_dynamics( par,csimpar,[]);
                E(i) = obs.re;
                I(i) = obs.ri;
                F(i) = obs.rf;
                
            end
            plot(contrast*100,E,'k-');
            plot(contrast*100,I,'r-');
            plot(contrast*100,F,'g-');
            xlim([0 120]);
            ylabel('Response (a.u.)');
            xlabel('Contrast (%)');
            ind = find(contrast>1,1);
            plot([100 100],ylim,'k--');
            plot(xlim,[E(ind) E(ind)],'k--');
            set(gca,'xtick',0:20:100);
            set(gca,'xticklabel',{'','20','','60','','100'});
            set(gca,'ytick',0:0.2:1.2);
            set(gca,'yticklabel',{'0','','0.4','','0.8','','1.2'});
            axis square
        end
    end
end

xlabel('Contrast');
ylabel('Response');


function [err,sim,E] = compute_statstate(par,fullpar,simpar,mea,verbose)
%COMPUTE_STATSTATE calls compute_dynamics for all 8 conditions
%  awake,anesthetized;without feedback,with feedback,center
%  only,center+surround

par = mergestruct(fullpar,par);

E = zeros(2,2,2);
F = zeros(2,2,2);
I = zeros(2,2,2);

for aw = 1:2
    if verbose
        fig = figure('Name',['Awake=' num2str(aw)]);
    else
        fig = 0;
    end
    for fb = 1:2
        for sr = 1:2
            par.aw = aw-1;
            par.fb = fb-1;
            par.sr = sr-1;
            obs(aw,fb,sr) = compute_dynamics( par,simpar, fig ); %#ok<AGROW>
            E(aw,fb,sr) = obs(aw,fb,sr).re;
            I(aw,fb,sr) = obs(aw,fb,sr).ri;
            F(aw,fb,sr) = obs(aw,fb,sr).rf;
        end % sr
    end % fb
end % aw

for aw = 1:2
    for fb = 1:2
        sim.si(aw,fb) = (E(aw,fb,1)-E(aw,fb,2))/(E(aw,fb,1)+0.001);
    end
end

% E(aw,fb,sr)
err = 0;

% set awake v1 center response with feedback to 1
err = err + (E(2,2,1) - 1)^2;

% set awake v1 center I response with feedback to 1
err = err + (I(2,2,1) - 1)^2;

% do not want oscillation
err = err + obs(2,2,1).se^2;
err = err + obs(2,2,2).se^2;

% anesthesia no effect of feedback center
err = err + (E(1,1,1)-E(1,2,1))^2;

% compute distance from measured si
for aw = 1:2
    for fb = 1:2
        err = err + (mea.si(aw,fb)-sim.si(aw,fb))^2;
    end
end

% awake no effect of feedback center
err = err + (E(2,1,1)-E(2,2,1))^2;

% awake and anesthetized surround responses should be similar
err = err + 0.2*(E(1,2,2)-E(2,2,2))^2 ;

% awake and anesthetized LVA surround responses should be similar
err = err + 0.2*(F(1,2,2)-F(2,2,2))^2 ;

% punish response to blank stimulus
par.aw = 1;
par.fb = 1;
par.sr = 1;
simpar.R = 0 * simpar.R;
zerocontrast = compute_dynamics( par,simpar);
err = err + zerocontrast.re^2;

% also for anesthesia
par.aw = 0;
par.fb = 1;
par.sr = 1;
simpar.R = 0 * simpar.R;
zerocontrast = compute_dynamics( par,simpar );
err = err + zerocontrast.re^2;

% punish negative exponents and negative weights
if par.w_er<0
    err = err *1.5;
end
if par.w_ir<0
    err = err *1.5;
end
if par.w_ee<0
    err = err *1.5;
end
if par.w_ie<0
    err = err *1.5;
end
if par.w_fe<0
    err = err *1.5;
end
if par.w_ei<0
    err = err *1.5;
end
if par.w_ii<0
    err = err *1.5;
end
if par.w_ef<0
    err = err *1.5;
end
if par.w_if<0
    err = err *1.5;
end

function obs = compute_dynamics(  par,simpar, h )

R = simpar.R;

if nargin<3 || isempty(h)
    h = 0;
end

nsamples = length(R);

% average firing rates
E = zeros(size(R)); % V1 excitatory neurons
I = zeros(size(R)); % V1 inhibitory neurons
F = zeros(size(R)); % LVA excitatory neurons
J = zeros(size(R)); % LVA inhibitory neurons

if isfield(par,'fi_par')
    fi = @(x) simpar.fi(x,par.fi_par * par.fi_para^(1-par.aw));
else
    fi = simpar.fi;
end

if isfield(par,'matchlva2v1') && par.matchlva2v1
    par = make_matchlva2v1(par);
end

if 1 % link all feedforward and feedback strengths
    par.w_ir = par.w_er;
    par.w_fe = par.w_er;
    par.w_je = par.w_er;
    par.w_ef = 0.5 * par.w_er;
    par.w_if = 0.5 * par.w_er;
end

delta_e = simpar.dt/par.tau_e;
delta_i = simpar.dt/par.tau_i;
delta_f = simpar.dt/par.tau_f;
delta_j = simpar.dt/par.tau_j;

dt = simpar.dt;

ge = par.ge * par.ge_a^(1-par.aw);
gi = par.gi * par.gi_a^(1-par.aw);
gf = par.gf * par.gf_a^(1-par.aw);
gj = par.gj * par.gj_a^(1-par.aw);

t_e = par.t_e * par.t_ea^(1-par.aw);
t_i = par.t_i * par.t_ia^(1-par.aw);
t_f = par.t_f * par.t_fa^(1-par.aw);
t_j = par.t_j * par.t_ja^(1-par.aw);

rt_delay_sample = round(par.rt_delay/dt);

for t = 3:nsamples
    % V1 excitatory neurons
    E(t) = (1-delta_e)*E(t-1) + ...
        delta_e * ge * fi(   (...
        - t_e ...
        + par.w_er * R(max(1,t-1-rt_delay_sample)) ...
        + par.w_ee * E(t-1) ...
        - par.w_ei * I(t-1) ...
        + par.w_ese * par.sr * E(t-1) ...
        - par.w_esi * par.sr * I(t-1) ...
        + par.w_ef * par.fb * F(t-1) ...
        + simpar.noise*(rand(1)-0.5) ...
        ));
    
    % V1 inhibitory neurons
    I(t) = (1-delta_i)*I(t-1) + ...
        delta_i * gi *  fi(( ...
        - t_i ...
        + par.w_ir * R(max(1,t-1-rt_delay_sample)) ...
        + par.w_ie * E(t-1) ...
        - par.w_ii * I(t-1) ...
        + par.w_ise * par.sr * E(t-1) ...
        - par.w_isi * par.sr * I(t-1) ...
        + par.w_if * par.fb * F(t-1) ...
        ));
    
    % LVA excitatory neurons
    F(t) = (1-delta_f)*F(t-1) + ...
        delta_f * gf * fi(  ...
        - t_f ...
        + par.w_fe * E(t-1) ...
        + par.w_ff * F(t-1) ...
        - par.w_fj * J(t-1) ...
        );
        
    % LVA inhibitory neurons
    J(t) = (1-delta_j)*J(t-1) + ...
        delta_j * gj * fi((  ...
        - t_j ...
        + par.w_je * E(t-1) ...
        + par.w_jf * F(t-1) ...
        - par.w_jj * J(t-1) ...
        ));
    
end % timestep t

if h~=0
    figure(h);
    
    t = (1:length(R))*dt -simpar.stim_onset;
    c = 'ck';
    
    subplot(2,2,par.sr+1);
    box off
    hold on
    plot(t,E,['-' c(par.fb+1)]);
    axis square
    ylim([0 1.4]);
    if par.sr==0
        ylabel('Pyr Response (a.u.)');
    else
        legend('No feedback','With feedback');
        legend boxoff
    end
    title(['fb=' num2str(par.fb)  ',' 'sr=' num2str(par.sr) ]);
    xlim( [-0.5 0.5]);
    set(gca,'xtick',-0.4:0.2:0.4);
    set(gca,'xticklabel',{});
    set(gca,'ytick',0:0.2:1.4);
    set(gca,'yticklabel',{'0','','0.4','','0.8','','1.2',''});
    
    c = 'mr';
    
    subplot(2,2,par.sr+1+2);
    hold on
    box off
    plot(t,I,['-' c(par.fb+1)]);
    axis square
    xlabel('Time (s)');
    if par.sr==0
        ylabel('Inh Response (a.u.)');
    else
        legend('No feedback','With feedback');
        legend boxoff
    end
    title(['fb=' num2str(par.fb)  ',' 'sr=' num2str(par.sr) ]);
    xlim( [-0.5 0.5]);
    ylim([0 1.4]);
    set(gca,'xtick',-0.4:0.2:0.4);
    set(gca,'xticklabel',{'-0.4','','0','','0.4'});
    set(gca,'ytick',0:0.2:1.4);
    set(gca,'yticklabel',{'0','','0.4','','0.8','','1.2',''});
    
end

% take second half of simulated time
obs.re = mean(E(round(end/2):end));
obs.ri = mean(I(round(end/2):end));
obs.rf = mean(F(round(end/2):end));
obs.se = std(E(round(end/2):end));

function par = make_matchlva2v1(par)
    par.w_fe = par.w_er;
    par.w_je = par.w_ir;
    par.w_jf = par.w_ie;
    par.w_fj = par.w_ei;
    par.t_f = par.t_e;
    par.t_j = par.t_i;
    par.gf = par.ge;
    par.gj = par.gi;
    par.tau_j = par.tau_i;
    par.w_jj = par.w_ii;
    par.w_ff = par.w_ee;

