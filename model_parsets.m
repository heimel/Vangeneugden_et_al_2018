function par = model_parsets(s)
%MODEL_PARSETS parameters set for vangeneugdenetal_model
% 
%  PAR = MODEL_PARSETS(S)
%
% 2018, Alexander Heimel

if nargin<1 || isempty(s)
    s = 1;
end

switch s
    case 0 % zero model with minimal parameters
        % e denotes V1 excitatory population
        % i denotes V1 inhibitory population
        % f denotes LVA excitatory population
        % j denotes LVA inhibitory population
        
        
        par.fi = @nonlinear; % transfer function

        par.rt_delay = 0.055; % s, delay from retina to V1 (from linear cutoff response Fig 3e)
        
        % decay times
        par.tau_e = 0.060; % s, from Ozeki et al.
        par.tau_i = 0.005; % s, inhibitory population faster
        par.tau_f = 0.060; % s
        par.tau_j = 0.005; % s, inhibitory population faster
        
        % feedforward excitation
        par.w_er = 1; % retinal input to V1 excitatory population
        par.w_ir = 1; % retinal input to V1 inhibitory population

        par.w_fe = 1; % V1 input to LVA excitatory population
        par.w_je = 1; % V1 input to LVA inhibitory population

        % local excitation
        par.w_ee = 0;
        par.w_ie = 0;
        par.w_ff = 0;
        par.w_jf = 0;
        
        % surround excitation
        par.w_ese = 0;
        par.w_ise = 0;
        par.w_fsf = 0;
        par.w_jsf = 0;
        
        % surround inhibition (inhibition from non-local interneurons)
        par.w_esi = 0;
        par.w_isi = 0;
        
        % local inhibition
        par.w_ei = 0;
        par.w_ii = 0;
        par.w_fj = 0;
        par.w_jj = 0;
        
        % feedback  excitation
        par.w_ef = 0.5; % arbitrary set to half of feedforward input
        par.w_if = par.w_ef; 
        
        % firing threshold
        par.t_e = 4;
        par.t_i = 4;
        par.t_f = 4;
        par.t_j = 4;
        
        % awake membrane potential shift
        par.t_ea = 1;
        par.t_ia = 1;
        par.t_fa = 1;
        par.t_ja = 1;
        
        % transfer function multiplier
        par.ge = 1;
        par.gi = 1;
        par.gf = 1;
        par.gj = 1;
        
        % nonlinearity awake multiplier
        par.ge_a = 1;
        par.gi_a = 1;
        par.gf_a = 1;
        par.gj_a = 1;
                
        par.matchlva2v1 = true;

    case 1 % solution with non-local surround inhibition
        % par.w_ise = 0; % no surround input to local inhibitory population
        par.w_esi = 0.62;
        par.w_ese = 0.18;
        par.w_ie = 0.491;
        par.w_ei = 1.14;
        
        % to normalize response
        par.ge = 245;
        par.gi = 8.07;  
    case 2    % solution with local surround inhibition
        par.w_ise = 1.0;
        par.w_esi = 0; % no inhibition from surrounding inhibitory population
        par.w_ese = 1.0;
        par.w_ie = 0;
        par.w_ei = 0.91;

        % to normalize response
        par.ge = 94;
        par.gi = 14;

    case 3 % solution with non-local surround inhibition and sigmoid
        par.fi = @sigmoid;
        par.w_er = 4.9 ;
        par.t_e = 10.0;
        par.t_i = 12.5;
        par.w_ie = 1.8;
        par.w_ei = 2.8;
        par.w_esi = 51;
        par.w_ese = 2.9;

        % to normalize response
        par.ge = 210;
        par.gi = 26;
    case 10 % used for track_solutions (mean of case 1 and 2)
        par.w_esi = mean([0.62 0]);
        par.w_ie = mean([0.491 0]);
        par.w_ei = mean([1.14 0.91]);
        par.w_ese = mean([0.18 1.0]);
        par.ge = mean([245 94]);
        par.gi = mean([8.07 4]);  
        
        
    case 123 % experimental model to play with parameters
        
          % par.w_ise = 0; % no surround input to local inhibitory population
        par.w_esi = 0.62;
        par.w_ese = 0.18;
        par.w_ie = 0.491;
        par.w_ei = 1.14;
        
        % to normalize response
        par.ge = 245;
        par.gi = 8.07;  
        
        par.w_ii = 0.6;
end


function x = steps(m) %#ok<DEFNU>
n = 3;
scale = 1.05;
x = [logspace(log10(1),log10(scale),n) logspace(log10(1/scale),log10(1),n)];
x(end) = [];
x = x * m;