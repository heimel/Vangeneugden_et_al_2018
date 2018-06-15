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
              
        % surround feedback  excitation
        par.w_efs = 0;
        par.w_ifs = 0;
        
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
        
        % surround retinal input
        par.w_ers = 0;
        par.w_irs = 0;
        
        par.matchlva2v1 = true;

    case 1 % solution with non-local surround inhibition
        par.w_esi = 0.494;
        par.w_ie = 0.491;
        par.w_ei = 1.14;
        par.w_ese = 0.103;
        par.ge = 245;
        par.gi = 8.07;  
        par.matchlva2v1 = true;
    case 2 % solution with local surround inhibition
        % set
        par.w_ise = 0.75;
        % optimized
        par.w_ie = 0.053;
        par.w_ei = 0.93;
        par.w_ese = 0.65;
        % normalized
        par.ge = 101;
        par.gi = 13.4;
        par.matchlva2v1 = true;

    case 10
        par.w_esi = mean([0.5 0]);
        par.w_ie = mean([0.49 0.053]);
        par.w_ei = mean([1.13 0.93]);
        par.w_ese = mean([0.103 0.65]);
        par.ge = mean([245 101]);
        par.gi = mean([8.07 13.4]);  
        par.matchlva2v1 = true;
 
end


function x = steps(m) %#ok<DEFNU>
n = 3;
scale = 1.05;
x = [logspace(log10(1),log10(scale),n) logspace(log10(1/scale),log10(1),n)];
x(end) = [];
x = x * m;