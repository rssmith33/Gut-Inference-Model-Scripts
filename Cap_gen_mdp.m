%% Function for constructing the generative model for gut inference

% Samuel Taylor and Ryan Smith 
% 1/5/2021

% p:        stores the field of parameters and initial values for the
%           parameters
% options:  options for fitting the model
function mdp = Cap_gen_mdp(p, options)
%% Set up gut inference model
%__________________________________________________________________________

%Prior beliefs (D) about initial states
%--------------------------------------------------------------------------
D{1} = [1 0 0]';           % context:         {'start' 'no vibration' 'vibration'}

% mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end
No    = [3];         % vibration signal (start, no vibration, vibration)
Ng    = numel(No);
for g = 1:Ng
    A{g} = zeros([No(g),Ns]);
end


% A  matrix (generative process)
%--------------------------------------------------------------------------          
A{1}(:,:) =   [1 0    0     ;  % start
               0 p.IP 1-p.IP;  % no vib
               0 1-p.IP p.IP]; % vib
        
for g = 1:Ng
    A{g} = double(A{g});
end

a = A; % only assigns later if a-learning is enabled

% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

% if b-learning is enabled
if options.b_mode > 0
    B{1} = [0     0 0; % start
            1-0.5 1 0; % no vib
            0.5   0 1];% vib
else
    B{1} = [0      0 0; % start
            1-p.pV 1 0; % no vib
            p.pV   0 1];% vib  
end

   
b{1} = [0     0 0; % start
       1-p.pV 1 0; % no vib
       p.pV   0 1];% vib

T=2;

% priors: (utility) C
%--------------------------------------------------------------------------
C{1}     = zeros(No(1),T); % No explicit preferred observations in this model

%% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.T = T;                      % number of moves
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states
mdp.erp = 1;                    % parameter for neural predictions 
                                % (1 indicates that it is not used) 

label.factor{1}   = 'Vib';   label.name{1}    = {'start','no vib','vib'};
label.modality{1} = 'Vib';  label.outcome{1} = {'start','no vib','vib'};

mdp.label = label;

% if learning
%--------------------------------------------------------------------------
    
% only if learning is enabled for b
if options.b_mode > 0
    mdp.b = b;
end

% only if learning is enabled for a
if options.a_mode > 0
    mdp.a = a; 
end

%--------------------------------------------------------------------------

% mdp.etaA    = .5;  % learning rate for A matrix
% mdp.etaB    = .5;  % learning rate for B matrix

if options.b_mode == 2      % if split learning rate for B
    mdp.etaBV   = p.etaBV;      % learning rate for B matrix when vibration occurs
    mdp.etaBNV  = p.etaBNV;     % learning rate for B matrix when no vibration occurs
elseif options.b_mode == 1  % if fit single learning rate for B
    mdp.etaB   = p.etaB;
end
    
if options.a_mode == 2      % if split learning rate for A
    mdp.etaAV   = p.etaAV;      % learning rate for A matrix when vibration occurs
    mdp.etaANV  = p.etaANV;     % learning rate for A matrix when no vibration occurs
elseif options.a_mode == 1  % if fit single learning rate A
    mdp.etaA   = p.etaA;
end

% interoceptive precision difference parameter (described in paper)
if options.IPdiff_on == 1 && options.a_mode == 0
    mdp.IPdiff = p.IPdiff;
end

mdp       = spm_MDP_check(mdp);

mdp.fit_options = options;

end