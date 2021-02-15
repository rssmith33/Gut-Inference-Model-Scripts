%% Model inversion script used to estimate model parameters

function [DCM] = Cap_inversion(DCM)

% MDP inversion using Variational Bayes
% FORMAT [DCM] = spm_dcm_mdp(DCM)

% If simulating - comment out section on line 196
% If not simulating - specify subject data file in this section 

%
% Expects:
%--------------------------------------------------------------------------
% DCM.MDP   % MDP structure specifying a generative model
% DCM.field % parameter (field) names to optimise
% DCM.U     % cell array of outcomes (stimuli)
% DCM.Y     % cell array of responses (action)
%
% Returns:
%--------------------------------------------------------------------------
% DCM.M     % generative model (DCM)
% DCM.Ep    % Conditional means (structure)
% DCM.Cp    % Conditional covariances
% DCM.F     % (negative) Free-energy bound on log evidence
% 
% This routine inverts (cell arrays of) trials specified in terms of the
% stimuli or outcomes and subsequent choices or responses. It first
% computes the prior expectations (and covariances) of the free parameters
% specified by DCM.field. These parameters are log scaling parameters that
% are applied to the fields of DCM.MDP. 
%
% If there is no learning implicit in multi-trial games, only unique trials
% (as specified by the stimuli), are used to generate (subjective)
% posteriors over choice or action. Otherwise, all trials are used in the
% order specified. The ensuing posterior probabilities over choices are
% used with the specified choices or actions to evaluate their log
% probability. This is used to optimise the MDP (hyper) parameters in
% DCM.field using variational Laplace (with numerical evaluation of the
% curvature).
%
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_mdp.m 7120 2017-06-20 11:30:30Z spm $

% OPTIONS
%--------------------------------------------------------------------------
    ALL = false;

% specify prior expectations and covariance
%--------------------------------------------------------------------------
    prior_variance = 2^-2;

    for i = 1:length(DCM.field)
        field = DCM.field{i};
        try
            param = DCM.MDP.(field);
            param = double(~~param);
        catch
            param = 1;
        end
        if ALL
            pE.(field) = zeros(size(param));
            pC{i,i}    = diag(param);
        else
            if strcmp(field,'IP')
                pE.(field) = log(0.95/(1-0.95));   % in logit-space - bounded between 0 and 1
                pC{i,i}    = prior_variance;
            elseif strcmp(field,'pV')
                pE.(field) = log(0.5/(1-0.5));     % in logit-space - bounded between 0 and 1
                pC{i,i}    = prior_variance;
            elseif strcmp(field,'etaA')
                pE.(field) = log(0.5/(1-0.5));     % in logit-space - bounded between 0 and 1
                pC{i,i}    = prior_variance;
            elseif strcmp(field,'etaAV')
                pE.(field) = log(0.5/(1-0.5));     % in logit-space - bounded between 0 and 1
                pC{i,i}    = prior_variance;
            elseif strcmp(field,'etaANV')
                pE.(field) = log(0.5/(1-0.5));     % in logit-space - bounded between 0 and 1
                pC{i,i}    = prior_variance;
            elseif strcmp(field,'etaB')
                pE.(field) = log(0.5/(1-0.5));     % in logit-space - bounded between 0 and 1
                pC{i,i}    = prior_variance;
            elseif strcmp(field,'etaBV')
                pE.(field) = log(0.5/(1-0.5));     % in logit-space - bounded between 0 and 1
                pC{i,i}    = prior_variance;
            elseif strcmp(field,'etaBNV')
                pE.(field) = log(0.5/(1-0.5));     % in logit-space - bounded between 0 and 1
                pC{i,i}    = prior_variance;
            elseif strcmp(field,'IPdiff')
                pE.(field) = log(0.2/(1-0.2));     % in logit-space - bounded between 0 and 1
                pC{i,i}    = prior_variance;
            else
                pE.(field) = 0;      
                pC{i,i}    = prior_variance;
            end
        end
    end

    pC      = spm_cat(pC);

% model specification
%--------------------------------------------------------------------------
    M.L     = @(P,M,U,Y)spm_mdp_L(P,M,U,Y);  % log-likelihood function
    M.pE    = pE;                            % prior means (parameters)
    M.pC    = pC;                            % prior variance (parameters)
    M.mdp   = DCM.MDP;                       % MDP structure

% SMT -------------+
    M.noprint = true;    
    M.nograph = true;
% end SMT ---------+

% Variational Laplace
%--------------------------------------------------------------------------

    [Ep,Cp,F] = spm_nlsi_Newton(M,DCM.U,DCM.Y);

% Store posterior densities and log evidnce (free energy)
%--------------------------------------------------------------------------
    DCM.M   = M;
    DCM.Ep  = Ep;
    DCM.Cp  = Cp;
    DCM.F   = F;


    return

function L = spm_mdp_L(P,M,U,Y)
% log-likelihood function
% FORMAT L = spm_mdp_L(P,M,U,Y)
% P    - parameter structure
% M    - generative model
% U    - inputs
% Y    - observed repsonses
%__________________________________________________________________________

    if ~isstruct(P); P = spm_unvec(P,M.pE); end

% multiply parameters in MDP
%--------------------------------------------------------------------------
    mdp   = M.mdp;
    field = fieldnames(M.pE);
    for i = 1:length(field)
        if strcmp(field{i},'IP')
            mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
        elseif strcmp(field{i},'pV')
            mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
        elseif strcmp(field{i},'etaA')
            mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
        elseif strcmp(field{i},'etaAV')
            mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
        elseif strcmp(field{i},'etaANV')
            mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
        elseif strcmp(field{i},'etaB')
            mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
        elseif strcmp(field{i},'etaBV')
            mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
        elseif strcmp(field{i},'etaBNV')
            mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
        elseif strcmp(field{i},'IPdiff')
            mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
        else
            mdp.(field{i}) = exp(P.(field{i}));
        end
    end


% discern whether learning is enabled - and identify unique trials if not
%--------------------------------------------------------------------------
    if any(ismember(fieldnames(mdp),{'a','b','d','c','d','e'}))
        j = 1:numel(U);
        k = 1:numel(U);
    else
% find unique trials (up until the last outcome)
%----------------------------------------------------------------------
        u       = spm_cat(U');
        [i,j,k] = unique(u(:,1:(end - 1)),'rows');
    end

    L     = 0;

    mdp_temp = Cap_gen_mdp(mdp, mdp.fit_options);

    [MDP(1:mdp.TpB)]   = deal(mdp_temp);
    
    for idx_trial = 1:size(MDP,2)
        MDP(idx_trial).o = [U{1}(:,:,idx_trial)];
        MDP(idx_trial).u = [Y{1}(:,:,idx_trial)];
    end

    for i = 1:mdp.TpB
        MDP(i).action = mdp(1).action(:,:,i);
        
        % if normal trial and using pVdiff
        if mdp(1).block(:, :, i) == 4 && mdp.fit_options.IPdiff_on == 1
            MDP(i).A{1}(:,:) = [1   0                           0;                          % start
                                0   (mdp.IP - mdp.IPdiff)       ((1-mdp.IP) + mdp.IPdiff);  % no vib
                                0   ((1-mdp.IP) + mdp.IPdiff)   (mdp.IP - mdp.IPdiff)];     % vib
        
        end
    end
 
% solve MDP and accumulate log-likelihood
%--------------------------------------------------------------------------
    MDP  = spm_MDP_VB_X_Cap(MDP);

    for j = 1:mdp.TpB
            L = L + log(MDP(j).X{1,1}(MDP(j).action,2) + eps);
    end
    
    clear('MDP')

% end

    fprintf('LL: %f \n',L)
