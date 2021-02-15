%% Function for fitting the gut inference model to subject data (after pre-processing)

% Samuel Taylor and Ryan Smith 
% 1/5/2021

% file:         full path of csv file containing subject data
% options:      model fitting options in a table
% subdat:       contains preprocessed subject data (obtained from Cap_fit)
function results = Cap_block_fit(file, options, subdat)
    TpB = size(subdat,1); %trials per block

    % Store actions and observations in a structure array
    sub.u   = ones(1, 2, TpB);
    sub.o   = ones(1, 2, TpB);
    action  = ones(1, 1, TpB);
    vib     = ones(1, 1, TpB);
    
    block   = zeros(1, 1, TpB);
    
    for i = 1:TpB  
        vib (1, 1, i) = subdat.Vib(i);
        sub.o (1,2,i) = subdat.Vib(i) + 2; 

        action(1, 1,i) = subdat.Press(i); 
        sub.u (1, :,i) = [1 1];
        
        block(1, 1, i) = subdat.Block(i);
    end     

    o_all = sub.o;
    u_all = sub.u;

    % Count correct and incorrect choices
    for i = 1:size(o_all,3)
        if vib(:,:,i) == (action(:,:,i))
            correct(i) = 1;
        else
            correct(i) = 0;
        end
    end

    % Store accuracy of subject
    accuracy = mean(correct)*100;

    % Store true/false positive/negatives
    TN = zeros(1,size(o_all,3));
    TP = zeros(1,size(o_all,3));
    FN = zeros(1,size(o_all,3));
    FP = zeros(1,size(o_all,3));

    for i = 1:size(o_all,3)
        if vib(:,:,i) == action(:,:,i)
            if vib(:,:,i) == 0
                TN(i) = 1;
            elseif vib(:,:,i) == 1
                TP(i) = 1;
            end
        elseif vib(:,:,i) == 0
            FP(i) = 1;
        elseif vib(:,:,i) == 1
            FN(i) = 1;
        end
    end

    %storing true/false negatives/positives
    TP_FP_FN_TN = array2table([TP; FP; FN; TN]', 'VariableNames', {'TP', 'FP', 'FN', 'TN'}); 

    %% Params
    %--------------------------------------------------------------------------
    IP      = 0.95;    % precision of tone (0-1)
    eta     = 0.5;    % learning rate (between 0-1; default = 1)
    etaV    = 0.5;    % learning rate (between 0-1; default = 1)
    etaNV   = 0.5;    % learning rate (between 0-1; default = 1)
    pV      = 0.5;    % prior bias (0-1; higher = prior favoring detecting vibration)
    IPdiff  = 0.25;

    %% Invert model and try to recover original parameters:
    %==========================================================================

    params = struct(     ...
        'IP', IP,        ...    # Represents IP in paper (Interoceptive Precision)
        'pV', pV,        ...    # Represents pV in paper (prior)
        'etaA', eta,     ...    # Single learning rate for A matrix
        'etaAV', etaV,   ...    # Learning rate for vibrations in A matrix
        'etaANV', etaNV, ...    # Learning rate for no vibrations in A matrix
        'etaB', eta,     ...    # Single learning rate for B matrix
        'etaBV', etaV,   ...    # Learning rate for vibrations in B matrix
        'etaBNV', etaNV, ...    # Learning rate for no vibrations in B matrix
        'IPdiff', IPdiff ...    # Difference in IP between normal and enhanced blocks
    );

    % Generate model structure from model options
    MDP = Cap_gen_mdp(params, options);

    MDP.TpB = TpB; % trials per block

    MDP.action = action + 2;
    MDP.vib = vib;
    
    MDP.block = block;


    DCM.MDP   = MDP;                  % MDP model

    % Specify set of params to fit, based on selected
    % options for the model.
    DCM.field = {'IP' 'pV'};

    if options.b_mode == 2
        DCM.field = [DCM.field {'etaBV' 'etaBNV'}];
    elseif options.b_mode == 1
        DCM.field = [DCM.field {'etaB'}];
    end

    if options.a_mode == 2
        DCM.field = [DCM.field {'etaAV' 'etaANV'}];
    elseif options.a_mode == 1
        DCM.field = [DCM.field {'etaA'}];
    end
    
    
    if options.IPdiff_on && options.a_mode == 0
        DCM.field = [DCM.field {'IPdiff'}];
    end

    DCM.U     = {o_all};              % trial specification (stimuli)
    DCM.Y     = {u_all};              % responses (action)

    % Perform model inversion
    DCM       = Cap_inversion(DCM);

    %--------------------------------------------------------------------------
    % re-transform values and compare prior with posterior estimates
    %--------------------------------------------------------------------------


    field = fieldnames(DCM.M.pE);

    prior = zeros(1, size(field, 1));
    posterior = zeros(1, size(field, 1));

    % Extract the fitted parameters, returning the parameters
    % to the appropriate space
    for i = 1:length(field)
        disp(field{i});
        if strcmp(field{i},'etaB')
            prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
            posterior(i) = 1/(1+exp(-DCM.Ep.(field{i}))); 
        elseif strcmp(field{i},'etaA')
            prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
            posterior(i) = 1/(1+exp(-DCM.Ep.(field{i}))); 
        elseif strcmp(field{i},'etaBV')
            prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
            posterior(i) = 1/(1+exp(-DCM.Ep.(field{i})));
        elseif strcmp(field{i},'etaBNV')
            prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
            posterior(i) = 1/(1+exp(-DCM.Ep.(field{i})));
        elseif strcmp(field{i},'etaAV')
            prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
            posterior(i) = 1/(1+exp(-DCM.Ep.(field{i})));
        elseif strcmp(field{i},'etaANV')
            prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
            posterior(i) = 1/(1+exp(-DCM.Ep.(field{i})));
        elseif strcmp(field{i},'IP')
            prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
            posterior(i) = 1/(1+exp(-DCM.Ep.(field{i})));
        elseif strcmp(field{i},'pV')
            prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
            posterior(i) = 1/(1+exp(-DCM.Ep.(field{i})));
        elseif strcmp(field{i},'IPdiff')
            prior(i) = 1/(1+exp(-DCM.M.pE.(field{i})));
            posterior(i) = 1/(1+exp(-DCM.Ep.(field{i})));
        else
            prior(i) = exp(DCM.M.pE.(field{i}));
            posterior(i) = exp(DCM.Ep.(field{i}));
        end
    end

    % Save all relevant model information and return the results
    prior       = array2table(prior,        'VariableNames', field);
    posterior   = array2table(posterior,    'VariableNames', field);

    [model_acc, P_avg, button_accuracy, BP_avg, nobutton_accuracy, NBP_avg] = Cap_acc(posterior, options, file);
    avg_delay = mean(subdat.Delay, 'omitnan');
    
    results = {{file} prior posterior DCM accuracy TP_FP_FN_TN avg_delay model_acc P_avg button_accuracy BP_avg nobutton_accuracy NBP_avg};
end