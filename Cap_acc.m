%% Function for computing the accuracy of a model fit on the subject's data
 
% Samuel Taylor and Ryan Smith 
% 1/5/2021

% params:       Values of fitted parameter estimates
% options:      model fitting options in a table
% data:         filepath of subject data

function [model_acc, P_avg, button_accuracy, BP_avg, nobutton_accuracy, NBP_avg] = Cap_acc(params, options, data)
    %--------------------------------------------------------------
    mdp = Cap_gen_mdp(params, options);

    %% Add Subj Data

    rawdat = readtable(data); %subject data

    % only grab trials of sufficient length (greater than 1.5s in len)
    subdat = rawdat(rawdat.Length > 1500, :); 

    TpB = size(subdat,1); %trials per block

    sub.u   = ones(1, 2, TpB);
    sub.o   = ones(1, 2, TpB);
    action  = ones(1, 1, TpB);
    vib     = ones(1, 1, TpB);
    block   = ones(1, 1, TpB);

    for i = 1:TpB  
        vib (1, 1, i) = subdat.Vib(i);
        sub.o (1,2,i) = subdat.Vib(i) + 2; 

        action(1, 1,i) = subdat.Press(i);
        sub.u (1, 1:2,i) = [subdat.Press(i) subdat.Press(i)] + 2; 
        sub.u (1, :,i) = [1 1];
        
        block(1, 1, i) = subdat.Block(i);
    end     

    mdp.o = [];
    mdp.vib = [];
    mdp.action = [];
    mdp.u = [];
    mdp.block = [];
    
    for i = 1:TpB
        MDP(i) = mdp;  
        MDP(i).action = action(:,:,i)+2;
        MDP(i).vib = vib(i);
        MDP(i).o = sub.o(:,:,i);
        MDP(i).u = sub.u(:,:,i);
        MDP(i).block = block(:, :, i);
        
        if MDP(i).block == 4 && MDP(i).fit_options.IPdiff_on == 1
            MDP(i).A{1}(:,:) = [1   0                           0;                                      % start
                                0   (params.IP - params.IPdiff)       ((1-params.IP) + params.IPdiff);  % no vib
                                0   ((1-params.IP) + params.IPdiff)   (params.IP - params.IPdiff)];     % vib
        end
    end

    % solve MDP and accumulate log-likelihood
    %-------------------------------------------------------------
    MDP_temp  = spm_MDP_VB_X_Cap(MDP);
    
%     assignin('base', 'genMDP', MDP_temp);

    P = 0;

    MDP = MDP_temp;
    
    model_correct = 0;

    good_answer = 0;
    
    button_correct = 0;
    press = 0;
    
    nobutton_correct = 0;
    nopress = 0;
    
    for i = 1:length(MDP)
        P(i) = MDP(i).X{1,1}(MDP(i).action,2);

        [~, arg] =  max(MDP(i).X{1,1}(:,2));
        
        if arg == MDP(i).action
            model_correct = model_correct + 1;
        end
        
        if arg == MDP(i).o(2)
            good_answer = good_answer + 1;
        end
        
        if MDP(i).action == 3
            BP(press + 1) = MDP(i).X{1,1}(3,2);
            
            if arg == MDP(i).action
                button_correct = button_correct + 1;
            end
            
            press = press + 1;
        elseif MDP(i).action == 2
            NBP(nopress + 1) = MDP(i).X{1,1}(2,2);
            
            if arg == MDP(i).action
                nobutton_correct = nobutton_correct + 1;
            end
            
            nopress = nopress + 1;
        end
    end

    P = P';
    P_avg = sum(P)/size(P,1);
    
    % Accuracy metrix for button press trials
    BP_avg = sum(BP) / press;
    button_accuracy = button_correct / press;
    
    % Accuracy metrics for non-button press trials
    NBP_avg = sum(NBP) / nopress;
    nobutton_accuracy = nobutton_correct / nopress;

    % Overall accuracy
    model_acc = model_correct / length(MDP);
end