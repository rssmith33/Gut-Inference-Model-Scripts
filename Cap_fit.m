%% Function for fitting the gut inference model to subject data

% Samuel Taylor and Ryan Smith
% 1/5/2021

% subject:      filename of subject data (without .csv)
% data_dir:     root directory containing the subject data
% results_dir:  directory to save output of fitted results
function Cap_fit(subject, data_dir, results_dir)
    tic;
    
    % Fitting options -- currently set to winning model
    
    % a_mode: Learning rates for A matrix (IP)
    %   0: no learning
    %   1: one fitted learning rate
    %   2: split learning rate (vibrations and no vibrations)
    
    % b_mode: Learning rates for B matrix (pV)
    %   0: no learning
    %   1: one fitted learning rate
    %   2: split learning rate (vibrations and no vibrations)
    
    % IPdiff_on: Toggle for including IPdiff
    %   0: Do not include IPDiff
    %   1: Include IPdiff
    %
    % Note that IPdiff_on must be 0 if a_mode > 0.
    cap_opts = struct('a_mode', [0], 'b_mode', [2], 'IPdiff_on', [1]);
    cap_opts = struct2table(cap_opts);
        
    %% Add Subj Data
    file = [data_dir '/' subject '.csv'];
    disp(file)

    rawdat = readtable(file); %subject data

    % only grab trials of sufficient length (greater than 1.5s in len)
    subdat = rawdat(rawdat.Length > 1500, :); 

    last_time = 0;

    % If fitting multiple models, iterate across all model combinations
    % (currently only setup to fit one model at a time)
    for opt_i = 1:size(cap_opts, 1)
        options = cap_opts(opt_i, :);
                
        % Check for conditions of redundant models (where parameter setting
        % precludes the possibility of another parameter being set
        % a particular value
        if(options.a_mode > 0 && options.IPdiff_on > 0), continue; end

        disp(options)    
        opt_names = options.Properties.VariableNames;
        
        % Construct model name
        model = '';
        for opt = 1:numel(opt_names)
            model = strcat(model, opt_names{opt}, num2str(options.(opt_names{opt})), '_');
        end
                
        disp([results_dir '/' model '_' subject '.out.mat']);
        
        % Determine where the order of blocks
        block_order = unique(subdat.Block, 'stable');
        output = {};
        if block_order(2) == 4
            sub_type = 'Normal_First';
        else
            sub_type = 'Enhanced_First';
        end
        
        % Treating entire task as one block
        for block=1:1
            % Fit model to task data
            output = [output;       
                  Cap_block_fit(file, options, subdat) block_order(2*block) sub_type ...
            ];
        end

        save([results_dir '/' model '_' subject '.out.mat'], 'output');  
        
        disp(toc - last_time);
        last_time = toc;
    end    
    
    toc;
    
%     disp(g);
end