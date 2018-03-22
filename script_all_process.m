close all; clear; clc;

%% Load .mat
SUB_dir = uigetdir();
tmp_list = dir( fullfile(SUB_dir) );

SUB_list = tmp_list(arrayfun(@(x) x.name(1), tmp_list) ~= '.'); % get rid of '.' and '..' folders
clear tmp_list

summary_variables = { 'SubID', 'Section', 'order',...
                      'TargetLocation', 'PertProb', 'GoCueDelay',...
                      'TrialFailed', 'Perturbed', 'RT', 't2pVel',...
                      'pVel', 'manMovementOnset', 'manRT', 'manMT',...
                      'man_t2pVel', 'iniAngDevi', 'maxPathDevi', 'pathLength'};
                  
summary_table = [];
for i = 1:length(SUB_list) % subject
    tmp_path = fullfile(SUB_dir, SUB_list(i).name);
    tmp_list = dir( fullfile(tmp_path, '*.mat') );
    
    for j = 1:length(tmp_list) % section
        full_filename = fullfile(tmp_path, tmp_list(j).name);
        load(full_filename);

        [trial_no, load_table, target_table, tp_table] = convert_data2tables(data);
        % comupter behavior variables
        trial_behavior = behavior( trial_no );
        
        % compile a summary table
        trial_behavior.SubID = repmat(i, height(trial_behavior), 1);
        trial_behavior.Section = repmat(j, height(trial_behavior), 1);
        if isempty(summary_table)
            summary_table = trial_behavior(:, summary_variables);
        else
            summary_table = [summary_table; trial_behavior(:, summary_variables)];
        end
    end
end

save(fullfile(SUB_dir, 'summary_table.mat'), 'summary_table');