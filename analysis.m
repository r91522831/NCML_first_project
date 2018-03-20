%%
close all; clear; clc;
load('/Users/yen-hsunwu/Dropbox (ASU)/mat files/summary_table.mat')

%%
trial_subid = summary_table{:, 'SubID'};
subid = [1, 2, 3];
trial_section = summary_table{:, 'Section'};
section = [1, 2, 3];
trial_target_location = summary_table{:, 'TargetLocation'};
target_location = {'L', 'M', 'R'};
trial_pert_lvl = cell2mat(summary_table{:, 'PertProb'});
pert_lvl = [75, 50, 25];
trial_godelay_lvl = cell2mat(summary_table{:, 'GoCueDelay'});
godelay_lvl = [80, 150, 220, 300];
trial_failed = cell2mat(summary_table{:, 'TrialFailed'});
trial_perturbed = summary_table{:, 'Perturbed'};
trial_no_response = ( summary_table{:, 'manMovementOnset'} == 0);

%% check the real perturbation probability
real_pertprob = zeros(3, 1);
for i = 1:3
    real_pertprob(i, 1) = sum( summary_table.Perturbed(trial_pert_lvl == pert_lvl(i), :) ) / length( summary_table.Perturbed(trial_pert_lvl == pert_lvl(i) , :) );
end

%%
% get subject #1, section 3, ...
%     target in Middle, ...
%     perturb probability 50, ...
%     go cue delay 150, ...
%     trial failed 0, ...
%     trial perturbed 1
%{
ind_target_location = cellfun(@eq, trial_target_location, repmat({'M'}, size(trial_target_location)));

demo_RT = summary_table.manRT( trial_subid == 1 & trial_section == 3 & ...
                               ind_target_location & ...
                               trial_pert_lvl == 50 & ...
                               trial_godelay_lvl == 150 & ...
                               trial_failed == 0 & ...
                               trial_perturbed == 1);
                           
ind_demo_RT = find( trial_subid == 1 & trial_section == 3 & ...
                    ind_target_location & ...
                    trial_pert_lvl == 50 & ...
                    trial_godelay_lvl == 150 & ...
                    trial_failed == 0 & ...
                    trial_perturbed == 1 );
%}
                
%%
% Compute mean and std of the variables in all trials, except trials with
% no response for individual subject.
% When a trial is failed due to not response, the manMovementOnset equals
% to zero!

tmp_table = [];
for i = subid
    for s = section
        if s ~= 1
            % process all subjects for the other section
            for pp = 1:length(pert_lvl)
                for gd = 1:length(godelay_lvl)
                    cond = ( trial_subid == subid(i) & trial_section == s & ...
                             trial_pert_lvl == pert_lvl(pp) & ...
                             trial_godelay_lvl == godelay_lvl(gd)& ...
                             ~trial_no_response );
                    cond_unpp = (cond & ~trial_perturbed & ~trial_failed);
                    
                    tmp_row = { i, s, pert_lvl(pp), godelay_lvl(gd), summary_table.manRT(cond), ...
                                                                     summary_table.pVel(cond), ...
                                                                     summary_table.man_t2pVel(cond), ...
                                                                     summary_table.manMT(cond), ...
                                                                     summary_table.pVel(cond_unpp), ...
                                                                     summary_table.man_t2pVel(cond_unpp), ...
                                                                     summary_table.manMT(cond_unpp), ...
                                                                     summary_table.iniAngDevi(cond_unpp) };

                    tmp_table = [tmp_table; tmp_row];
                end
            end
        else
            % process all subjects for the section 1
            cond = ( trial_subid == subid(i) & trial_section == 1 & ~trial_no_response);
            cond_unpp = cond;
            tmp_row = {i, s, 0, 0, summary_table.manRT(cond), ...
                                   summary_table.pVel(cond), ...
                                   summary_table.man_t2pVel(cond), ...
                                   summary_table.manMT(cond), ...
                                   summary_table.pVel(cond_unpp), ...
                                   summary_table.man_t2pVel(cond_unpp), ...
                                   summary_table.manMT(cond_unpp), ...
                                   summary_table.iniAngDevi(cond_unpp) };
            tmp_table = [tmp_table; tmp_row];
        end
    end
end
trial_var_names = { 'SubID', 'Section', 'PertProb', 'GoDelay', 'RT', ...
                    'pVel_all', ...
                    't2pVel_all', 'MT_all', ...
                    'pVel_np', ...
                    't2pVel_np', ...
                    'MT_np', ...
                    'iniAngDevi_np' };
                
trial_variables = array2table(tmp_table, 'VariableNames', trial_var_names);

%% plot for section 1
figure(1)
avg_var = zeros(3, 8);
sd_var = avg_var;
for i = subid
    cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
             cell2mat(trial_variables{:, 'Section'}) == 1 );
    avg_var(i, :) = cellfun(@mean, trial_variables{cond, 5:12});
    sd_var(i, :) = cellfun(@std, trial_variables{cond, 5:12});
end
for j = 1:size(avg_var, 2)
    subplot(2, 4, j)
    errorbar(avg_var(:, j), sd_var(:, j), 'o')
    ylabel(trial_var_names{4 + j});
end
savefig('figure1_baseline.fig')
%% plot for section 2 pertprob
figure(2)
avg_var = [];
sd_var = [];
real_picked = cell(3, 8);
for i = subid
    for pp = pert_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 2 & ...
                 cell2mat(trial_variables{:, 'PertProb'}) == pp );
        % merge trials
        tmp_picked = trial_variables{cond, 5:12};
        for c = 1:size(tmp_picked, 2)
            real_picked{i, c} = cell2mat(tmp_picked(:, c));
        end
        avg_var = [avg_var; i, pp, cellfun(@mean, real_picked(i, :))];
        sd_var = [sd_var; i, pp, cellfun(@std, real_picked(i, :))];
    end
end
avg_var = sortrows(avg_var, 2);
sd_var = sortrows(sd_var, 2);

for j = 3:size(avg_var, 2)
    subplot(2, 4, j-2)
    hold on
    for i = 1:3:size(avg_var, 1)
        errorbar((avg_var(i + 0, 2) - 2), avg_var(i + 0, j), sd_var(i + 0, j), 'or')
        errorbar((avg_var(i + 1, 2) + 0), avg_var(i + 1, j), sd_var(i + 1, j), '*b')
        errorbar((avg_var(i + 2, 2) + 2), avg_var(i + 2, j), sd_var(i + 2, j), 'xk')
    end
    ylabel(trial_var_names{4 + j - 2});
    hold off
end
savefig('figure2_pertb2.fig')
%% plot for section 2 godelay
figure(3)
avg_var = [];
sd_var = [];
real_picked = cell(3, 8);
for i = subid
    for gd = godelay_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 2 & ...
                 cell2mat(trial_variables{:, 'GoDelay'}) == gd );
        % merge trials
        tmp_picked = trial_variables{cond, 5:12};
        for c = 1:size(tmp_picked, 2)
            real_picked{i, c} = cell2mat(tmp_picked(:, c));
        end
        avg_var = [avg_var; i, gd, cellfun(@mean, real_picked(i, :))];
        sd_var = [sd_var; i, gd, cellfun(@std, real_picked(i, :))];
    end
end
avg_var = sortrows(avg_var, 2);
sd_var = sortrows(sd_var, 2);

for j = 3:size(avg_var, 2)
    subplot(2, 4, j-2)
    hold on
    for i = 1:3:size(avg_var, 1)
        errorbar((avg_var(i + 0, 2) - 20), avg_var(i + 0, j), sd_var(i + 0, j), 'or')
        errorbar((avg_var(i + 1, 2) + 00), avg_var(i + 1, j), sd_var(i + 1, j), '*b')
        errorbar((avg_var(i + 2, 2) + 20), avg_var(i + 2, j), sd_var(i + 2, j), 'xk')
    end
    ylabel(trial_var_names{4 + j - 2});
    hold off
end
savefig('figure3_godelay2.fig')
%% plot for section 3 pertprob
figure(4)
avg_var = [];
sd_var = [];
real_picked = cell(3, 8);
for i = subid
    for pp = pert_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 3 & ...
                 cell2mat(trial_variables{:, 'PertProb'}) == pp );
        % merge trials
        tmp_picked = trial_variables{cond, 5:12};
        for c = 1:size(tmp_picked, 2)
            real_picked{i, c} = cell2mat(tmp_picked(:, c));
        end
        avg_var = [avg_var; i, pp, cellfun(@mean, real_picked(i, :))];
        sd_var = [sd_var; i, pp, cellfun(@std, real_picked(i, :))];
    end
end
avg_var = sortrows(avg_var, 2);
sd_var = sortrows(sd_var, 2);

for j = 3:size(avg_var, 2)
    subplot(2, 4, j-2)
    hold on
    for i = 1:3:size(avg_var, 1)
        errorbar((avg_var(i + 0, 2) - 2), avg_var(i + 0, j), sd_var(i + 0, j), 'or')
        errorbar((avg_var(i + 1, 2) + 0), avg_var(i + 1, j), sd_var(i + 1, j), '*b')
        errorbar((avg_var(i + 2, 2) + 2), avg_var(i + 2, j), sd_var(i + 2, j), 'xk')
    end
    ylabel(trial_var_names{4 + j - 2});
    hold off
end
savefig('figure4_pertb3.fig')
%% plot for section 3 godelay
figure(5)
avg_var = [];
sd_var = [];
real_picked = cell(3, 8);
for i = subid
    for gd = godelay_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 3 & ...
                 cell2mat(trial_variables{:, 'GoDelay'}) == gd );
        % merge trials
        tmp_picked = trial_variables{cond, 5:12};
        for c = 1:size(tmp_picked, 2)
            real_picked{i, c} = cell2mat(tmp_picked(:, c));
        end
        avg_var = [avg_var; i, gd, cellfun(@mean, real_picked(i, :))];
        sd_var = [sd_var; i, gd, cellfun(@std, real_picked(i, :))];
    end
end
avg_var = sortrows(avg_var, 2);
sd_var = sortrows(sd_var, 2);

for j = 3:size(avg_var, 2)
    subplot(2, 4, j-2)
    hold on
    for i = 1:3:size(avg_var, 1)
        errorbar((avg_var(i + 0, 2) - 20), avg_var(i + 0, j), sd_var(i + 0, j), 'or')
        errorbar((avg_var(i + 1, 2) + 00), avg_var(i + 1, j), sd_var(i + 1, j), '*b')
        errorbar((avg_var(i + 2, 2) + 20), avg_var(i + 2, j), sd_var(i + 2, j), 'xk')
    end
    ylabel(trial_var_names{4 + j - 2});
    hold off
end

savefig('figure5_godelay3.fig')

%%
%{
trial_variables.Properties.VariableNames = { 'SubID', 'Section', 'RT', ...
                                             'pVel_all_avg', 'pVel_all_sd', ...
                                             't2pVel_all_avg', 't2pVel_all_sd', ...
                                             'MT_all_avg', 'MT_all_sd', ...
                                             'pVel_np_avg','pVel_np_sd', ...
                                             't2pVel_np_avg', 't2pVel_np_sd', ...
                                             'MT_np_avg', 'MT_np_sd', ...
                                             'iniAngDevi_np_avg', 'iniAngDevi_np_sd' };


trials = array2table(trials);
trials.Properties.VariableNames = { 'Time', 'TPID', 'Repetition', 'order',...
                                      'TargetLocation', 'PertProb', 'GoCueDelay',...
                                      'TrialFailed', 'Score',...
                                      'EventTime', 'KinARM' };                    
                    
                    
                    
                    

for i  = 1:length(subid)
RT{i, 1} = summary_table.manRT( trial_subid == 1 & trial_section == 3 & ...
                               ind_target_location & ...
                               trial_pert_lvl == 50 & ...
                               trial_godelay_lvl == 150 & ...
                               trial_failed == 0 & ...
                               trial_perturbed == 1);
end
%}





