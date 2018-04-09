%%
close all; clear; clc;
load('/Users/yen-hsunwu/Dropbox (ASU)/summary_table.mat')

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
trial_no_response = ( summary_table{:, 'manMoveOnset'} == 0 );

%% Get section III of all subjects
% exclude trials with no response for individual subject.
% When a trial is failed due to not response, the manMoveOnset equals
% to zero!

tmp_vars_traj = {'events', 'pos', 'vel', 'acc', 'force'};
tmp_vars = {'CueRT', 'manCueRT', 'RT', 'manRT', 'manMT', 'pVel', 'man_t2pVel', 'pAcc', 'man_t2pAcc', 'pForce', 'man_t2pForce', 'maxPathDevi', 'pathLength', 'iniAngDevi'};
sub_array_traj = cell(length(subid), length(section));
sub_array = cell(length(subid), length(section));
for i = subid
    for s = section
        if s ~= 1
        
            tmp_cond_traj = cell(length(godelay_lvl), length(pert_lvl));
            tmp_cond = cell(length(godelay_lvl), length(pert_lvl));
            for pp = 1:length(pert_lvl) % pertb has 3 level
                for gd = 1:length(godelay_lvl) % go delay has 4 level

                    cond = ( trial_subid == subid(i) & trial_section == s & ...
                             trial_pert_lvl == pert_lvl(pp) & ...
                             trial_godelay_lvl == godelay_lvl(gd)& ...
                             ~trial_no_response );
                    cond_unpb_all = (cond & ~trial_perturbed);
                    cond_pbed_all = (cond & trial_perturbed);
                    cond_unpb = (cond & ~trial_perturbed & ~trial_failed);
                    cond_pbed = (cond & trial_perturbed & ~trial_failed);


                    tmp_cond_traj{pp, gd} = { s, pert_lvl(pp), godelay_lvl(gd), summary_table(cond, tmp_vars_traj),...
                                                                                summary_table(cond_pbed_all, tmp_vars_traj),...
                                                                                summary_table(cond_unpb_all, tmp_vars_traj),...
                                                                                summary_table(cond_pbed, tmp_vars_traj),...
                                                                                summary_table(cond_unpb, tmp_vars_traj) };
                                                                       
                    tmp_cond{pp, gd} = { s, pert_lvl(pp), godelay_lvl(gd), summary_table(cond, tmp_vars),...
                                                                           summary_table(cond_pbed_all, tmp_vars),...
                                                                           summary_table(cond_unpb_all, tmp_vars),...
                                                                           summary_table(cond_pbed, tmp_vars),...
                                                                           summary_table(cond_unpb, tmp_vars) };
                end
            end
        else
            cond = ( trial_subid == subid(i) & trial_section == 1 & ~trial_no_response);
            cond_success = ( cond & ~trial_failed );
            tmp_cond_traj = { s, summary_table(cond, tmp_vars_traj), summary_table(cond_success, tmp_vars_traj) };
            tmp_cond = { s, summary_table(cond, tmp_vars), summary_table(cond_success, tmp_vars) };
        end
        
        sub_array_traj{i, s} = tmp_cond_traj;
        sub_array{i, s} = tmp_cond;
    end
end

%% align trajectory data at go cue (ind_go_on), normalize data [ind_move_onset, ind_move_end]
% and compute mean trajectories



%% compute mean variables
sub_avg = cell(length(subid), length(section));
for i = subid
    for s = section        
        if s ~= 1
            tmp_cond = cell(length(godelay_lvl), length(pert_lvl));
            for pp = 1:length(pert_lvl) % pertb has 3 level
                for gd = 1:length(godelay_lvl) % go delay has 4 level
                    tmp_avg = cell(1, 5);
                    for k = 1:5
                        tmp_table = sub_array{i, s}{pp, gd}{1, 3 + k};
                        tmp_avg{1, k} = array2table(nanmean(table2array(tmp_table), 1), 'VariableNames', tmp_table.Properties.VariableNames);
                    end
                    
                    tmp_cond{pp, gd} = [s, pert_lvl(pp), godelay_lvl(gd), tmp_avg];
                end
            end
        else
            tmp_avg = cell(1, 2);
            for k = 1:2
                tmp_table = sub_array{i, s}{1, 1 + k};
                tmp_avg{1, k} = array2table(nanmean(table2array(tmp_table), 1), 'VariableNames', tmp_table.Properties.VariableNames);
            end
            tmp_cond = [s, tmp_avg];
        end
        
        sub_avg{i, s} = tmp_cond;
    end
end

%% reshape each condition into individual variables
sub_reshape = cell(length(subid), 1);
for i = subid
    s = 3;

    tmp_var = cell(1, length(tmp_vars));
    for v = 1:14
        tmp_pick = cell(1, 5); % 5 is length( [all, unpp_all, pped_all, unpp, pped] )
        
        for k = 1:(length(tmp_conds) - 3)
            tmp_table = cell(length(godelay_lvl), length(pert_lvl));
            for pp = 1:length(pert_lvl) % pertb has 3 level
                for gd = 1:length(godelay_lvl) % go delay has 4 level
                    tmp_conds = sub_avg{i, s}{pp, gd};
                
                    tmp_table{pp, gd} = tmp_conds{1, 3 + k}{1, v};
                    tmp_pick{1, k} = tmp_table;
                end
            end
        end
        tmp_var{1, v} = cell2table(tmp_pick, 'VariableNames', {'all', 'pb_all', 'unpb_all', 'pb', 'unpb'});
    end

    sub_reshape{i, 1} = cell2table(tmp_var, 'VariableNames', tmp_vars);
end

%% plot bar plots
sub_no = length(sub_reshape);
xlebels = {'80ms', '150ms', '220ms', '300ms'};
legends = {'pp75%', 'pp50%', 'pp25%'};
colormap('gray')
for i = 1:sub_no
    h(1) = figure(1);
    subplot(sub_no, 2, (2 * i - 1))
    bar(cell2mat(sub_reshape{i, 1}{1, 'CueRT'}{1, 'all'}{:})', 'FaceColor', 'flat')
    ylim([.5, 1]); xticklabels(xlebels); ylabel('CueRT (s)');
    subplot(sub_no, 2, (2 * i))
    bar(cell2mat(sub_reshape{i, 1}{1, 'RT'}{1, 'all'}{:})')
    ylim([0, .5]); xticklabels(xlebels); ylabel('RT (s)');
    
    h(2) = figure(2);
    subplot(sub_no, 2, (2 * i - 1))
    bar(cell2mat(sub_reshape{i, 1}{1, 'manCueRT'}{1, 'all'}{:})')
    ylim([.5, 1]); xticklabels(xlebels); ylabel('manual CueRT (s)');
    subplot(sub_no, 2, (2 * i))
    bar(cell2mat(sub_reshape{i, 1}{1, 'manRT'}{1, 'all'}{:})')
    ylim([0, .5]); xticklabels(xlebels); ylabel('manual (s)');
    
    h(3) = figure(3);
    subplot(sub_no, 2, (2 * i - 1))
    bar(cell2mat(sub_reshape{i, 1}{1, 'manMT'}{1, 'all'}{:})')
    ylim([.5, 1]); xticklabels(xlebels); ylabel('MT (s)');
    subplot(sub_no, 2, (2 * i))
    bar(cell2mat(sub_reshape{i, 1}{1, 'iniAngDevi'}{1, 'all'}{:})')
    ylim([-8, 4]); xticklabels(xlebels); ylabel('Ini Ang Devi (deg)');
    
    h(4) = figure(4);
    subplot(sub_no, 2, (2 * i - 1))
    bar(cell2mat(sub_reshape{i, 1}{1, 'man_t2pVel'}{1, 'all'}{:})')
    ylim([0, .5]); xticklabels(xlebels); ylabel('t2pVel (s)');
    subplot(sub_no, 2, (2 * i))
    bar(cell2mat(sub_reshape{i, 1}{1, 'pVel'}{1, 'all'}{:})')
    ylim([0, .5]); xticklabels(xlebels); ylabel('pVel (m/s)');
    
    h(5) = figure(5);
    subplot(sub_no, 2, (2 * i - 1))
    bar(cell2mat(sub_reshape{i, 1}{1, 'man_t2pAcc'}{1, 'all'}{:})')
    ylim([0, .5]); xticklabels(xlebels); ylabel('t2pAcc (s)');
    subplot(sub_no, 2, (2 * i))
    bar(cell2mat(sub_reshape{i, 1}{1, 'pAcc'}{1, 'all'}{:})')
    ylim([0, .005]); xticklabels(xlebels); ylabel('pAcc (m^2/s)');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h(6) = figure(6);
    subplot(sub_no, 2, (2 * i - 1))
    bar(cell2mat(sub_reshape{i, 1}{1, 'man_t2pForce'}{1, 'pb'}{:})')
    ylim([0, .5]); xticklabels(xlebels); ylabel('pb t2pForce (s)');
    subplot(sub_no, 2, (2 * i))
    bar(cell2mat(sub_reshape{i, 1}{1, 'pForce'}{1, 'pb'}{:})')
    ylim([0, 8]); xticklabels(xlebels); ylabel('pb pFroce (N)');

    h(7) = figure(7);
    subplot(sub_no, 2, (2 * i - 1))
    bar(cell2mat(sub_reshape{i, 1}{1, 'maxPathDevi'}{1, 'pb'}{:})')
    ylim([0, .05]); xticklabels(xlebels); ylabel('pb max devi (m)');
    subplot(sub_no, 2, (2 * i))
    bar(cell2mat(sub_reshape{i, 1}{1, 'maxPathDevi'}{1, 'unpb'}{:})')
    ylim([0, .05]); xticklabels(xlebels); ylabel('unpb max devi (m)'); legend(legends)
end

for i = 1:7
    savefig(h(i), ['figure', num2str(i), '.fig']);
end








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
% When a trial is failed due to not response, the manMoveOnset equals
% to zero!

%{
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
                    cond_pped = (cond & trial_perturbed & ~trial_failed);
                    
                    tmp_row = { i, s, pert_lvl(pp), godelay_lvl(gd), summary_table.manRT(cond_pped), ...
                                                                     summary_table.pVel(cond), ...
                                                                     summary_table.man_t2pVel(cond), ...
                                                                     summary_table.manMT(cond), ...
                                                                     summary_table.pVel(cond_unpp), ...
                                                                     summary_table.man_t2pVel(cond_unpp), ...
                                                                     summary_table.manMT(cond_unpp), ...
                                                                     summary_table.iniAngDevi(cond_unpp)...
                                                                     summary_table.maxPathDevi(cond_pped), ...
                                                                     summary_table.pathLength(cond_pped) };

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
                                   summary_table.iniAngDevi(cond_unpp), ...
                                   summary_table.maxPathDevi(cond_unpp), ...
                                   summary_table.pathLength(cond_unpp) };
            tmp_table = [tmp_table; tmp_row];
        end
    end
end
trial_var_names = { 'SubID', 'Section', 'PertProb', 'GoDelay', 'RT_p', ...
                    'pVel_all', ...
                    't2pVel_all', 'MT_all', ...
                    'pVel_np', ...
                    't2pVel_np', ...
                    'MT_np', ...
                    'iniAngDevi_np', 'maxPathDevi_p', 'pathLength_p' };
                
trial_variables = array2table(tmp_table, 'VariableNames', trial_var_names);
%}





%% correlation analysis: scatter plot, Pearson corr, Spearman corr
avg_var = [];
sd_var = [];
raw_var = [];
real_picked = cell(3, 10);
for i = subid
    for gd = godelay_lvl
        for pp = pert_lvl
            cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                     cell2mat(trial_variables{:, 'Section'}) == 3 & ...
                     cell2mat(trial_variables{:, 'PertProb'}) == pp & ...
                     cell2mat(trial_variables{:, 'GoDelay'}) == gd );
            % merge trials
            tmp_picked = trial_variables{cond, 5:14};
            for c = 1:size(tmp_picked, 2)
                real_picked{i, c} = cell2mat(tmp_picked(:, c));
            end
            raw_var = [raw_var; i, pp, gd, real_picked(i, :)];
            avg_var = [avg_var; i, pp, gd, cellfun(@mean, real_picked(i, :))];
            sd_var = [sd_var; i, pp, gd, cellfun(@std, real_picked(i, :))];
        end
    end
end

tmp_var_names = { 'SubID', 'PertProb', 'GoDelay', 'RT', ...
                  'pVel_all', ...
                  't2pVel_all', 'MT_all', ...
                  'pVel_np', ...
                  't2pVel_np', ...
                  'MT_np', ...
                  'iniAngDevi_np', 'maxPathDevi_p', 'pathLength_p' };

raw_var = array2table(raw_var, 'VariableNames', tmp_var_names);
avg_var = array2table(avg_var, 'VariableNames', tmp_var_names);
sd_var = array2table(sd_var, 'VariableNames', tmp_var_names);

raw_4corr = raw_var(:, {'SubID', 'PertProb', 'GoDelay', 'RT', 'maxPathDevi_p'});
avg_4corr = avg_var(:, {'SubID', 'PertProb', 'GoDelay', 'RT', 'maxPathDevi_p'});
sd_4corr = sd_var(:, {'SubID', 'PertProb', 'GoDelay', 'RT', 'maxPathDevi_p'});

raw_4corr = sortrows(raw_4corr, {'PertProb', 'GoDelay'});
avg_4corr = sortrows(avg_4corr, {'PertProb', 'GoDelay'});
sd_4corr = sortrows(sd_4corr, {'PertProb', 'GoDelay'});

r = zeros(height(raw_4corr), 5);
for i = 1:height(raw_4corr)
    r(i ,1) = raw_4corr{i, 'SubID'}{:};
    r(i ,2) = raw_4corr{i, 'PertProb'}{:};
    r(i, 3) = raw_4corr{i, 'GoDelay'}{:};
    [r(i, 4), ~] = corr(raw_4corr{i, 'RT'}{:}, raw_4corr{i, 'maxPathDevi_p'}{:},'type', 'Pearson');
    [r(i, 5), ~] = corr(raw_4corr{i, 'RT'}{:}, raw_4corr{i, 'maxPathDevi_p'}{:},'type', 'Spearman');  
end
tmp_var_names = { 'SubID', 'PertProb', 'GoDelay', 'PearsonRT2Err', 'SpearmanRT2Err' };
r = array2table(r, 'VariableNames', tmp_var_names);

r = sortrows(r, {'SubID', 'PertProb', 'GoDelay'});

for i = 1:height(raw_4corr)
    scatter(raw_4corr{i, 'RT'}{:}, raw_4corr{i, 'maxPathDevi_p'}{:});
    pause
end


% % % 
% % % r = zeros(height(avg_4corr)/3, 4);
% % % j = 1;
% % % for i = 1:3:height(avg_4corr)
% % %     r(j ,1) = avg_4corr{i, 'PertProb'};
% % %     r(j, 2) = avg_4corr{i, 'GoDelay'};
% % %     [r(j, 3), ~] = corr(avg_4corr{i:i+2, 'RT'}, avg_4corr{i:i+2, 'maxPathDevi_p'},'type', 'Pearson');
% % %     [r(j, 4), ~] = corr(avg_4corr{i:i+2, 'RT'}, avg_4corr{i:i+2, 'maxPathDevi_p'},'type', 'Spearman');
% % %     j = j + 1;
% % % end
% % % tmp_var_names = { 'PertProb', 'GoDelay', 'PearsonRT2Err', 'SpearmanRT2Err' };
% % % r = array2table(r, 'VariableNames', tmp_var_names);








%% plot for section 1
figure(1)
avg_var = zeros(3, 10);
sd_var = avg_var;
for i = subid
    cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
             cell2mat(trial_variables{:, 'Section'}) == 1 );
    avg_var(i, :) = cellfun(@mean, trial_variables{cond, 5:14});
    sd_var(i, :) = cellfun(@std, trial_variables{cond, 5:14});
end
for j = 1:size(avg_var, 2)
    subplot(2, 5, j)
    errorbar(avg_var(:, j), sd_var(:, j), 'o')
    ylabel(trial_var_names{4 + j});
end
savefig('figure1_baseline.fig')
%% plot for section 2 pertprob
figure(2)
avg_var = [];
sd_var = [];
real_picked = cell(3, 10);
for i = subid
    for pp = pert_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 2 & ...
                 cell2mat(trial_variables{:, 'PertProb'}) == pp );
        % merge trials
        tmp_picked = trial_variables{cond, 5:14};
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
    subplot(2, 5, j-2)
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
real_picked = cell(3, 10);
for i = subid
    for gd = godelay_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 2 & ...
                 cell2mat(trial_variables{:, 'GoDelay'}) == gd );
        % merge trials
        tmp_picked = trial_variables{cond, 5:14};
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
    subplot(2, 5, j-2)
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
real_picked = cell(3, 10);
for i = subid
    for pp = pert_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 3 & ...
                 cell2mat(trial_variables{:, 'PertProb'}) == pp );
        % merge trials
        tmp_picked = trial_variables{cond, 5:14};
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
    subplot(2, 5, j-2)
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
real_picked = cell(3, 10);
for i = subid
    for gd = godelay_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 3 & ...
                 cell2mat(trial_variables{:, 'GoDelay'}) == gd );
        % merge trials
        tmp_picked = trial_variables{cond, 5:14};
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
    subplot(2, 5, j-2)
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

%% plot for section 2 & 3 godelay
figure(6)
avg_var = [];
sd_var = [];
real_picked = cell(3, 10);
for i = subid
    for gd = godelay_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 2 & ...
                 cell2mat(trial_variables{:, 'GoDelay'}) == gd );
        % merge trials
        tmp_picked = trial_variables{cond, 5:14};
        for c = 1:size(tmp_picked, 2)
            real_picked{i, c} = cell2mat(tmp_picked(:, c));
        end
        avg_var = [avg_var; i, gd, cellfun(@mean, real_picked(i, :))];
        sd_var = [sd_var; i, gd, cellfun(@std, real_picked(i, :))];
    end
end
avg_var = sortrows(avg_var, 2);
sd_var = sortrows(sd_var, 2);

subplot(2, 3, 1)
hold on
for i = 1:3:size(avg_var, 1)
    errorbar((avg_var(i + 0, 2) - 10), avg_var(i + 0, 3), sd_var(i + 0, 3), 'or')
    errorbar((avg_var(i + 1, 2) + 00), avg_var(i + 1, 3), sd_var(i + 1, 3), '*b')
    errorbar((avg_var(i + 2, 2) + 10), avg_var(i + 2, 3), sd_var(i + 2, 3), 'xk')
end
ylabel([trial_var_names{5}, ' (s)']);
xlabel('Go delay (ms)');
hold off

subplot(2, 3, 4)
hold on
for i = 1:3:size(avg_var, 1)
    errorbar((avg_var(i + 0, 2) - 10), avg_var(i + 0, 11) * 100, sd_var(i + 0, 11) * 100, 'or')
    errorbar((avg_var(i + 1, 2) + 00), avg_var(i + 1, 11) * 100, sd_var(i + 1, 11) * 100, '*b')
    errorbar((avg_var(i + 2, 2) + 10), avg_var(i + 2, 11) * 100, sd_var(i + 2, 11) * 100, 'xk')
end
ylabel('maxPathDevi perturbed (cm)');
xlabel('Go delay (ms)');
hold off

% plot for section 3 godelay
figure(6)
avg_var = [];
sd_var = [];
real_picked = cell(3, 10);
for i = subid
    for gd = godelay_lvl
            cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                     cell2mat(trial_variables{:, 'Section'}) == 3 & ...
                     cell2mat(trial_variables{:, 'GoDelay'}) == gd );
            % merge trials
            tmp_picked = trial_variables{cond, 5:14};
            for c = 1:size(tmp_picked, 2)
                real_picked{i, c} = cell2mat(tmp_picked(:, c));
            end
            avg_var = [avg_var; i, gd, cellfun(@mean, real_picked(i, :))];
            sd_var = [sd_var; i, gd, cellfun(@std, real_picked(i, :))];
    end
end
avg_var = sortrows(avg_var, 2);
sd_var = sortrows(sd_var, 2);

subplot(2, 3, 3)
hold on
for i = 1:3:size(avg_var, 1)
    errorbar((avg_var(i + 0, 2) - 10), avg_var(i + 0, 3), sd_var(i + 0, 3), 'or')
    errorbar((avg_var(i + 1, 2) + 00), avg_var(i + 1, 3), sd_var(i + 1, 3), '*b')
    errorbar((avg_var(i + 2, 2) + 10), avg_var(i + 2, 3), sd_var(i + 2, 3), 'xk')
end
ylabel([trial_var_names{5}, ' (s)']);
xlabel('Go delay (ms)');
hold off

subplot(2, 3, 6)
hold on
for i = 1:3:size(avg_var, 1)
    errorbar((avg_var(i + 0, 2) - 10), avg_var(i + 0, 11) * 100, sd_var(i + 0, 11) * 100, 'or')
    errorbar((avg_var(i + 1, 2) + 00), avg_var(i + 1, 11) * 100, sd_var(i + 1, 11) * 100, '*b')
    errorbar((avg_var(i + 2, 2) + 10), avg_var(i + 2, 11) * 100, sd_var(i + 2, 11) * 100, 'xk')
end
ylabel('maxPathDevi perturbed (cm)');
xlabel('Go delay (ms)');
hold off

savefig('figure6_godelay.fig')

%% plot for section 2 & 3 pertprob
figure(7)
avg_var = [];
sd_var = [];
real_picked = cell(3, 10);
for i = subid
    for pp = pert_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 2 & ...
                 cell2mat(trial_variables{:, 'PertProb'}) == pp );
        % merge trials
        tmp_picked = trial_variables{cond, 5:14};
        for c = 1:size(tmp_picked, 2)
            real_picked{i, c} = cell2mat(tmp_picked(:, c));
        end
        avg_var = [avg_var; i, pp, cellfun(@mean, real_picked(i, :))];
        sd_var = [sd_var; i, pp, cellfun(@std, real_picked(i, :))];
    end
end
avg_var = sortrows(avg_var, 2);
sd_var = sortrows(sd_var, 2);

subplot(2, 3, 1)
hold on
for i = 1:3:size(avg_var, 1)
    errorbar((avg_var(i + 0, 2) - 5), avg_var(i + 0, 3), sd_var(i + 0, 3), 'or')
    errorbar((avg_var(i + 1, 2) + 0), avg_var(i + 1, 3), sd_var(i + 1, 3), '*b')
    errorbar((avg_var(i + 2, 2) + 5), avg_var(i + 2, 3), sd_var(i + 2, 3), 'xk')
end
ylabel([trial_var_names{5}, ' (s)']);
xlabel('PP (%)');
ylim([-0.2, 0.6])
hold off

subplot(2, 3, 4)
hold on
for i = 1:3:size(avg_var, 1)
    errorbar((avg_var(i + 0, 2) - 5), avg_var(i + 0, 11) * 100, sd_var(i + 0, 11) * 100, 'or')
    errorbar((avg_var(i + 1, 2) + 0), avg_var(i + 1, 11) * 100, sd_var(i + 1, 11) * 100, '*b')
    errorbar((avg_var(i + 2, 2) + 5), avg_var(i + 2, 11) * 100, sd_var(i + 2, 11) * 100, 'xk')
end
ylabel('maxPathDevi perturbed (cm)');
xlabel('PP (%)');

avg_var = [];
sd_var = [];
real_picked = cell(3, 10);
for i = subid
    for pp = pert_lvl
        cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                 cell2mat(trial_variables{:, 'Section'}) == 3 & ...
                 cell2mat(trial_variables{:, 'PertProb'}) == pp );
        % merge trials
        tmp_picked = trial_variables{cond, 5:14};
        for c = 1:size(tmp_picked, 2)
            real_picked{i, c} = cell2mat(tmp_picked(:, c));
        end
        avg_var = [avg_var; i, pp, cellfun(@mean, real_picked(i, :))];
        sd_var = [sd_var; i, pp, cellfun(@std, real_picked(i, :))];
    end
end
avg_var = sortrows(avg_var, 2);
sd_var = sortrows(sd_var, 2);

subplot(2, 3, 3)
hold on
for i = 1:3:size(avg_var, 1)
    errorbar((avg_var(i + 0, 2) - 5), avg_var(i + 0, 3), sd_var(i + 0, 3), 'or')
    errorbar((avg_var(i + 1, 2) + 0), avg_var(i + 1, 3), sd_var(i + 1, 3), '*b')
    errorbar((avg_var(i + 2, 2) + 5), avg_var(i + 2, 3), sd_var(i + 2, 3), 'xk')
end
ylabel([trial_var_names{5}, ' (s)']);
xlabel('PP (%)');
ylim([-0.2, 0.6])
hold off

subplot(2, 3, 6)
hold on
for i = 1:3:size(avg_var, 1)
    errorbar((avg_var(i + 0, 2) - 5), avg_var(i + 0, 11) * 100, sd_var(i + 0, 11) * 100, 'or')
    errorbar((avg_var(i + 1, 2) + 0), avg_var(i + 1, 11) * 100, sd_var(i + 1, 11) * 100, '*b')
    errorbar((avg_var(i + 2, 2) + 5), avg_var(i + 2, 11) * 100, sd_var(i + 2, 11) * 100, 'xk')
end
ylabel('maxPathDevi perturbed (cm)');
xlabel('PP (%)');
hold off
savefig('figure7_pertb.fig')


%% plot for section 3 pertprob x godelay
figure(8)
avg_var = [];
sd_var = [];
real_picked = cell(3, 10);
for i = subid
    for gd = godelay_lvl
        for pp = pert_lvl
            cond = ( cell2mat(trial_variables{:, 'SubID'}) == subid(i) & ...
                     cell2mat(trial_variables{:, 'Section'}) == 3 & ...
                     cell2mat(trial_variables{:, 'PertProb'}) == pp & ...
                     cell2mat(trial_variables{:, 'GoDelay'}) == gd );
            % merge trials
            tmp_picked = trial_variables{cond, 5:14};
            for c = 1:size(tmp_picked, 2)
                real_picked{i, c} = cell2mat(tmp_picked(:, c));
            end
            avg_var = [avg_var; i, pp, gd, cellfun(@mean, real_picked(i, :))];
            sd_var = [sd_var; i, pp, gd, cellfun(@std, real_picked(i, :))];
        end
    end
end
avg_var = sortrows(avg_var, [2, 3]);
sd_var = sortrows(sd_var, [2, 3]);



x = {'P(pert) 25%', 'P(pert) 50%', 'P(pert) 75%'};
% x = {'Go delay 80 ms', 'Go delay 150 ms', 'Go delay 220 ms', 'Go delay 300 ms'};
offset = 10;
condition = 3;
for j = 1:12:size(avg_var, 1)
    subplot(2, condition, floor(j / 12) + 1)
    hold on
    for i = j:3:(j + 11)
        errorbar((avg_var(i + 0, 2) + avg_var(i + 0, 3) - offset), avg_var(i + 0, 4), sd_var(i + 0, 4), 'or')
        errorbar((avg_var(i + 1, 2) + avg_var(i + 1, 3) + 0     ), avg_var(i + 1, 4), sd_var(i + 1, 4), '*b')
        errorbar((avg_var(i + 2, 2) + avg_var(i + 2, 3) + offset), avg_var(i + 2, 4), sd_var(i + 2, 4), 'xk')
    end
    ylabel([trial_var_names{5}, ' (s)']);
    xlabel(x{floor(j / 12) + 1});
    ylim([-0.2, 0.6]);
    hold off
    
    subplot(2, condition, floor(j / 12) + 1 + condition)
    hold on
    for i = j:3:(j + 11)
        errorbar((avg_var(i + 0, 2) + avg_var(i + 0, 3) - offset), avg_var(i + 0, 12) * 100, sd_var(i + 0, 12) * 100, 'or')
        errorbar((avg_var(i + 1, 2) + avg_var(i + 1, 3) + 0     ), avg_var(i + 1, 12) * 100, sd_var(i + 1, 12) * 100, '*b')
        errorbar((avg_var(i + 2, 2) + avg_var(i + 2, 3) + offset), avg_var(i + 2, 12) * 100, sd_var(i + 2, 12) * 100, 'xk')
    end
    ylabel('maxPathDevi perturbed (cm)');
    xlabel(x{floor(j / 12) + 1});
    ylim([-0.2, 3.5]);
    hold off
end

savefig('figure8_section3_pertb_x_godelay.fig')

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
