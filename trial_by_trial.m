close all; clear; clc;

%% Load .mat
[filename, filepath, ~] = uigetfile('summary.mat');
full_filename = fullfile(filepath, filename);
load(full_filename);






%% 
logic_M = ~cellfun(@isempty, strfind(trial_behavior{:, 'TargetLocation'}, 'M'));
logic_R = ~cellfun(@isempty, strfind(trial_behavior{:, 'TargetLocation'}, 'R'));
logic_L = ~cellfun(@isempty, strfind(trial_behavior{:, 'TargetLocation'}, 'L'));
logic_25 = (cell2mat(trial_behavior{:, 'PertProb'}) == 25);
logic_50 = (cell2mat(trial_behavior{:, 'PertProb'}) == 50);
logic_75 = (cell2mat(trial_behavior{:, 'PertProb'}) == 75);


trial_behavior_M = trial_behavior(logic_M, :);
trial_behavior_M_25 = trial_behavior( all([logic_M, logic_25], 2), : );
trial_behavior_M_50 = trial_behavior( all([logic_M, logic_50], 2), : );
trial_behavior_M_75 = trial_behavior( all([logic_M, logic_75], 2), : );

%{
% trial_behavior_M_25_80 = trial_behavior( all([logic_M, logic_25, logic_80], 2), : );
% trial_behavior_M_25_150 = trial_behavior( all([logic_M, logic_25, logic_150], 2), : );
% trial_behavior_M_25_220 = trial_behavior( all([logic_M, logic_25, logic_220], 2), : );
% trial_behavior_M_25_300 = trial_behavior( all([logic_M, logic_25, logic_300], 2), : );
% trial_behavior_M_50_80 = trial_behavior( all([logic_M, logic_50, logic_80], 2), : );
% trial_behavior_M_50_150 = trial_behavior( all([logic_M, logic_50, logic_150], 2), : );
% trial_behavior_M_50_220 = trial_behavior( all([logic_M, logic_50, logic_220], 2), : );
% trial_behavior_M_50_300 = trial_behavior( all([logic_M, logic_50, logic_300], 2), : );
% trial_behavior_M_75_80 = trial_behavior( all([logic_M, logic_75, logic_80], 2), : );
% trial_behavior_M_75_150 = trial_behavior( all([logic_M, logic_75, logic_150], 2), : );
% trial_behavior_M_75_220 = trial_behavior( all([logic_M, logic_75, logic_220], 2), : );
% trial_behavior_M_75_300 = trial_behavior( all([logic_M, logic_75, logic_300], 2), : );
%}

trial_behavior_L = trial_behavior(logic_R, :);
trial_behavior_L_25 = trial_behavior( all([logic_L, logic_25], 2), : );
trial_behavior_L_50 = trial_behavior( all([logic_L, logic_50], 2), : );
trial_behavior_L_75 = trial_behavior( all([logic_L, logic_75], 2), : );
trial_behavior_R = trial_behavior(logic_L, :);
trial_behavior_R_25 = trial_behavior( all([logic_R, logic_25], 2), : );
trial_behavior_R_50 = trial_behavior( all([logic_R, logic_50], 2), : );
trial_behavior_R_75 = trial_behavior( all([logic_R, logic_75], 2), : );

%%
figure(2)
plot_variable = 'MT'; str_variable = sprintf('%s (s)', plot_variable);
subplot(3, 3, 2)
logic_failed = logical(cell2mat(trial_behavior_M_25{:, 'TrialFailed'}));
logic_80 = (cell2mat(trial_behavior_M_25{:, 'GoCueDelay'}) == 80);
logic_150 = (cell2mat(trial_behavior_M_25{:, 'GoCueDelay'}) == 150);
logic_220 = (cell2mat(trial_behavior_M_25{:, 'GoCueDelay'}) == 220);
logic_300 = (cell2mat(trial_behavior_M_25{:, 'GoCueDelay'}) == 300);

plot( find( all([logic_failed, logic_80], 2) ), trial_behavior_M_25{ all([logic_failed, logic_80], 2), plot_variable }, 'or',...
      find( all([~logic_failed, logic_80], 2) ), trial_behavior_M_25{ all([~logic_failed, logic_80], 2), plot_variable }, 'xr',...
      find( all([logic_failed, logic_150], 2) ), trial_behavior_M_25{ all([logic_failed, logic_150], 2), plot_variable }, 'ob',...
      find( all([~logic_failed, logic_150], 2) ), trial_behavior_M_25{ all([~logic_failed, logic_150], 2), plot_variable }, 'xb',...
      find( all([logic_failed, logic_220], 2) ), trial_behavior_M_25{ all([logic_failed, logic_220], 2), plot_variable }, 'ok',...
      find( all([~logic_failed, logic_220], 2) ), trial_behavior_M_25{ all([~logic_failed, logic_220], 2), plot_variable }, 'xk',...
      find( all([logic_failed, logic_300], 2) ), trial_behavior_M_25{ all([logic_failed, logic_300], 2), plot_variable }, 'om',...
      find( all([~logic_failed, logic_300], 2) ), trial_behavior_M_25{ all([~logic_failed, logic_300], 2), plot_variable }, 'xm' );
ylabel(str_variable); xlabel('Middle Target, 25%');
ylim([0, 0.8]);
subplot(3, 3, 5)
logic_failed = logical(cell2mat(trial_behavior_M_50{:, 'TrialFailed'}));
logic_80 = (cell2mat(trial_behavior_M_50{:, 'GoCueDelay'}) == 80);
logic_150 = (cell2mat(trial_behavior_M_50{:, 'GoCueDelay'}) == 150);
logic_220 = (cell2mat(trial_behavior_M_50{:, 'GoCueDelay'}) == 220);
logic_300 = (cell2mat(trial_behavior_M_50{:, 'GoCueDelay'}) == 300);

plot( find( all([logic_failed, logic_80], 2) ), trial_behavior_M_50{ all([logic_failed, logic_80], 2), plot_variable }, 'or',...
      find( all([~logic_failed, logic_80], 2) ), trial_behavior_M_50{ all([~logic_failed, logic_80], 2), plot_variable }, 'xr',...
      find( all([logic_failed, logic_150], 2) ), trial_behavior_M_50{ all([logic_failed, logic_150], 2), plot_variable }, 'ob',...
      find( all([~logic_failed, logic_150], 2) ), trial_behavior_M_50{ all([~logic_failed, logic_150], 2), plot_variable }, 'xb',...
      find( all([logic_failed, logic_220], 2) ), trial_behavior_M_50{ all([logic_failed, logic_220], 2), plot_variable }, 'ok',...
      find( all([~logic_failed, logic_220], 2) ), trial_behavior_M_50{ all([~logic_failed, logic_220], 2), plot_variable }, 'xk',...
      find( all([logic_failed, logic_300], 2) ), trial_behavior_M_50{ all([logic_failed, logic_300], 2), plot_variable }, 'om',...
      find( all([~logic_failed, logic_300], 2) ), trial_behavior_M_50{ all([~logic_failed, logic_300], 2), plot_variable }, 'xm' );
ylabel(str_variable); xlabel('Middle Target, 50%');
ylim([0, 0.8]);

subplot(3, 3, 8)
logic_failed = logical(cell2mat(trial_behavior_M_75{:, 'TrialFailed'}));
logic_80 = (cell2mat(trial_behavior_M_75{:, 'GoCueDelay'}) == 80);
logic_150 = (cell2mat(trial_behavior_M_75{:, 'GoCueDelay'}) == 150);
logic_220 = (cell2mat(trial_behavior_M_75{:, 'GoCueDelay'}) == 220);
logic_300 = (cell2mat(trial_behavior_M_75{:, 'GoCueDelay'}) == 300);

plot( find( all([logic_failed, logic_80], 2) ), trial_behavior_M_75{ all([logic_failed, logic_80], 2), plot_variable }, 'or',...
      find( all([~logic_failed, logic_80], 2) ), trial_behavior_M_75{ all([~logic_failed, logic_80], 2), plot_variable }, 'xr',...
      find( all([logic_failed, logic_150], 2) ), trial_behavior_M_75{ all([logic_failed, logic_150], 2), plot_variable }, 'ob',...
      find( all([~logic_failed, logic_150], 2) ), trial_behavior_M_75{ all([~logic_failed, logic_150], 2), plot_variable }, 'xb',...
      find( all([logic_failed, logic_220], 2) ), trial_behavior_M_75{ all([logic_failed, logic_220], 2), plot_variable }, 'ok',...
      find( all([~logic_failed, logic_220], 2) ), trial_behavior_M_75{ all([~logic_failed, logic_220], 2), plot_variable }, 'xk',...
      find( all([logic_failed, logic_300], 2) ), trial_behavior_M_75{ all([logic_failed, logic_300], 2), plot_variable }, 'om',...
      find( all([~logic_failed, logic_300], 2) ), trial_behavior_M_75{ all([~logic_failed, logic_300], 2), plot_variable }, 'xm' );
ylabel(str_variable); xlabel('Middle Target, 70%');
ylim([0, 0.8]);

subplot(3, 3, 1)
logic_failed = logical(cell2mat(trial_behavior_L_25{:, 'TrialFailed'}));
logic_80 = (cell2mat(trial_behavior_L_25{:, 'GoCueDelay'}) == 80);
logic_150 = (cell2mat(trial_behavior_L_25{:, 'GoCueDelay'}) == 150);
logic_220 = (cell2mat(trial_behavior_L_25{:, 'GoCueDelay'}) == 220);
logic_300 = (cell2mat(trial_behavior_L_25{:, 'GoCueDelay'}) == 300);

plot( find( all([logic_failed, logic_80], 2) ), trial_behavior_L_25{ all([logic_failed, logic_80], 2), plot_variable }, 'or',...
      find( all([~logic_failed, logic_80], 2) ), trial_behavior_L_25{ all([~logic_failed, logic_80], 2), plot_variable }, 'xr',...
      find( all([logic_failed, logic_150], 2) ), trial_behavior_L_25{ all([logic_failed, logic_150], 2), plot_variable }, 'ob',...
      find( all([~logic_failed, logic_150], 2) ), trial_behavior_L_25{ all([~logic_failed, logic_150], 2), plot_variable }, 'xb',...
      find( all([logic_failed, logic_220], 2) ), trial_behavior_L_25{ all([logic_failed, logic_220], 2), plot_variable }, 'ok',...
      find( all([~logic_failed, logic_220], 2) ), trial_behavior_L_25{ all([~logic_failed, logic_220], 2), plot_variable }, 'xk',...
      find( all([logic_failed, logic_300], 2) ), trial_behavior_L_25{ all([logic_failed, logic_300], 2), plot_variable }, 'om',...
      find( all([~logic_failed, logic_300], 2) ), trial_behavior_L_25{ all([~logic_failed, logic_300], 2), plot_variable }, 'xm' );
ylabel(str_variable); xlabel('Left Target, 25%');
ylim([0, 0.8]);
subplot(3, 3, 4)
logic_failed = logical(cell2mat(trial_behavior_L_50{:, 'TrialFailed'}));
logic_80 = (cell2mat(trial_behavior_L_50{:, 'GoCueDelay'}) == 80);
logic_150 = (cell2mat(trial_behavior_L_50{:, 'GoCueDelay'}) == 150);
logic_220 = (cell2mat(trial_behavior_L_50{:, 'GoCueDelay'}) == 220);
logic_300 = (cell2mat(trial_behavior_L_50{:, 'GoCueDelay'}) == 300);

plot( find( all([logic_failed, logic_80], 2) ), trial_behavior_L_50{ all([logic_failed, logic_80], 2), plot_variable }, 'or',...
      find( all([~logic_failed, logic_80], 2) ), trial_behavior_L_50{ all([~logic_failed, logic_80], 2), plot_variable }, 'xr',...
      find( all([logic_failed, logic_150], 2) ), trial_behavior_L_50{ all([logic_failed, logic_150], 2), plot_variable }, 'ob',...
      find( all([~logic_failed, logic_150], 2) ), trial_behavior_L_50{ all([~logic_failed, logic_150], 2), plot_variable }, 'xb',...
      find( all([logic_failed, logic_220], 2) ), trial_behavior_L_50{ all([logic_failed, logic_220], 2), plot_variable }, 'ok',...
      find( all([~logic_failed, logic_220], 2) ), trial_behavior_L_50{ all([~logic_failed, logic_220], 2), plot_variable }, 'xk',...
      find( all([logic_failed, logic_300], 2) ), trial_behavior_L_50{ all([logic_failed, logic_300], 2), plot_variable }, 'om',...
      find( all([~logic_failed, logic_300], 2) ), trial_behavior_L_50{ all([~logic_failed, logic_300], 2), plot_variable }, 'xm' );
ylabel(str_variable); xlabel('Left Target, 50%');
ylim([0, 0.8]);

subplot(3, 3, 7)
logic_failed = logical(cell2mat(trial_behavior_L_75{:, 'TrialFailed'}));
logic_80 = (cell2mat(trial_behavior_L_75{:, 'GoCueDelay'}) == 80);
logic_150 = (cell2mat(trial_behavior_L_75{:, 'GoCueDelay'}) == 150);
logic_220 = (cell2mat(trial_behavior_L_75{:, 'GoCueDelay'}) == 220);
logic_300 = (cell2mat(trial_behavior_L_75{:, 'GoCueDelay'}) == 300);

plot( find( all([logic_failed, logic_80], 2) ), trial_behavior_L_75{ all([logic_failed, logic_80], 2), plot_variable }, 'or',...
      find( all([~logic_failed, logic_80], 2) ), trial_behavior_L_75{ all([~logic_failed, logic_80], 2), plot_variable }, 'xr',...
      find( all([logic_failed, logic_150], 2) ), trial_behavior_L_75{ all([logic_failed, logic_150], 2), plot_variable }, 'ob',...
      find( all([~logic_failed, logic_150], 2) ), trial_behavior_L_75{ all([~logic_failed, logic_150], 2), plot_variable }, 'xb',...
      find( all([logic_failed, logic_220], 2) ), trial_behavior_L_75{ all([logic_failed, logic_220], 2), plot_variable }, 'ok',...
      find( all([~logic_failed, logic_220], 2) ), trial_behavior_L_75{ all([~logic_failed, logic_220], 2), plot_variable }, 'xk',...
      find( all([logic_failed, logic_300], 2) ), trial_behavior_L_75{ all([logic_failed, logic_300], 2), plot_variable }, 'om',...
      find( all([~logic_failed, logic_300], 2) ), trial_behavior_L_75{ all([~logic_failed, logic_300], 2), plot_variable }, 'xm' );
ylabel(str_variable); xlabel('Left Target, 70%');
ylim([0, 0.8]);

subplot(3, 3, 3)
logic_failed = logical(cell2mat(trial_behavior_R_25{:, 'TrialFailed'}));
logic_80 = (cell2mat(trial_behavior_R_25{:, 'GoCueDelay'}) == 80);
logic_150 = (cell2mat(trial_behavior_R_25{:, 'GoCueDelay'}) == 150);
logic_220 = (cell2mat(trial_behavior_R_25{:, 'GoCueDelay'}) == 220);
logic_300 = (cell2mat(trial_behavior_R_25{:, 'GoCueDelay'}) == 300);

plot( find( all([logic_failed, logic_80], 2) ), trial_behavior_R_25{ all([logic_failed, logic_80], 2), plot_variable }, 'or',...
      find( all([~logic_failed, logic_80], 2) ), trial_behavior_R_25{ all([~logic_failed, logic_80], 2), plot_variable }, 'xr',...
      find( all([logic_failed, logic_150], 2) ), trial_behavior_R_25{ all([logic_failed, logic_150], 2), plot_variable }, 'ob',...
      find( all([~logic_failed, logic_150], 2) ), trial_behavior_R_25{ all([~logic_failed, logic_150], 2), plot_variable }, 'xb',...
      find( all([logic_failed, logic_220], 2) ), trial_behavior_R_25{ all([logic_failed, logic_220], 2), plot_variable }, 'ok',...
      find( all([~logic_failed, logic_220], 2) ), trial_behavior_R_25{ all([~logic_failed, logic_220], 2), plot_variable }, 'xk',...
      find( all([logic_failed, logic_300], 2) ), trial_behavior_R_25{ all([logic_failed, logic_300], 2), plot_variable }, 'om',...
      find( all([~logic_failed, logic_300], 2) ), trial_behavior_R_25{ all([~logic_failed, logic_300], 2), plot_variable }, 'xm' );
ylabel(str_variable); xlabel('Right Target, 25%');
ylim([0, 0.8]);
subplot(3, 3, 6)
logic_failed = logical(cell2mat(trial_behavior_R_50{:, 'TrialFailed'}));
logic_80 = (cell2mat(trial_behavior_R_50{:, 'GoCueDelay'}) == 80);
logic_150 = (cell2mat(trial_behavior_R_50{:, 'GoCueDelay'}) == 150);
logic_220 = (cell2mat(trial_behavior_R_50{:, 'GoCueDelay'}) == 220);
logic_300 = (cell2mat(trial_behavior_R_50{:, 'GoCueDelay'}) == 300);

plot( find( all([logic_failed, logic_80], 2) ), trial_behavior_R_50{ all([logic_failed, logic_80], 2), plot_variable }, 'or',...
      find( all([~logic_failed, logic_80], 2) ), trial_behavior_R_50{ all([~logic_failed, logic_80], 2), plot_variable }, 'xr',...
      find( all([logic_failed, logic_150], 2) ), trial_behavior_R_50{ all([logic_failed, logic_150], 2), plot_variable }, 'ob',...
      find( all([~logic_failed, logic_150], 2) ), trial_behavior_R_50{ all([~logic_failed, logic_150], 2), plot_variable }, 'xb',...
      find( all([logic_failed, logic_220], 2) ), trial_behavior_R_50{ all([logic_failed, logic_220], 2), plot_variable }, 'ok',...
      find( all([~logic_failed, logic_220], 2) ), trial_behavior_R_50{ all([~logic_failed, logic_220], 2), plot_variable }, 'xk',...
      find( all([logic_failed, logic_300], 2) ), trial_behavior_R_50{ all([logic_failed, logic_300], 2), plot_variable }, 'om',...
      find( all([~logic_failed, logic_300], 2) ), trial_behavior_R_50{ all([~logic_failed, logic_300], 2), plot_variable }, 'xm' );
ylabel(str_variable); xlabel('Right Target, 50%');
ylim([0, 0.8]);

subplot(3, 3, 9)
logic_failed = logical(cell2mat(trial_behavior_R_75{:, 'TrialFailed'}));
logic_80 = (cell2mat(trial_behavior_R_75{:, 'GoCueDelay'}) == 80);
logic_150 = (cell2mat(trial_behavior_R_75{:, 'GoCueDelay'}) == 150);
logic_220 = (cell2mat(trial_behavior_R_75{:, 'GoCueDelay'}) == 220);
logic_300 = (cell2mat(trial_behavior_R_75{:, 'GoCueDelay'}) == 300);

plot( find( all([logic_failed, logic_80], 2) ), trial_behavior_R_75{ all([logic_failed, logic_80], 2), plot_variable }, 'or',...
      find( all([~logic_failed, logic_80], 2) ), trial_behavior_R_75{ all([~logic_failed, logic_80], 2), plot_variable }, 'xr',...
      find( all([logic_failed, logic_150], 2) ), trial_behavior_R_75{ all([logic_failed, logic_150], 2), plot_variable }, 'ob',...
      find( all([~logic_failed, logic_150], 2) ), trial_behavior_R_75{ all([~logic_failed, logic_150], 2), plot_variable }, 'xb',...
      find( all([logic_failed, logic_220], 2) ), trial_behavior_R_75{ all([logic_failed, logic_220], 2), plot_variable }, 'ok',...
      find( all([~logic_failed, logic_220], 2) ), trial_behavior_R_75{ all([~logic_failed, logic_220], 2), plot_variable }, 'xk',...
      find( all([logic_failed, logic_300], 2) ), trial_behavior_R_75{ all([logic_failed, logic_300], 2), plot_variable }, 'om',...
      find( all([~logic_failed, logic_300], 2) ), trial_behavior_R_75{ all([~logic_failed, logic_300], 2), plot_variable }, 'xm' );
ylabel(str_variable); xlabel('Right Target, 75%');
ylim([0, 0.8]);


%{
subplot(3, 1, 2)
logic_failed = logical(cell2mat(trial_behavior_R{:, 'Triallogic_failed'}));
plot(find(logic_failed), trial_behavior_R{logic_failed, 'RT'}, 'or', find(~logic_failed), trial_behavior_R{~logic_failed, 'RT'}, 'xb');
ylabel('RT (s)');
xlabel('Right Target')

subplot(3, 1, 3)
logic_failed = logical(cell2mat(trial_behavior_L{:, 'Triallogic_failed'}));
plot(find(logic_failed), trial_behavior_R{logic_failed, 'RT'}, 'or', find(~logic_failed), trial_behavior_R{~logic_failed, 'RT'}, 'xb');
ylabel('RT (s)');
xlabel('Left Target')

%%
figure(2)
subplot(3, 1, 1)
logic_failed = logical(cell2mat(trial_behavior_M{:, 'Triallogic_failed'}));
plot(find(logic_failed), trial_behavior_M{logic_failed, 'MT'}, 'or', find(~logic_failed), trial_behavior_M{~logic_failed, 'MT'}, 'xb');
ylabel('MT(s)');
xlabel('Middle Target')

subplot(3, 1, 2)
logic_failed = logical(cell2mat(trial_behavior_R{:, 'Triallogic_failed'}));
plot(find(logic_failed), trial_behavior_R{logic_failed, 'MT'}, 'or', find(~logic_failed), trial_behavior_R{~logic_failed, 'MT'}, 'xb');
ylabel('MT (s)');
xlabel('Right Target')

subplot(3, 1, 3)
logic_failed = logical(cell2mat(trial_behavior_L{:, 'Triallogic_failed'}));
plot(find(logic_failed), trial_behavior_R{logic_failed, 'MT'}, 'or', find(~logic_failed), trial_behavior_R{~logic_failed, 'MT'}, 'xb');
ylabel('MT (s)');
xlabel('Left Target')
%}
%%
trial_no = sortrows(trial_no, 1);