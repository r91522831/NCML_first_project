function [trials, load_table, target_table, tp_table] = convert_data2tables(data, side)
if nargin < 2
    side = 'Right'; % default right handed
end

% convert_data2tables
%   Detailed explanation goes here

%
% get the load table, target table, and trial protocol table
% they are all the same for all trials
load_table = [ data.c3d(1).LOAD_TABLE.A, data.c3d(1).LOAD_TABLE.B,...
               data.c3d(1).LOAD_TABLE.C, data.c3d(1).LOAD_TABLE.D ];
load_table( ~any(load_table, 2), : ) = []; % remove rows with all zeros
load_table = array2table(load_table);
load_table.Properties.VariableNames = {'A', 'B', 'C', 'D'};

target_table = [ data.c3d(1).TARGET_TABLE.X, data.c3d(1).TARGET_TABLE.Y,...
                 data.c3d(1).TARGET_TABLE.Major, data.c3d(1).TARGET_TABLE.Minor, data.c3d(1).TARGET_TABLE.Rotation];
target_table( ~any(target_table, 2), : ) = []; % remove rows with all zeros
target_table = array2table(target_table);
target_table.Properties.VariableNames = {'X', 'Y', 'Major', 'Minor','Rotation'};

% trial protocol table:
% { probability no shrink target, probability perturbation,...
%   time to show target (fix), time to show target (random),...
%   max time for reach, time to stay in target,...
%   time after each trial, perturbation cue dispaly time,...
%   probability perturb in CCW, contour target display time,...
%   max time allowed to be out of target, max time to initiate movement,...
%   reach time lower bound, upper bound,...
%   is movement time constrainted? }
tp_table = [ data.c3d(1).TP_TABLE.Larger_Target_Prob, data.c3d(1).TP_TABLE.Perturb_Prob,...
             data.c3d(1).TP_TABLE.Ready_Delay__Fix_, data.c3d(1).TP_TABLE.Ready_Delay__Rand_,...
             data.c3d(1).TP_TABLE.Max_Reach_Time, data.c3d(1).TP_TABLE.Target_Delay,...
             data.c3d(1).TP_TABLE.Post_Trial_Delay, data.c3d(1).TP_TABLE.Cue_Delay,...
             data.c3d(1).TP_TABLE.Direction_Prob, data.c3d(1).TP_TABLE.Go_Delay,...
             data.c3d(1).TP_TABLE.Max_Recover_Time, data.c3d(1).TP_TABLE.Max_React_Time,...
             data.c3d(1).TP_TABLE.Lower_Reach_Time, data.c3d(1).TP_TABLE.Upper_Reach_Time,...
             data.c3d(1).TP_TABLE.Constrained_MT ];
tp_table( ~any(tp_table, 2), : ) = []; % remove rows with all zeros
tp_table = array2table(tp_table);
tp_table.Properties.VariableNames = { 'LTargetProb', 'PerturbProb', 'ReadyDelayFix', 'ReadyDelayRandom',...
                                      'MaxReachTime', 'TargetFinishDelay', 'NextTrialDelay', 'LocationDelay',...
                                      'DirectionProb', 'GoCueDelay', 'MaxRecoverTime', 'MaxReactTime',...
                                      'MTLowB', 'MTUpB', 'IsMTConstrained' };

%%
% data format:
% { system time string, trial protocol ID, trial protocol repetition,...
%   trial sequence, target location, perturbation probability,...
%   go delay, is trial failed?, score, events timing, KinARM data,...
%    }
trials = cell(length(data.c3d), 11);
for i = 1:length(data.c3d)
    tmp_tpID = data.c3d(i).TRIAL.TP;
    % system time string: Time that data logging started for this trial.
    trials{i, 1} = data.c3d(i).TRIAL.TIME;
    % trial protocol ID: Trial protocol (TP table row) for this trial.
    trials{i, 2} = tmp_tpID;
    % trial protocol repetition: Trial is the Nth instance of this TP in the set.
    trials{i, 3} = data.c3d(i).TRIAL.TP_NUM;
    % trial sequence: Trial is the Nth trial of the set.
    trials{i, 4} = data.c3d(i).TRIAL.TRIAL_NUM;
    % trial condition:
    % target location: 
    switch data.c3d(i).TP_TABLE.Cue_Location(tmp_tpID)
        case 4
            trials{i, 5} = 'M';
        case 5
            trials{i, 5} = 'R';
        case 6
            trials{i, 5} = 'L';
        otherwise
            trials{i, 5} = [];
    end
    % perturb probability
    trials{i, 6} = tp_table{tmp_tpID, 'PerturbProb'};
    % time from perturb cue to go signal
    trials{i, 7} = tp_table{tmp_tpID, 'GoCueDelay'};
    % failed trial: 1 if trial is an error trial; 0 otherwise.
    trials{i, 8} = data.c3d(i).TRIAL.IS_ERROR;
    % score: score result of the trial
    if isfield(data.c3d(i), 'score')
        trials{i, 9} = data.c3d(i).score(end, 1);
    else
        trials{i, 9} = [];
    end
    % events timing
    trials{i, 10} = [cellstr(char(data.c3d(i).EVENTS.LABELS)), num2cell(data.c3d(i).EVENTS.TIMES')];
    % KinARM data
    tmp_KinTable = [ data.c3d(i).([side, '_HandX']), data.c3d(i).([side, '_HandY']),...
                     data.c3d(i).([side, '_HandXVel']), data.c3d(i).([side, '_HandYVel']),...
                     data.c3d(i).([side, '_HandXAcc']), data.c3d(i).([side, '_HandYAcc']),...
                     data.c3d(i).([side, '_Hand_ForceCMD_X']), data.c3d(i).([side, '_Hand_ForceCMD_Y']),...
                     data.c3d(i).([side, '_L1Ang']), data.c3d(i).([side, '_L2Ang']),...
                     data.c3d(i).([side, '_L1Vel']), data.c3d(i).([side, '_L2Vel']),...
                     data.c3d(i).([side, '_L1Acc']), data.c3d(i).([side, '_L2Acc']),...
                     data.c3d(i).([side, '_M1TorCMD']), data.c3d(i).([side, '_M2TorCMD']) ];
    tmp_KinTable = array2table(tmp_KinTable);
    tmp_KinTable.Properties.VariableNames = { 'X', 'Y', 'XVel', 'YVel', 'XAcc', 'YAcc', 'XForceCMD', 'YForceCMD', ...
                                              'L1Ang', 'L2Ang', 'L1Vel', 'L2Vel', 'L1Acc', 'L2Acc', 'M1TorCMD', 'M2TorCMD' };
    trials{i, 11} = tmp_KinTable;
    
end
trials = array2table(trials);
trials.Properties.VariableNames = { 'Time', 'TPID', 'Repetition', 'order',...
                                    'TargetLocation', 'PertProb', 'GoCueDelay',...
                                    'TrialFailed', 'Score',...
                                    'EventTime', 'KinARM' };
end

