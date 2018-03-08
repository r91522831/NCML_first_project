close all; clear; clc;

%%
% User selection of the main folder
SUB_dir = uigetdir;
tmp_list = dir( fullfile(SUB_dir) );

SUB_list = tmp_list(arrayfun(@(x) x.name(1), tmp_list) ~= '.'); % get rid of '.' and '..' folders
clear tmp_list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%
[filename, filepath, ~] = uigetfile('*.zip');
full_filename = fullfile(filepath, filename);

data = zip_load(full_filename);		% Loads all c3d_files into a new structure called 'data'.
data = KINARM_add_hand_kinematics(data);					% Add hand velocity, acceleration and commanded forces to the data structure
clear filename filepath full_filename

%%
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
%   trial sequence, events timing, is trial failed?, score, KinARM data,...
%    }
trial_no = cell(length(data.c3d), 11);
for i = 1:length(data.c3d)
    tmp_tpID = data.c3d(i).TRIAL.TP;
    % system time string: Time that data logging started for this trial.
    trial_no{i, 1} = data.c3d(i).TRIAL.TIME;
    % trial protocol ID: Trial protocol (TP table row) for this trial.
    trial_no{i, 2} = tmp_tpID;
    % trial protocol repetition: Trial is the Nth instance of this TP in the set.
    trial_no{i, 3} = data.c3d(i).TRIAL.TP_NUM;
    % trial sequence: Trial is the Nth trial of the set.
    trial_no{i, 4} = data.c3d(i).TRIAL.TRIAL_NUM;
    % trial condition:
    % target location: 
    switch data.c3d(i).TP_TABLE.Cue_Location(tmp_tpID)
        case 4
            trial_no{i, 5} = 'M';
        case 5
            trial_no{i, 5} = 'R';
        case 6
            trial_no{i, 5} = 'L';
        otherwise
            trial_no{i, 5} = [];
    end
    % perturb probability
    trial_no{i, 6} = tp_table{tmp_tpID, 'PerturbProb'};
    % time from perturb cue to go signal
    trial_no{i, 7} = tp_table{tmp_tpID, 'GoCueDelay'};
    
    % events timing
    trial_no{i, 8} = [data.c3d(i).EVENTS.LABELS', num2cell(data.c3d(i).EVENTS.TIMES')];
    % failed trial: 1 if trial is an error trial; 0 otherwise.
    trial_no{i, 9} = data.c3d(i).TRIAL.IS_ERROR;
    % score: score result of the trial
    trial_no{i, 10} = data.c3d(i).score(end, 1);
    % KinARM data
    trial_no{i, 11} = [ data.c3d(i).Right_HandX, data.c3d(i).Right_HandY,...
                       data.c3d(i).Right_L1Ang, data.c3d(i).Right_L2Ang,...
                       data.c3d(i).Right_L1Vel, data.c3d(i).Right_L2Vel,...
                       data.c3d(i).Right_L1Acc, data.c3d(i).Right_L2Acc,...
                       data.c3d(i).Right_M1TorCMD, data.c3d(i).Right_M2TorCMD ];
    
end

%%
trial_no = sortrows(trial_no, 1);