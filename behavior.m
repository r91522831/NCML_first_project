function [ behave ] = behavior( data )
%behavior Summary of this function goes here
%   Detailed explanation goes here

%% Start process
behave = data(:, 2:9);
behave.RT = zeros(height(behave), 1);
behave.MT = zeros(height(behave), 1);
behave.pVel = zeros(height(behave), 1);
behave.t2pVel = zeros(height(behave), 1);

for i = 1:height(data)
    tmp_event_time = data.EventTime{i};
    % find the event time for movement initiation
    ind_move_ini = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'MOVE_OUT_SP')));
    if isempty(ind_move_ini)
        fprintf('The subject did not move out of the starting area in trial %d.\n', i);
        clear ind_move_ini
        continue;
    end
    time_move_ini = tmp_event_time{ind_move_ini, 2};
    
    ind_cue_on = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'CUE_ON')));
    ind_location_on = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'LOCATION_ON')));
    ind_go_on = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'TARGET_GO')));
    time_cue_on = tmp_event_time{ind_cue_on, 2};
    time_location_on = tmp_event_time{ind_location_on, 2};
    time_go_on = tmp_event_time{ind_go_on, 2};
    % compute reaction time in seconds
    behave{i, 'RT'} = time_move_ini - time_go_on;
    
    ind_reach = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'REACH_TARGET')));
    if isempty(ind_reach)
        fprintf('The subject did not reach the target area in trial %d.\n', i);
        clear ind_move_ini
        continue;
    end
    time_reach = tmp_event_time{ind_reach, 2};
    % compute movement time in seconds
    behave{i, 'MT'} = time_reach - time_move_ini;
    
    % compute peak velocity
    tmp_kinematic = data.KinARM{i};
    dt = 0.001; % 1 ms, sampling rate is 1K Hz
    cutoff = 10; % 10 Hz
    tmp_kine_filter = filtmat_class( dt, cutoff, tmp_kinematic{:, :});
    [~, tmp_grad] = gradient(tmp_kine_filter);
    tmp_derivative = tmp_grad ./ dt;
    tmp_vel = sqrt(tmp_derivative(:, 1).^2 + tmp_derivative(:, 2).^2);
    [peak_vel, ind_peak] = max(tmp_vel);
    
    behave{i, 'pVel'} = peak_vel;
    behave{i, 't2pVel'} = ind_peak * dt - time_move_ini;
end
behave = sortrows(behave, 'order');

end

