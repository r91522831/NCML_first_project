function [ behave ] = behavior( data )
%behavior Summary of this function goes here
%   Detailed explanation goes here

%% Start process
behave = data(:, 2:9);
behave.RT = zeros(height(behave), 1);
behave.MT = zeros(height(behave), 1);
behave.pVel = zeros(height(behave), 1);
behave.t2pVel = zeros(height(behave), 1);
behave.Perturbed = zeros(height(behave), 1);
behave.manMovementOnset = zeros(height(behave), 1);
behave.manRT = zeros(height(behave), 1);
behave.manMT = zeros(height(behave), 1);
behave.man_t2pVel = zeros(height(behave), 1);
behave.iniAngDevi = zeros(height(behave), 1); % initial angle deviation from ideal straight line. (CCW is positive)


for i = 1:height(data)
    tmp_event_time = data.EventTime{i};
    % find the event time for movement initiation by the MOVE_OUT_SP event
    ind_move_ini = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'MOVE_OUT_SP')));
    if isempty(ind_move_ini)
        fprintf('The subject did not move out of the starting area in trial %d.\n', i);
        clear ind_move_ini
        continue;
    end
    time_move_ini = tmp_event_time{ind_move_ini, 2};
    
    ind_cue_on = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'CUE_ON')), 1);
    time_cue_on = tmp_event_time{ind_cue_on, 2};
%     ind_location_on = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'LOCATION_ON')), 1);
    
%     time_location_on = tmp_event_time{ind_location_on, 2};

    ind_go_on = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'TARGET_GO')), 1);
    time_go_on = tmp_event_time{ind_go_on, 2};
    
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
    % compute the manual movement onset by 1% pVel backward from pVel
    for j = ind_peak:-1:1
        if tmp_vel(j) < 0.01 * tmp_vel(ind_peak)
            ind_mov = j;
            break;
        end
    end
    time_movement_onset = ind_mov * dt;
    behave{i, 'manMovementOnset'} = time_movement_onset;
    
    % compute time in seconds based on logged KinARM event
    behave{i, 'RT'} = time_move_ini - time_go_on;
    behave{i, 't2pVel'} = ind_peak * dt - time_move_ini;
    % compute time in seconds based on velocity profile
    behave{i, 'manRT'} = time_movement_onset - time_go_on;
    behave{i, 'man_t2pVel'} = ind_peak * dt - time_movement_onset;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind_pert_on = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'PERTURB_ON')), 1);
    if isempty(ind_pert_on)
        behave{i, 'Perturbed'} = 0;
    else
        behave{i, 'Perturbed'} = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ind_reach = find(~cellfun(@isempty, strfind(tmp_event_time(:, 1), 'REACH_TARGET')), 1);
    if isempty(ind_reach)
        fprintf('The subject did not reach the target area in trial %d.\n', i);
        clear ind_reach
        continue;
    end
    time_reach = tmp_event_time{ind_reach, 2};
    % compute movement time in seconds based on logged KinARM event
    behave{i, 'MT'} = time_reach - time_move_ini;
    % compute movement time in seconds based on velocity profile
    behave{i, 'manMT'} = time_reach - time_movement_onset;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % compute initial movement angle deviate from ideal straight initial-end targets line
    ind_cue = floor(time_cue_on / dt);
    ini_position = mean(tmp_kine_filter(ind_cue:ind_mov, 1:2), 1);
    pVel_postion = tmp_kine_filter(ind_peak, 1:2);
    fin_position = mean(tmp_kine_filter((end-1000):end, 1:2), 1);
    ideal_vector = fin_position - ini_position;
    movement_ini_vector = pVel_postion - ini_position;
    
    tmp_sgn = sign(cross([ideal_vector, 0], [movement_ini_vector, 0]));
    behave{i, 'iniAngDevi'} = tmp_sgn(1, 3) * acosd(dot(ideal_vector, movement_ini_vector) / (norm(ideal_vector) * norm(movement_ini_vector)));
    
end
behave = sortrows(behave, 'order');

end

