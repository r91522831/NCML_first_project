function [ behave ] = behavior( data, subid )
%behavior Summary of this function goes here
%   Detailed explanation goes here

%% Start process
dt = 0.001; % 1 ms, sampling rate is 1K Hz
cutoff = 10; % 10 Hz

behave = data(:, 2:9);
behave.events = cell(height(behave), 1);
behave.pos = cell(height(behave), 1);
behave.vel = cell(height(behave), 1);
behave.acc = cell(height(behave), 1);
behave.force = cell(height(behave), 1);
behave.CueRT = NaN(height(behave), 1);
behave.RT = NaN(height(behave), 1);
behave.MT = NaN(height(behave), 1);
behave.pVel = NaN(height(behave), 1);
behave.t2pVel = NaN(height(behave), 1);
behave.pAcc = NaN(height(behave), 1);
behave.t2pAcc = NaN(height(behave), 1);
behave.pForce = NaN(height(behave), 1);
behave.t2pForce = NaN(height(behave), 1);
behave.Perturbed = false(height(behave), 1);
behave.manMoveOnset = NaN(height(behave), 1);
behave.manMoveEnd = NaN(height(behave), 1);
behave.manCueRT = NaN(height(behave), 1);
behave.manRT = NaN(height(behave), 1);
behave.manMT = NaN(height(behave), 1);
behave.man_t2pVel = NaN(height(behave), 1);
behave.man_t2pAcc = NaN(height(behave), 1);
behave.man_t2pForce = NaN(height(behave), 1);
behave.maxPathDevi = NaN(height(behave), 1); % maximun trajectory deviation from ideal straight line.
behave.pathLength = NaN(height(behave), 1); % trajectory length.
behave.iniAngDevi = NaN(height(behave), 1); % initial angle deviation from ideal straight line. (CCW is positive)


for i = 1:height(data)
    tmp_event_time = data.EventTime{i};
    % find the event time for movement initiation by the MOVE_OUT_SP event
    tmp_ind = contains(tmp_event_time(:, 1), 'MOVE_OUT_SP');
    if ~any(tmp_ind)
        fprintf('Subject %d did not move out of the starting area in trial %d.\n', subid, i);
        clear ind_move_ini
        continue;
    end
    time_move_out = tmp_event_time{tmp_ind, 2};
    ind_move_out = floor(time_move_out / dt);
    
    time_cue_on = tmp_event_time{ contains(tmp_event_time(:, 1), 'CUE_ON'), 2 };
    ind_cue_on = floor(time_cue_on / dt);
%     time_location_on = tmp_event_time{ contains(tmp_event_time(:, 1), 'LOCATION_ON'), 2 };
%     ind_location_on = floor(time_location_on / dt);

    time_go_on = tmp_event_time{ contains(tmp_event_time(:, 1), 'TARGET_GO'), 2 };
    ind_go_on = floor(time_go_on / dt);
    
    % Get all endpoint kinematics, X, Y, XVel, YVel, XAcc, YAcc
    tmp_kinematics = data.KinARM{i}{:, {'X', 'Y', 'XVel', 'YVel', 'XAcc', 'YAcc'}};
    tmp_kine_filter = filtmat_class( dt, cutoff, tmp_kinematics);
    
    tmp_posXY = tmp_kine_filter(:, 1:2);
    tmp_velXY = tmp_kine_filter(:, 3:4);
%     tmp_accXY = tmp_kine_filter(:, 5:6);
    % compute pos
    tmp_pos = sqrt(sum(tmp_posXY.^2, 2));
    behave{i, 'pos'} = {[tmp_pos, tmp_posXY]};
    % compute velocity
    tmp_vel = sqrt(sum(tmp_velXY.^2, 2));
    behave{i, 'vel'} = {tmp_vel};
    % find peak velocity
    [peak_vel, ind_pVel] = max(tmp_vel);
    behave{i, 'pVel'} = peak_vel;
    time_pVel = ind_pVel * dt;
    % compute acceleration
    tmp_acc = gradient(tmp_vel);
    behave{i, 'acc'} = {tmp_acc};
    [peak_acc, ind_pAcc] = max(tmp_acc);
    behave{i, 'pAcc'} = peak_acc;
    time_pAcc = ind_pAcc * dt;
    
    % compute the manual movement onset by 5% pVel backward from MOVE_OUT_SP
    for j = ind_move_out:-1:1
        if tmp_vel(j) < 0.05 * tmp_vel(ind_pVel)
            ind_move_onset = j;
            break;
        end
    end
    time_move_onset = ind_move_onset * dt;
    behave{i, 'manMoveOnset'} = time_move_onset;
    
    % compute time in seconds based on logged KinARM event
    behave{i, 'RT'} = time_move_out - time_go_on;
    behave{i, 't2pVel'} = time_pVel - time_move_out;
    behave{i, 't2pAcc'} = time_pAcc - time_move_out;
    behave{i, 'CueRT'} = time_move_out - time_cue_on;
    % compute time in seconds based on velocity profile
    behave{i, 'manRT'} = time_move_onset - time_go_on;
    behave{i, 'man_t2pVel'} = time_pVel - time_move_onset;
    behave{i, 'man_t2pAcc'} = time_pAcc - time_move_onset;
    
    % comupte time in seconds from cue on to move onset
    behave{i, 'manCueRT'} = time_move_onset - time_cue_on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp_ind = contains(tmp_event_time(:, 1), 'PERTURB_ON');
    behave{i, 'Perturbed'} = any(tmp_ind);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp_ind = contains(tmp_event_time(:, 1), 'REACH_TARGET');
    if ~any(tmp_ind)
        fprintf('Subject %d did not reach the target area in trial %d.\n', subid, i);
        clear ind_reach
        continue;
    elseif sum(tmp_ind) > 1
        fprintf('Subject %d reached the target area more than once in trial %d.\n', subid, i);
    end
    time_reach = tmp_event_time{tmp_ind, 2};
    ind_reach = floor(time_reach / dt);
    
    % compute the manual movement end by 2% pVel forward from REACH_TARGET
    for j = ind_reach:1:length(tmp_vel)
        if tmp_vel(j) < 0.02 * tmp_vel(ind_pVel)
            ind_move_end = j;
            break;
        end
    end
    time_move_end = ind_move_end * dt;
    behave{i, 'manMoveEnd'} = time_move_end;
    
    % compute movement time in seconds based on logged KinARM event
    behave{i, 'MT'} = time_reach - time_move_out;
    % compute movement time in seconds based on velocity profile
    behave{i, 'manMT'} = time_move_end - time_move_onset;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % pertubed trial has endpoint force information
    if behave{i, 'Perturbed'}
        time_pertb = tmp_event_time{contains(tmp_event_time(:, 1), 'PERTURB_ON'), 2};
        ind_pertb = floor(time_pertb / dt); % next to the ind_move_out
        % comupte the hand force
        tmp_kinetic = data.KinARM{i}{:, {'XForceCMD', 'YForceCMD'}};
        tmp_force = sqrt(sum(tmp_kinetic.^2, 2));
        behave{i, 'force'} = {tmp_force};
        
        [peak_force, ind_pForce] = max(tmp_force);
        behave{i, 'pForce'} = peak_force;
        time_pForce = ind_pForce * dt;
        
        behave{i, 't2pForce'} = time_pAcc - time_move_out;
        behave{i, 'man_t2pForce'} = time_pForce - time_move_onset;
    end
           
    % compute movement path deviate from ideal straight initial-end targets line
    ini_position = tmp_posXY(ind_move_onset, :);
    fin_position = tmp_posXY(ind_move_end, :);
    ideal_vector = fin_position - ini_position;
    dist = zeros(length(ind_move_onset:ind_move_end), 1);
    for pt = ind_move_onset:ind_move_end
        tmp_vector = fin_position - tmp_posXY(pt, :);
        dist(pt, 1) = norm(cross([ideal_vector, 0], [tmp_vector, 0]))/norm(ideal_vector);
    end    
    behave{i, 'maxPathDevi'} = max(dist);
    
    % compute movement path length
    behave{i, 'pathLength'} = sum(tmp_vel(ind_move_onset:ind_move_end, :) .* dt);
    
    % compute initial movement angle deviate from ideal straight initial-end targets line
    % Since the angle comuptation is very sensitive, the final position is
    % defined at the end of the reach instead of the reach event
    ind_cue = floor(time_cue_on / dt);
    ini_position = mean(tmp_posXY(ind_cue:ind_move_onset, :), 1);
    fin_position = mean(tmp_posXY((end-1000):end, :), 1);
    ideal_vector = fin_position - ini_position;
    pVel_postion = tmp_posXY(ind_pVel, :);
    movement_ini_vector = pVel_postion - ini_position;
    
    tmp_sgn = sign(cross([ideal_vector, 0], [movement_ini_vector, 0]));
    behave{i, 'iniAngDevi'} = tmp_sgn(1, 3) * acosd(dot(ideal_vector, movement_ini_vector) / (norm(ideal_vector) * norm(movement_ini_vector)));
 
    
    behave{i, 'events'} = {[ind_cue_on, ind_go_on, ind_move_onset, ind_move_out, ind_reach, ind_move_end]};
end
end

