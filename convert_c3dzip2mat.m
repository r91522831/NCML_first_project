function [ data ] = convert_c3dzip2mat( filepath, filename )
% convert_c3dzip2mat Summary of this function goes here
%   Detailed explanation goes here
full_filename = fullfile(filepath, filename);

data = zip_load(full_filename);             % Loads all c3d_files into a new structure called 'data'.
data.c3d = KINARM_add_hand_kinematics(data.c3d);	% Add hand velocity, acceleration and commanded forces to the data structure

end

