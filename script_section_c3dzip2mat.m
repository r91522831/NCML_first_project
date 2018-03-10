 close all; clear; clc;

%%
[filename, filepath, ~] = uigetfile('*.zip');

data = convert_c3dzip2mat( filepath, filename );

full_filename = fullfile(filepath, filename);
[upper_dir, subj_dir, ~] = fileparts(fileparts(filepath));
[~, filename, ~] = fileparts(full_filename);
save_path = fullfile(fileparts(upper_dir), 'mat files', subj_dir);
if ~exist(save_path, 'dir'); mkdir(save_path); end
save_filename = fullfile(save_path, filename);

save(save_filename, 'data')