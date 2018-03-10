close all; clear; clc;

%%
% User selection of the main folder
SUB_dir = uigetdir;
SUB_list = dir( fullfile(SUB_dir, '*.zip') );

for j = 1:length(SUB_list)
    % Load c3d .zip file and convert data into a .mat file
    data = convert_c3dzip2mat( SUB_dir, SUB_list(j).name );
    
    [~, subj_dir, ~] = fileparts(SUB_dir);
    [~, filename, ~] = fileparts(SUB_list(j).name);
    
    save_path = fullfile( fileparts(fileparts(SUB_dir)), 'mat files', subj_dir );
    if ~exist(save_path, 'dir'); mkdir(save_path); end
    save_filename = fullfile(save_path, filename);
    
    save(save_filename, 'data')
    
end
