close all; clear; clc;

%%
% User selection of the main folder
SUB_dir = uigetdir;
tmp_list = dir( fullfile(SUB_dir) );

SUB_list = tmp_list(arrayfun(@(x) x.name(1), tmp_list) ~= '.'); % get rid of '.' and '..' folders
clear tmp_list

for i = 1:length(SUB_list)
    tmp_path = fullfile(SUB_dir, SUB_list(i).name);
    tmp_list = dir( fullfile(tmp_path, '*.zip') );
    
    for j = 1:length(tmp_list)
        % Load c3d .zip file and convert data into a .mat file
        data = convert_c3dzip2mat( tmp_path, tmp_list(j).name );
        
        [~, subj_dir, ~] = fileparts(tmp_path);
        [~, filename, ~] = fileparts(tmp_list(j).name);
        
        save_path = fullfile(fileparts(SUB_dir), 'mat files', subj_dir);
        if ~exist(save_path, 'dir'); mkdir(save_path); end
        save_filename = fullfile(save_path, filename);
        
        save(save_filename, 'data')

    end
end