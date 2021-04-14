% path to main folder
dir_current = fileparts(which('current_directory.m'));

% path to data
dir_data = fullfile(dir_current,'\data');
if exist(dir_data,'dir') ~= 7
    dir_data=0;
end
% path to figures
dir_fig = fullfile(dir_current,'\figures');
if exist(dir_fig,'dir') ~= 7
    disp('create folder \figure')
    mkdir(dir_fig);
end
% path to tables
dir_tab = fullfile(dir_current,'\tables');
if exist(dir_tab,'dir') ~= 7
    disp('create folder \figure')
    mkdir(dir_tab);
end

