% --------------------------------------------------------------------------------------------------------------------------------- %
% FIGURE C3 - Monte Carlo Simulations: Impact of Economy Size on the Equilibrium Distribution of Cutoffs (Constrained/Truncated DA) %
% --------------------------------------------------------------------------------------------------------------------------------- %

clear

if exist('set_path.m', 'file') == 2
    set_path
else 
    error('ERROR: specify path to folder containing replication files in change_directory.m and execute script')  
end

% Parameters
M = 500; % 500 MC samples
J = 6; % 6 schools
A = 4; % Maximum size of ROL = 4
unit_cost = 0; % No cost of ranking more schools

%%% Figure C3(a): 100 students
propor = 1; % 100 students
MC_data_simulation; % simulate data
writetable(FIG.DA_CONS.Cutoffs,fullfile(dir_fig,'xfigure_C3a.xlsx'),'Sheet',1,'Range','A1');
clearvars -except dir_fig M J A unit_cost;

%%% Figure C3(b): 200 students
propor = 2; % 200 students
MC_data_simulation; % simulate data
writetable(FIG.DA_CONS.Cutoffs,fullfile(dir_fig,'xfigure_C3b.xlsx'),'Sheet',1,'Range','A1');
clearvars -except dir_fig M J A unit_cost;

%%% Figure C3(c): 500 students
propor = 5; % 500 students
MC_data_simulation; % simulate data
writetable(FIG.DA_CONS.Cutoffs,fullfile(dir_fig,'xfigure_C3c.xlsx'),'Sheet',1,'Range','A1');
clearvars -except dir_fig M J A unit_cost;

%%% Figure C3(d): 5,000 students
propor = 50; % 5,000 students
MC_data_simulation; % simulate data
writetable(FIG.DA_CONS.Cutoffs,fullfile(dir_fig,'xfigure_C3d.xlsx'),'Sheet',1,'Range','A1');
