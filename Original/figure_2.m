% ---------------------------------------------------------------------------------------------------------------------------------- %
% FIGURE 2 - Monte Carlo Simulations: Impact of Economy Size on the Equilibrium Distribution of Cutoffs (Constrained / Truncated DA) %
% ---------------------------------------------------------------------------------------------------------------------------------- %

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

%%% Figure 1(a): 100 students, 6 schools

propor = 1; % 500 students
MC_data_simulation; % simulate data
writetable(FIG.DA_CONS.Cutoffs,fullfile(dir_fig,'figure_2a.xlsx'),'Sheet',1,'Range','A1');
clearvars -except M J A dir_fig unit_cost;

%%% Figure 1(b): 500 students, 6 schools

propor = 5; % 500 students
MC_data_simulation; % simulate data
writetable(FIG.DA_CONS.Cutoffs,fullfile(dir_fig,'figure_2b.xlsx'),'Sheet',1,'Range','A1');
