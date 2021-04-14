% --------------------------------------------------------------------------------------------------------- %
% FIGURE C2 - Monte Carlo Simulations: Equilibrium Distribution of School Cutoffs (6 schools, 500 students) %
% --------------------------------------------------------------------------------------------------------- %

clear

if exist('set_path.m', 'file') == 2
    set_path
else 
    error('ERROR: specify path to folder containing replication files in change_directory.m and execute script')  
end

% Parameters
propor = 5; % 500 students
M = 500; % 500 MC samples
J = 6; % 6 schools

%%% Figure C2(a): Constrained/truncated DA

A = 4; % Maximum size of ROL = 4
unit_cost = 0; % No cost of ranking more schools
MC_data_simulation; % simulate data
writetable(FIG.DA_CONS.Cutoffs,fullfile(dir_fig,'xfigure_C2a.xlsx'),'Sheet',1,'Range','A1');
clearvars -except dir_fig M J propor;

%%% Figure C2(b): Unconstrained DA with cost

A = 6; % Maximum size of ROL = 6
unit_cost = 1e-6; % Cost of ranking more school
MC_data_simulation; % simulate data
writetable(FIG.DA_CONS.Cutoffs,fullfile(dir_fig,'xfigure_C2b.xlsx'),'Sheet',1,'Range','A1');
