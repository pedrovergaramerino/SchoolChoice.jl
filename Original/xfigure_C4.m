% ----------------------------------------------------------------------------------------------------------------------------------------- %
% FIGURE C4 - Monte Carlo Simulations: Impact of the Marginal Cost of Applying to Schools on Equilibrium Outcomes (500 Students, 6 Schools) %
% ----------------------------------------------------------------------------------------------------------------------------------------- %

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
A = 6; % Maximum size of ROL = 4

comparative_statics = NaN(6,4);

%%% c = 0
unit_cost = 0; % No cost of ranking more schools
MC_data_simulation; % simulate data
comparative_statics(1,:) = [PARAM.unit_cost, STATS.size_ROL, STATS.pct_assigned_pref_feas, STATS.pct_topk];
clearvars -except comparative_statics M propor J A dir_fig;

%%% c = 1e-6
unit_cost = 1e-6; % No cost of ranking more schools
MC_data_simulation; % simulate data
comparative_statics(2,:) = [PARAM.unit_cost, STATS.size_ROL, STATS.pct_assigned_pref_feas, STATS.pct_topk];
clearvars -except comparative_statics M propor J A dir_fig;

%%% c = 1e-3
unit_cost = 1e-3; % No cost of ranking more schools
MC_data_simulation; % simulate data
comparative_statics(3,:) = [PARAM.unit_cost, STATS.size_ROL, STATS.pct_assigned_pref_feas, STATS.pct_topk];
clearvars -except comparative_statics M propor J A dir_fig;

%%% c = 1e-2
unit_cost = 1e-2; % No cost of ranking more schools
MC_data_simulation; % simulate data
comparative_statics(4,:) = [PARAM.unit_cost, STATS.size_ROL, STATS.pct_assigned_pref_feas, STATS.pct_topk];
clearvars -except comparative_statics M propor J A dir_fig;

%%% c = 0.1
unit_cost = 0.1; % No cost of ranking more schools
MC_data_simulation; % simulate data
comparative_statics(5,:) = [PARAM.unit_cost, STATS.size_ROL, STATS.pct_assigned_pref_feas, STATS.pct_topk];
clearvars -except comparative_statics M propor J A dir_fig;

%%% c = 1
unit_cost = 1; % No cost of ranking more schools
MC_data_simulation; % simulate data
comparative_statics(6,:) = [PARAM.unit_cost, STATS.size_ROL, STATS.pct_assigned_pref_feas, STATS.pct_topk];
clearvars -except comparative_statics dir_fig;

%%% Export data

fig_C4 = table((1:6)',comparative_statics(:,1),comparative_statics(:,2),comparative_statics(:,3),comparative_statics(:,4), ...
                'VariableNames',{'id', 'application_cost' 'm_ROL_length' 'pct_stable' 'pct_WTT'});

writetable(fig_C4,fullfile(dir_fig,'xfigure_C4.xlsx'),'Sheet',1,'Range','A1');
