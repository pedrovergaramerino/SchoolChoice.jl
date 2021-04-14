% ------------------------------------------------------- %
% TABLE C1 -  Monte Carlo Simulations: Summary Statistics %
% ------------------------------------------------------- %

clear

if exist('set_path.m', 'file') == 2
    set_path
else 
    error('ERROR: specify path to folder containing replication files in change_directory.m and execute script')  
end

%%% Column (1): Constrained/truncated DA
propor = 5; % 500 students
M = 500; % 500 MC samples
J = 6; % 6 schools
A = 4; % Maximum size of ROL = 4
unit_cost = 0; % No cost of ranking more schools
MC_data_simulation; % simulate data

col1.params = PARAM;
col1.stats = STATS;
clearvars -except col1 dir_tab

%%% Column (2): Unconstrained DA with cost

% Simulate data

propor = 5; % 500 students
M = 500; % 500 MC samples
J = 6; % 6 schools
A = 6; % Maximum size of ROL = 4
unit_cost = 1e-6; % No cost of ranking more schools
MC_data_simulation; % simulate data

col2.params = PARAM;
col2.stats = STATS;
clearvars -except col1 col2 dir_tab

%%% TABLE
tab_C1 = fullfile(dir_tab,'xtable_C1.txt');
f =  fopen(tab_C1, 'wt', 'native', 'UTF-8');

fprintf(f,'\n');
fprintf(f,'TABLE C1 -  Monte Carlo Simulations: Summary Statistics\n\n');
fprintf(f,'--------------------------------------------------------------------------------------------\n');
fprintf(f,'                                                     Data generating process                \n');
fprintf(f,'                                           -------------------------------------------------\n');
fprintf(f,'                                     Constrained/truncated DA   Unconstrained DA with cost  \n');
fprintf(f,'                                              (1)                        (2)                \n');
fprintf(f,'--------------------------------------------------------------------------------------------\n');
fprintf(f,'                                    Panel A. Outcomes\n\n');
fprintf(f,'Average length of submitted ROLs           %3.2f                        %3.2f\n',col1.stats.size_ROL,col2.stats.size_ROL);
fprintf(f,'                                          (%4.3f)                     (%4.3f)\n',col1.stats.size_ROL_sd,col2.stats.size_ROL_sd);
fprintf(f,'Assigned to a school                       %4.3f                       %4.3f\n',col1.stats.pct_assigned,col2.stats.pct_assigned);
fprintf(f,'                                          (%4.3f)                     (%4.3f)\n',col1.stats.pct_assigned_sd,col2.stats.pct_assigned_sd);
fprintf(f,'Weakly truth-telling                       %4.3f                       %4.3f\n',col1.stats.pct_topk,col2.stats.pct_topk);
fprintf(f,'                                          (%4.3f)                     (%4.3f)\n',col1.stats.pct_topk_sd,col2.stats.pct_topk_sd);
fprintf(f,'Assigned to favorite feasible school       %4.3f                       %4.3f\n',col1.stats.pct_assigned_pref_feas,col2.stats.pct_assigned_pref_feas);
fprintf(f,'                                          (%4.3f)                     (%4.3f)\n',col1.stats.pct_assigned_pref_feas_sd,col2.stats.pct_assigned_pref_feas_sd);
fprintf(f,'\n');
fprintf(f,'                                    Panel B. Parameters\n\n');
fprintf(f,'Number of students                         %3.0f                         %3.0f\n',col1.params.I,col2.params.I);
fprintf(f,'Number of schools                           %1.0f                           %1.0f\n',col1.params.J,col2.params.J);
fprintf(f,'Number of simulated samples                %3.0f                         %3.0f\n',col1.params.M,col2.params.M);
fprintf(f,'Maximum possible length of ROL              %1.0f                           %1.0f\n',col1.params.A,col2.params.A);
fprintf(f,'Marginal application cost (c)               %1.0f                         %g\n',col1.params.unit_cost,col2.params.unit_cost);
fprintf(f,'--------------------------------------------------------------------------------------------\n');
