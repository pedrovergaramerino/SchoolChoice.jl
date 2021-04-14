% --------------------------------------------------------------------------------------- %
% TABLE C2 - Monte Carlo Results: Unconstrained DA (500 Students, 6 Schools, 500 Samples) %
% --------------------------------------------------------------------------------------- %

clear

if exist('set_path.m', 'file') == 2
    set_path
else 
    error('ERROR: specify path to folder containing replication files in change_directory.m and execute script')  
end

% 1. Simulate data

propor = 5; % 500 students
M = 500; % 500 MC samples
J = 6; % 6 schools
A = 6; % Maximum size of ROL = 6
unit_cost = 0; % No cost of ranking more schools
MC_data_simulation; % simulate data

% 2. Estimate preferences

MC_preference_estimation;
% Save main results
Panel = RESULTS;
Panel.true_values = SETTINGS.theta_true;
Panel.stats = STATS;
clearvars -except Panel dir_tab

%%% TABLE

tab_C2 = fullfile(dir_tab,'xtable_C2.txt');
f =  fopen(tab_C2, 'wt', 'native', 'UTF-8');

parameters_lab = {'School 2','School 3','School 4','School 5','School 6','Own ability x school quality','Distance'};

fprintf(f,'\n');
fprintf(f,'TABLE C2 - Monte Carlo Simulations: Unconstrained DA (500 Students, 6 Schools, 500 Samples)\n\n');
fprintf(f,'----------------------------------------------------------------------------------------------------------\n');
fprintf(f,'                                                              Identifying assumptions                      \n');
fprintf(f,'                                           ----------------------------------------------------------------\n');
fprintf(f,'                                                  Weak                                     Stability and \n');
fprintf(f,'                                                  Truth-             Stability              undominated   \n');
fprintf(f,'                                                 telling                                    strategies   \n');
fprintf(f,'                                           ------------------     ------------------     ------------------\n');
fprintf(f,'                             True value    Mean    SD     CP      Mean    SD     CP      Mean    SD     CP\n');
fprintf(f,'                                (1)        (2)     (3)    (4)     (5)     (6)    (7)     (8)     (9)    (10)\n');
fprintf(f,[repmat('-',1,108) '\n']);
fprintf(f,'PARAMETERS\n');
for pp=1:7
fprintf(f,'%-30s % 3.2f      % 3.2f   %3.2f   %3.2f    % 3.2f   %3.2f   %3.2f    % 3.2f   %3.2f   %3.2f\n', ...
    parameters_lab{pp}, Panel.true_values(pp), ...
    Panel.TT_ML.m_est(pp), Panel.TT_ML.sd_est(pp), Panel.TT_ML.m_cp(pp), ...
    Panel.ST_ML.m_est(pp), Panel.ST_ML.sd_est(pp), Panel.ST_ML.m_cp(pp), ...
    Panel.ST_MEI.m_est(pp), Panel.ST_MEI.sd_est(pp), Panel.ST_MEI.m_cp(pp) ...
    );
end
fprintf(f,'\n');
fprintf(f,'SUMMARY STATISTICS (AVERAGED ACROSS MONTE CARLO SAMPLES)\n');
fprintf(f,'Average length of submitted ROLs                 %3.2f\n',Panel.stats.mean_nb_ranked);
fprintf(f,'Fraction of weakly truth-telling students        %3.2f\n',Panel.stats.pct_topk);
fprintf(f,'Fraction assigned to favorite feasible school    %3.2f\n',Panel.stats.pct_assigned_pref_feas);
fprintf(f,'\n');
fprintf(f,'MODEL SELECTION TESTS\n');
fprintf(f,'Truth-Telling (H0) vs. Stability (H1):           H0 rejected in %1.0f%% of samples (at 5%% significance level)\n',Panel.TT_vs_ST.m_TT_vs_ST_Hausman_rejectTT*100);
fprintf(f,'Stability (H0) vs. Undominated strategies (H1):  H0 rejected in %1.0f%% of samples (at 5%% significance level)\n',Panel.ST_vs_US.m_reject_Stability*100);
fprintf(f,'----------------------------------------------------------------------------------------------------------\n');

fclose(f);