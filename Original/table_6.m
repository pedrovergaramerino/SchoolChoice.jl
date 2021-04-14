% ------------------------------------------------------------------------------------- %
% TABLE 6 - Goodness-of-Fit Measures Based on Different Sets of Identifying Assumptions %
% ------------------------------------------------------------------------------------- %

clear

if exist('set_path.m', 'file') == 2
    set_path
else 
    error('ERROR: specify path to folder containing replication files in change_directory.m and execute script')  
end

if dir_data == 0
    error('Cannot execute table_6.m. This program requires proprietary data - check ReadMe file');
end

% Estimate preferences under alternative assumptions
SETTINGS.estimates = {'TT_ML' 'ST_ML' 'ST_MEI'};
SETTINGS.tests = {};
SETTINGS.mci = {};

Paris_preference_estimation; % estimate preferences

Goodness_of_fit; % run goodness of fit

%%% TABLE

tab_6 = fullfile(dir_tab,'table_6.txt');
f =  fopen(tab_6, 'wt', 'native', 'UTF-8');

fprintf(f,'\n');
fprintf(f,'TABLE 6 - Goodness-of-Fit Measures Based on Different Sets of Identifying Assumptions\n\n');
fprintf(f,'----------------------------------------------------------------------------------------------------\n');
fprintf(f,'                                                                  Estimates from                    \n');
fprintf(f,'                                      --------------------------------------------------------------\n');
fprintf(f,'                                                      Weak          Stability of    Stability and   \n');
fprintf(f,'                                                 Truth-telling      the matching     undominated    \n');
fprintf(f,'                                                                                     strategies     \n');
fprintf(f,'                                                       (1)              (2)              (3)        \n');
fprintf(f,'----------------------------------------------------------------------------------------------------\n');
fprintf(f,'\n');
fprintf(f,'      Panel A. Simulated vs. observed assignment (300 simulated samples) \n\n');
fprintf(f,'Mean predicted fraction of students                   %4.3f            %4.3f            %4.3f\n', ...
    GOF.TT_ML.prob_assignment,GOF.ST_ML.prob_assignment,GOF.ST_MEI.prob_assignment);
fprintf(f,'  assigned to observed assignment                    (%4.3f)          (%4.3f)          (%4.3f)\n', ...
    GOF.TT_ML.prob_assignment_sd,GOF.ST_ML.prob_assignment_sd,GOF.ST_MEI.prob_assignment_sd);
fprintf(f,'\n');
fprintf(f,'      Panel B. Predicted vs. observed partial preference order of given schools) \n\n');
fprintf(f,'Mean predicted probability that a student\n');
fprintf(f,'  prefers the top-ranked school to the                %4.3f            %4.3f            %4.3f\n', ...
    GOF.TT_ML.pred_ordering_top2_choices,GOF.ST_ML.pred_ordering_top2_choices,GOF.ST_MEI.pred_ordering_top2_choices);
fprintf(f,'  2nd-ranked in her submitted ROL\n');
fprintf(f,'\n');
fprintf(f,'Mean predicted probability that a student''s partial\n');
fprintf(f,'  preference order among the schools in her ROL       %4.3f            %4.3f            %4.3f\n', ...
    GOF.TT_ML.pred_ordering_all_choices,GOF.ST_ML.pred_ordering_all_choices,GOF.ST_MEI.pred_ordering_all_choices);
fprintf(f,'  coincides with the submitted rank order\n');
fprintf(f,'----------------------------------------------------------------------------------------------------\n');