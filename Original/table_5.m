% ---------------------------------------------------------------------------- %
% TABLE 5 - Estimation Results under Different Sets of Identifying Assumptions %
% ---------------------------------------------------------------------------- %

clear

if exist('set_path.m', 'file') == 2
    set_path
else 
    error('ERROR: specify path to folder containing replication files in change_directory.m and execute script')  
end

if dir_data == 0
    error('Cannot execute table_5.m. This program requires proprietary data - check ReadMe file');
end

% Estimate preferences under alternative assumptions
SETTINGS.estimates = {'TT_ML' 'ST_ML' 'ST_ME' 'ST_MEI'};
SETTINGS.tests = {'TT_vs_ST' 'ST_vs_US'};
SETTINGS.mci = {'MCI_ST_MEI'};

Paris_preference_estimation

%%% TABLE

tab_5 = fullfile(dir_tab,'table_5.txt');
f =  fopen(tab_5, 'wt', 'native', 'UTF-8');

parameters_lab1 = {'School 2','School 3','School 4','School 5','School 6','School 7','School 8','School 9','School 10','School 11', ...
                  'Closest school','High school co-located','Student French score [x10]' ...
                  'Student math score [x10]','High SES','Scaling parameter '};
parameters_lab2 = {'','','','','','','','','','', ...
                  '','  with middle school','  x school French score [x10]' ...
                  '  x school math score [x10]','  x fraction high SES in school',''};

fprintf(f,'\n');
fprintf(f,'TABLE 5 - Estimation Results under Different Sets of Identifying Assumptions\n\n');
fprintf(f,'----------------------------------------------------------------------------------------\n');
fprintf(f,'                                                      Identifying assumptions           \n');
fprintf(f,'                                    ----------------------------------------------------\n');
fprintf(f,'                                         Weak                            Stability and  \n');
fprintf(f,'                                        Truth-          Stability         undominated   \n');
fprintf(f,'                                        telling                           strategies    \n');
fprintf(f,'                                    ---------------   ---------------    ---------------\n');
fprintf(f,'                                          (1)               (2)               (3)  \n');
fprintf(f,'----------------------------------------------------------------------------------------\n');
fprintf(f,'\n');
fprintf(f,'                                  Panel A. School fixed effects\n\n');
for pp = 1:16
    ll=pp+2;
fprintf(f,'%-40s % 3.2f             % 3.2f             % 3.2f\n', ...
    parameters_lab1{pp}, cell2mat(RESULTS.TT_ML(ll,2)), cell2mat(RESULTS.ST_ML(ll,2)), cell2mat(RESULTS.ST_MEI(ll,2)));  
fprintf(f,'%-36s [% .2f;% .2f]     [% .2f;% .2f]     [% .2f;% .2f]\n\n', ...
    parameters_lab2{pp}, cell2mat(RESULTS.TT_ML(ll,5)),cell2mat(RESULTS.TT_ML(ll,6)), ...
        cell2mat(RESULTS.ST_ML(ll,5)),cell2mat(RESULTS.ST_ML(ll,6)), ...
        cell2mat(RESULTS.ST_MEI(ll,5)),cell2mat(RESULTS.ST_MEI(ll,6)));
if pp == 10
    fprintf(f,'\n');
    fprintf(f,'                                  Panel B.  Covariates\n\n');
end
end
fprintf(f,'\n');
fprintf(f,'Number of students                        %.0f              %0.f              %0.f\n',NB.stu,NB.assigned,NB.stu);
fprintf(f,'----------------------------------------------------------------------------------------\n');
fprintf(f,'\nTests:\n');
fprintf(f,'- TT vs. ST (Hausman): reject TT = %1.0f (test-statistic = %.1f; critical value = %.1f; p-value: %f)\n', ...
    TESTS.TT_vs_ST.Hausman_rejectTT2,TESTS.TT_vs_ST.Hausman_test2,TESTS.TT_vs_ST.crit_Hausman2,TESTS.TT_vs_ST.pvalue2);
fprintf(f,'- ST vs. US (Bugni, Canay & Shi, 2015): reject ST = %1.0f (test-statistic = %.1f; critical value = %.1f)\n', ...
    TESTS.ST_vs_US.BCS_MEI_Stability_reject,TESTS.ST_vs_US.BCS_MEI_Tn,TESTS.ST_vs_US.BCS_MEI_crit_val);
