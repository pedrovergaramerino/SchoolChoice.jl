% ---------------------------------------------------------- %
% TABLE D2 - Goodness of Fit: Observed vs. Simulated Cutoffs %
% ---------------------------------------------------------- %

clear

if exist('set_path.m', 'file') == 2
    set_path
else 
    error('ERROR: specify path to folder containing replication files in change_directory.m and execute script')  
end

if dir_data == 0
    error('Cannot execute table_D2.m. This program requires proprietary data - check ReadMe file');
end

% Estimate preferences under alternative assumptions
SETTINGS.estimates = {'TT_ML' 'ST_ML' 'ST_MEI'};
SETTINGS.tests = {};
SETTINGS.mci = {};

Paris_preference_estimation; % estimate preferences

Goodness_of_fit; % run goodness of fit


%%% TABLE

tab_D2 = fullfile(dir_tab,'xtable_D2.txt');
f =  fopen(tab_D2, 'wt', 'native', 'UTF-8');

fprintf(f,'\n');
fprintf(f,'TABLE D2 - Goodness of Fit: Observed vs. Simulated Cutoffs\n\n');
fprintf(f,'----------------------------------------------------------------------------------------------------\n');
fprintf(f,'                                               Cutoffs in simulated samples with estimates from     \n');
fprintf(f,'                                           ---------------------------------------------------------\n');
fprintf(f,'                               Observed               Weak          Stability of    Stability and   \n');
fprintf(f,'                                cutoffs          Truth-telling      the matching     undominated    \n');
fprintf(f,'                                                                                     strategies     \n');
fprintf(f,'                                 (1)                   (2)              (3)              (4)        \n');
fprintf(f,'----------------------------------------------------------------------------------------------------\n');
for ss=1:11
fprintf(f,'School %2.0f                       %4.3f                 %4.3f            %4.3f             %4.3f\n', ...
    ss,SCH.cutoffs(ss),GOF.TT_ML.m_cutoffs(ss),GOF.ST_ML.m_cutoffs(ss),GOF.ST_MEI.m_cutoffs(ss));
fprintf(f,'                                                     (%4.3f)          (%4.3f)           (%4.3f)\n', ...
    GOF.TT_ML.sd_cutoffs(ss),GOF.ST_ML.sd_cutoffs(ss),GOF.ST_MEI.sd_cutoffs(ss));  
end
fprintf(f,'----------------------------------------------------------------------------------------------------\n');
