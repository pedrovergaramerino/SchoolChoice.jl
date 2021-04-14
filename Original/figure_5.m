% ---------------------------------------------------------- %
% FIGURE 5 - Goodness of Fit: Observed vs. Simulated Cutoffs %
% ---------------------------------------------------------- %

clear

if exist('set_path.m', 'file') == 2
    set_path
else 
    error('ERROR: specify path to folder containing replication files in change_directory.m and execute script')  
end

if dir_data == 0
    error('Cannot execute figure_5.m. This program requires proprietary data - check ReadMe file');
end

% Estimate preferences under alternative assumptions
SETTINGS.estimates = {'TT_ML' 'ST_ML' 'ST_MEI'};
SETTINGS.tests = {};
SETTINGS.mci = {};

Paris_preference_estimation; % estimate preferences

Goodness_of_fit; % run goodness of fit

% Export figures

fig_5 = table((1:11)',SCH.cutoffs, GOF.TT_ML.m_cutoffs, GOF.ST_ML.m_cutoffs, GOF.ST_MEI.m_cutoffs, 'VariableNames',{'id', 'observed' 'TT' 'ST' 'MEI'});

writetable(fig_5,fullfile(dir_fig,'figure_5.xlsx'),'Sheet',1,'Range','A1');

