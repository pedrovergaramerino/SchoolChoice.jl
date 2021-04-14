function [test, Tn, crit_val_TR1]  = Stability_test(is_eq,is_ineq,theta_0,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,z_eq,z_ineq,Kappa,conf_level,eta,draws,is_WTT)
% STABILITY_TEST   Tests whether theta_0 belongs to the confidence set of theta, using 
%                  Test RS in Bugni, Canay and Shi (JoE 2015). This test is used to test 
%                  stability + undominated strategy against undominated strategies only
%
% Inputs:
% - is_eq: 1 to compute moment equalities / 0 otherwise
% - is_ineq: 1 to compute moment inequalities / 0 otherwise
% - theta_0: value of parameters
% - is_assigned: (I,1) vector indicating whether student is assigned to a school
% - covariates: (I,N_C) matrix of covariates
% - nb_pairs: number of pairs of schools
% - pairwise: matrix of pairwise comparison of schools
% - feasible: (I,J) matrix indicating feasible schools for each student
% - e_eq: (I,N_eq) matrix of empirical probabilities for moment equalities
% - elbnd_ineq: (I,N_ineq) matrix of empirical lower bounds for moment inequalities
% - eubnd_ineq: (I,N_ineq) matrix of empirical upper bounds for moment inequalities
% - z_eq: (I,N_eq*N_Z) matrix of empirical probabilities for moment equalities, interacted with instruments Z
% - z_ineq: (I,N_ineq*N_Z) matrix of empirical probabilities for moment inequalities, interacted with instruments Z
% - Kappa: tuning parameter
% - conf_level: confidence level adjusted by eta
% - eta: tuning parameter eta
% - draws: (I,N_draws) matrix of random normal draws to compute asymptotic approximation
% - is_WTT: whether to normalize estimates as willingness-to-travel (WTT)
%
% Output:
% - test: returns 0 if theta_0 is not in confidence set of theta (stability rejected) / 1 if theta_0 is in the confidence set (stability not rejected) 
% - Tn: test statistic Tn
% - crit_val_TR1: critical values based on Test R1
%
% See also: QSTAT

% By Gabrielle Fack, Julien Grenet and Yinghua He
% 4 Oct. 2018

% ------------------ %
% TEST STATISTIC Tn  %
% ------------------ %

% Tn = value of QStat(0) at theta_0
Tn = QStat(0,is_eq,is_ineq,theta_0,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,z_eq,z_ineq,[],[],is_WTT);

% ------------------------------------- %
% COMPUTE CRITICAL VALUES BASES ON TR1  %
% ------------------------------------- %

TR1 = NaN(1,size(draws,2));

parfor bb = 1:size(draws,2)
    
    % Random draws of N(0,1)
    draws_bb = draws(:,bb);
    
    % Test R1 statistic
    TR1(bb) = QStat(1,is_eq,is_ineq,theta_0,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,z_eq,z_ineq,Kappa,draws_bb,is_WTT);
    
end

% Critical value
crit_val_TR1 = prctile(TR1,conf_level)+eta;

% Inequality constraint: Tn - crit_val_TR1 <= 0
test = (Tn <= crit_val_TR1);

end