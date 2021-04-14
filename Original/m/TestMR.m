function out = TestMR(type,is_eq,is_ineq, position,theta_est,theta_i_0,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,z_eq,z_ineq,Kappa,conf_level,eta,draws,is_WTT)
% TESTMR   Tests whether the value theta(i)0 of an individual coordinate theta(i) 
%          of theta belongs to the marginal confidence interval (mci) of that coordinate.         
%          This test is used to construct marginal confidence intervals for
%          each coordinate of theta. It implements the Minimum Resampling Test (Test MR)
%          in Bugni, Canay & Shi (QE 2017), based on the pseudo code in the 
%          working paper version (Algorithm A.1). Critical values for the test statistic 
%          Tn are based on Test R1 (TR1)
%
% Inputs:
% - type: 
%   'cp': if objective is to test whether a particular value is in theta(i)'s mci
%   'lb': if objective is to find lower bound of theta(i)'s mci
%   'ub': if objective is to find lower bound of theta(i)'s mci
% - is_eq: 1 to compute moment equalities / 0 otherwise
% - is_ineq: 1 to compute moment inequalities / 0 otherwise
% - position: coordinate of individual coordinate theta(i) whose mci is computed
% - theta_est: vector of point estimates for theta
% - theta_i_0: tentative value of the i-th coordinate of theta whose marginal confidence interval is being constructed
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
% - Kappa: tuning parameter Kappa
% - conf_level: confidence level adjusted by eta
% - eta: tuning parameter eta
% - draws: (I,N_draws) matrix of random normal draws to compute asymptotic approximation
% - is_WTT: whether to normalize estimates as willingness-to-travel (WTT)
%
% Output:
% - out: 
%   if type='cp': whether theta(i)_0 is in theta(i)'s mci
%   if type='lb': value of objective function for finding lower bound, evaluated at theta(i)_0
%   if type='ub': value of objective function for finding upper bound, evaluated at theta(i)_0
%
% See also: QSTAT

% By Gabrielle Fack, Julien Grenet and Yinghua He
% 4 Oct. 2018
                               
% Notations:
% - theta is the generic notation for theta
% - theta_i is the generic notation for theta(i), the value of the i-th coordinate of theta
% - theta_mi is the generic notation for theta(-i), the vector of values of all coordinates of theta other than theta(i)
% - theta_est: is the vector of point estimates for theta
% - theta_i_0 is the tentative value theta(i)0 of the i-th coordinate of theta whose marginal confidence interval is being constructed
% - theta_mi_hat is the mimizer theta(-i)hat of the test statistic Tn when the value of theta(i) is fixed at theta(i)0

% Position of free parameters
free_pos = [1:size(theta_est,1)];
free_pos(position) = [];

% ----------------------------------------------------------------------------------------------- %
%                                     TUNING PARAMETERS                                           %
% ----------------------------------------------------------------------------------------------- %

% -------------------- %
% Optimization options %
% -------------------- %

optfminunc_TMR = optimoptions(@fminunc,'display','off','Algorithm','quasi-newton','MaxIter',5e+4,'MaxFunEvals',5e+4,'TolFun',1e-2,'TolX',1e-2);
                                        
% ----------------------------------------------------------------------------------------------- %
%                   COMPUTE TEST STATIC (Tn) AND ASSOCIATED MINIMIZER (THETA(-i))                 %
% ----------------------------------------------------------------------------------------------- %

% For a tentative value of theta(i), the test statistic Tn is the min of QStat(0) over theta(-i) while keeping theta(i) fixed
% N.B.: Tn is the same for all i in theta(i)

% --------------------------- %
% Objective function QStat(0) %
% --------------------------- %

% Objective function: value of QStat (type=0) for a given theta(i)0 (which is fixed) and theta(-i) (over which QStat will be minimized)
function QStat_0 = obj_fun_Tn(theta_i, theta_mi)
    % Reconstruct full theta vector
    % Start with empty theta
    theta_full = NaN(size(theta_est,1),1);
    % Add theta(i), which is fixed
    theta_full(position) = theta_i; 
    % Add free parameters: theta(-i)
    theta_full(free_pos) = theta_mi;
    % Compute Qstat (type = 0) for theta = (theta(i),theta(-i))
    QStat_0 = QStat(0,is_eq,is_ineq,theta_full,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,z_eq,z_ineq,[],[],is_WTT);
end
% Create handle function
fun = @(x)obj_fun_Tn(theta_i_0,x);

% ------------------------------------------------------------------------------------ %
% Find Tn = minimand of QStat(0) over theta(-i) and associated minimizer theta(-i)hat  %
% ------------------------------------------------------------------------------------ %

% Minimize Tn over theta(-i)

% Starting value for theta(-i): estimate of theta(-i)
theta_mi_0 = theta_est;
theta_mi_0(position) = []; % delete theta(i)

if size(theta_est,1)==1 % If the model has only one free parameter
    Tn = QStat(0,extra_term,is_eq,is_ineq,theta_i_0,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,z_eq,z_ineq,[],[],is_WTT);
else % If the model has more than one free parameter
    % Tn: minimand of QStat(0) over theta(-i) [line 22 of BCS]
    % theta_mi_hat: minizer of QStat(0) over theta(-i) [line 23 of BCS]
    [theta_mi_hat, Tn, ~] = fminunc(fun,theta_mi_0,optfminunc_TMR);
end

% ----------------------------------------------------------------------------------------------- %
%             COMPUTE  TMR(b) USING ASYMPTOTIC APPROXIMATION (N_draws NORMAL DRAWS)               %
% ----------------------------------------------------------------------------------------------- %

% ------------------------------------ %
% Objective function for TR1: QStat(1) %
% ------------------------------------ %

% Objective function for TR1: if theta(-i)hat is point identified, TR1 is simply obtained by plugging theta(i) and theta(-i)hat into QStat(1)

% ----------------------------- %
% Compute TR1 over normal draws %
% ----------------------------- %

% Only compute TR1 if Tn>0
if Tn>0
    
% # of draws    
ndraws = size(draws,2);
   
% Initialize values of TR1 and TMR over normal draws     
TR1 = NaN(1,ndraws);
TMR = NaN(1,ndraws);

% Empty theta
theta = NaN(size(theta_est,1),1);
% Add theta(i)_0, which is fixed
theta(position) = theta_i_0;  
% Add free parameters: theta(-i)
theta(free_pos) = theta_mi_hat;  
    
% Loop over normal draws
for bb = 1:ndraws
    
    % Extract relevant draws
    draws_bb = draws(:,bb);
        
    % ------------------ %
    % TEST STATISTIC TMR %
    % ------------------ %

    % If the model is point identified (ME+MI): theta(-i)hat is a vector
    % TR1 is obtained by evaluating QStat(1) at (theta(i)est, theta(-i)hat)
    if (is_eq == 1 || is_ineq == 1)
        % If the model has only one free parameter
        if  size(theta_est,1)==1
            TR1(bb) = QStat(1,is_eq,is_ineq,theta_i_0,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,z_eq,z_ineq,Kappa,draws_bb,is_WTT);
        % If the model has more than one free parameter
        else
            TR1(bb) = QStat(1,is_eq,is_ineq,theta,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,z_eq,z_ineq,Kappa,draws_bb,is_WTT);
        end
    end
    % If the model is not point identified (MI only), print error
    if (is_eq == 0 && is_ineq == 1)
        error('ERROR: the model is not point identified')  
    end  

    TMR(bb) = TR1(bb);
end

% ----------------------------------------------------------------------------------------------- %
%                               COMPUTE CRITICAL VALUE FOR TEST MR                                %
% ----------------------------------------------------------------------------------------------- %

% Critical value for TMR
crit_val = prctile(TMR,conf_level) + eta;

% test takes value 1 if theta(i) belongs to the confidence set of theta and 0 otherwise
test = (Tn <= crit_val);

% ----------------------------------------------------------------------------------------------- %
%                                        OUTPUT                                                   %
% ----------------------------------------------------------------------------------------------- %

% --------------------------------------------------------------------------------------------- %
% Simple test: whether a particular value is in the marginal confidence interval of theta(i)est %
% --------------------------------------------------------------------------------------------- %

if strcmp(type,'cp') == 1
    out = test;
end

% ----------------------------------------------------------- %
% Upper bound of marginal confidence interval for theta(i)est %
% ----------------------------------------------------------- %

% Objective function: find the largest theta(i) that belongs to the confidence set of theta
if strcmp(type,'ub') == 1
    out = -theta_i_0 + (10^6)*max(Tn-crit_val, 0);
end

% ----------------------------------------------------------- %
% Lower bound of marginal confidence interval for theta(i)est %
% ----------------------------------------------------------- %

% Objective function: find the smaller theta(i) that belongs to the confidence set of theta
if strcmp(type,'lb') == 1
    out = theta_i_0 + (10^6)*max(Tn-crit_val, 0);
end

end % end of case: Tn>0

if Tn==0 % Case: Tn=0
    % Tn=0 -> theta(i) belongs to confidence set of theta
    % Upper bound:  continue searching upwards
    if strcmp(type,'ub') == 1
        out = -theta_i_0 ;
    % Lower bound: continue searching downwards
    elseif strcmp(type,'lb') == 1
        out = theta_i_0;
    elseif strcmp(type,'cp') == 1
        out = 1;       
    end
end

end