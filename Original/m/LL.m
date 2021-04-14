function out = LL(theta,choices,choice_set,covariates,is_WTT)
% LL    Log-likelihood
%
% Inputs:
% - theta: parameter estimate at which to evaluate the function
% - choices: (I,J) matrix of choices (1 if school is chosen, 0 otherwise)
% - choice_set: (I,J) choice set matrix (1 if school in choice set, 0 otherwise)
% - covariates: (I,N_C) matrix of covariates
% - is_WTT: whether to normalize estimates as willingness-to-travel (WTT)
%
%  Output:
% - out: log-likelihood function evaluated at theta

% By Gabrielle Fack, Julien Grenet and Yinghua He
% 4 Oct. 2018

% -------------------- %
% CHOICE PROBABILITIES %
% -------------------- %

% # of parameters
nb_params = size(theta,1);

% # of (potentially exploded) observations
Nobs = size(choices,1);

% # of schools
J = size(choice_set,2);

% WTT: distance parameter = -1
if is_WTT == 1
    if abs(theta(end)) <= 1e-5
        disp('Warning: sigma = zero')
        theta(end) = 0.01;
    end    
    theta(1:end-1) = theta(1:end-1)./theta(end);
    theta(end) = -1./theta(end);
end   

% V_ij: Fixed component of utility (I x J matrix)
V_ij = repmat([0, theta(1:J-1)'], [Nobs 1]);
replication_factor = size(choices,1)/size(covariates,1);
if nb_params>J-1
    for pp=J:nb_params          
        V_ij = V_ij + theta(pp).*repmat(covariates(:,:,pp+1),replication_factor,1);
    end       
end

% exp(V_ij)
expV = exp(V_ij);

% Probability that each school is preferred in the choice set (I x J matrix)
prob = expV.*choice_set./repmat((expV.*choice_set)*ones(J,1), [1 J]);

% Log of choice probability for chosen school
log_prob = log(prob(choices==1));

% ------ %
% OUTPUT %
% ------ %
   
% Compute log-likelihood by summing over individuals (take minus for minimization)
out = -sum(log_prob);   