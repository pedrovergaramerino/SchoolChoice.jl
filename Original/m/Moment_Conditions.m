function [m_eq, m_ineq] = Moment_Conditions(speed,is_eq,is_ineq,theta,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,is_WTT)
% MOMENT_CONDITIONS   Computes moment equalities and inequalities at a given value of theta
%
% Inputs:
% - speed: 'fast' for fast version (not robust to large V_ij) / 'slow' for % slow version (robust to large V_ij)
% - is_eq: 1 to compute moment equalities / 0 otherwise
% - is_ineq: 1 to compute moment inequalities / 0 otherwise
% - theta: vector of parameters at which the moment conditions are evaluated
% - is_assigned: (I,1) vector indicating whether student is assigned to a school
% - covariates: (I,N_C) matrix of covariates
% - nb_pairs: number of pairs of schools
% - pairwise: matrix of pairwise comparison of schools
% - feasible: (I,J) matrix indicating feasible schools for each student
% - e_eq: (I,N_eq) matrix of empirical probabilities for moment equalities
% - elbnd_ineq: (I,N_eq) matrix of empirical lower bounds for moment inequalities
% - eubnd_ineq: (I,N_ineq) of empirical upper bounds for moment inequalities
% - is_WTT: whether estimates are normalized to be in willingness-to-travel (WTT) form
%
%  Output:
% - m_eq: (I,N_eq) matrix of moment inequalities 
% - m_ineq: (I,N_ineq) matrix of moment inequalities

% By Gabrielle Fack, Julien Grenet and Yinghua He
% 4 Oct. 2018

% ---------- %
% PARAMETERS %
% ---------- %

% # of parameters
nb_params = size(theta,1);  

% Add index to matrix of pairwise comparison of schools
pairwise2 = [(1:nb_pairs)', pairwise];

% I: # of students
I = size(is_assigned,1);

% J: # of schools
J = size(feasible,2);

% ----------------------------------------------------------------------------------------------- %
%                                     PREDICTED PROBABILITIES                                     %
% ----------------------------------------------------------------------------------------------- %

% Compute predicted probabilities that will be used to form moment conditions:
% p_eq: predicted probabilities for moment equalities. 2 components:
%       - predicted school assignment probability ((I,J) matrix)
%       - predicted distance/own score x school score at assigned school ((I,2) matrix)
% p_ineq: predicted choice probabilities for every pairwise comparison of schools ((I,K) matrix, where K = J!/(2!(J-2)!))

% --------------------------------------------------- %
% P_eq: PREDICTED PROBABILITIES FOR MOMENT EQUALITIES %
% --------------------------------------------------- %

% Predicted school assignment probabilities ((I,J) matrix) %
% -------------------------------------------------------- %

% (I,J) matrix where (i,j) is the predicted assignment of student i to school j
% - if the school is feasible, it is the probability that j is student i's preferred school among the feasible
% - if school is not feasible, the probability is zero

if is_WTT == 1
    if abs(theta(end)) <= 1e-5
        disp('Warning: sigma = zero')
        theta(end)=0.01;
    end    
    theta(1:end-1) = theta(1:end-1)./theta(end);
    theta(end) = -1./theta(end);
end   

% FAST VERSION (not robust to large V_ij)
if strcmp(speed,'fast') == 1
    
    % V_ij: Fixed component of utility ((I,J) matrix)
    V_ij = repmat([0, theta(1:J-1)'],I,1);
    
    if nb_params>J-1
        for pp=J:nb_params
            V_ij = V_ij + theta(pp).*covariates(:,:,pp+1);
        end
    end 
    
    % exp(V_ij)
    expV = exp(V_ij);

    % Predicted assignment among feasible schools: J-1 theoretical moments
    % Note: NaNs are possible if no school is feasible
    predicted = (expV.*feasible)./repmat((expV.*feasible )*ones(J,1),1,J);

end

% SLOW VERSION (robust to large V_ij (avoids infty/infty) problem)
if strcmp(speed,'slow') == 1
    predicted = NaN(I,J);

    % V_ij: Fixed component of utility (N x J matrix)
    theta_fe = [0, theta(1:J-1)'];
    if is_ineq == 1
        p_ineq = zeros(I,nb_pairs);
    end
    % Instead of computing V_ij for all j=1...J, compute V_ik - Vij for all k=1...J, j=1,...J
    % Then P(i chooses j) = 1 / [ 1 + sum_k exp(V_ik - V_ij) ]
    for jj = 1:J
        V_diff = repmat(theta_fe, [I 1])-repmat(theta_fe(1,jj),[I J]);
        if nb_params>J-1
            for pp=J:nb_params
                V_diff = V_diff + theta(pp).*(covariates(:,:,pp+1) - repmat(covariates(:,jj,pp+1),[1 J]));
            end
        end               
        expV_diff = exp(V_diff);
        expV_diff(isinf(expV_diff))=10.^20;  % As precaution when exp(V) is too large for machine 
        % Predicted assignment among feasible schools: J-1 theoretical moments
        predicted(:,jj) = feasible(:,jj)./sum(expV_diff.*feasible,2);
        % Predicted assignment for non-feasible schools is set to zero
        predicted(feasible(:,jj)==0,jj) = 0;
        % Predicted choice probability for pairwise comparisons (for moment inequalities) 
        if is_ineq == 1
           % find corresponding pair
           ind_cc = pairwise2(pairwise2(:,3)==jj,2)';
           num_cc = pairwise2(pairwise2(:,3)==jj,1)';
           p_ineq(:,num_cc) = 1./(1+expV_diff(:,ind_cc));          
        end        
    end
end

% Predicted distance/own score x school score of assigned school ((I,2) matrix) %
% ----------------------------------------------------------------------------- %

p_eq = predicted; 
if nb_params>J-1
    for pp=J:nb_params
        p_eq = [p_eq, sum(predicted.*covariates(:,:,pp+1),2)];
    end            
end

p_eq(is_assigned==0,:) = NaN;

% ------------------------------------------------------- %
% P_ineq: PREDICTED PROBABILITIES FOR MOMENT INEQUALITIES %
% ------------------------------------------------------- %

% FAST VERSION
% Note: to speed up computation (p_ineq in the slow version is above)
if strcmp(speed,'fast')==1 && is_ineq ==1
    p_ineq = zeros(I,nb_pairs);
    for cc=1:nb_pairs
        % For each pairwise comparison, compute probability of second school being preferred to the first
        p_ineq(:,cc) = expV(:,pairwise(cc,2))./(expV(:,pairwise(cc,1)) + expV(:,pairwise(cc,2)));
    end   
end

% ----------------------------------------------------------------------------------------------- %
%                                       MOMENT CONDITIONS                                         %
% ----------------------------------------------------------------------------------------------- %

% Moment conditions = empirical probabilities - predicted probabilities

% ----------------- %
% MOMENT EQUALITIES %
% ----------------- %

% Note: # of observatins for moment equalities is # of assigned students

if is_eq == 0
    % output empty matrix if no moment equalities
    m_eq=[];
elseif is_eq == 1                         
    % If ME+MI -> keep redundant moment
    if is_ineq == 1
        m_eq = e_eq(:,1:end) - p_eq(:,1:end);
    end
    % If ME only -> remove redundant moment
    if is_ineq == 0
        m_eq = e_eq(:,2:end) - p_eq(:,2:end);
    end
end

% ------------------- %
% MOMENT INEQUALITIES %
% ------------------- %

% Note: # of observatins for moment inequalities is # of students (whether assigned or not)

if is_ineq == 1
    m_ineq = [p_ineq - elbnd_ineq , eubnd_ineq - p_ineq];
elseif is_ineq == 0
    % output empty matrix if no moment inequalities are supplied
    m_ineq=[];
end