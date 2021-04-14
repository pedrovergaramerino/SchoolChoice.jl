function Q = QStat(type,is_eq,is_ineq,theta,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,z_eq,z_ineq,Kappa,draws,is_WTT)
% QSTAT   Returns test statistics for Minimum Resampling Test (Test MR)
%         in Bugni, Canay and Shi (QE 2017), based on the pseudo code in the 
%         working paper version (Algorithm A.1). Critical values for the 
%         test statistic Tn are based on Test R1 (TR1)
%
% Inputs:
% - type: 0 for test statistic / 1 for Test R1 (TR1)
% - is_eq: 1 to compute moment equalities / 0 otherwise
% - is_ineq: 1 to compute moment inequalities / 0 otherwise
% - theta: parameter estimate at which to evaluate the function
% - is_assigned: (I,1) vector indicating whether student is assigned to a school
% - covariates: (I,N_C) matrix of covariates
% - nb_pairs: number of pairs of schools
% - pairwise: matrix of pairwise comparison of schools
% - feasible: (I,J) matrix indicating feasible schools for each student
% - e_eq: (I,N_eq) matrix of empirical probabilities for moment equalities
% - elbnd_ineq: (I,N_ineq) matrix of empirical lower bounds for moment inequalities
% - eubnd_ineq: (I,N_ineq) matrix of empirical upper bounds for moment inequalities
% - z_eq: (I,N_eq*N_Z) matrix of empirical probabilities for moment equalities, interacted with instruments Z (not used)
% - z_ineq: (I,N_ineq*N_Z) matrix of empirical probabilities for moment inequalities, interacted with instruments Z
% - Kappa: tuning parameter Kappa
% - draws: (I,N_draws) matrix of random normal draws to compute asymptotic approximation
% - is_WTT: whether to normalize estimates as willingness-to-travel (WTT)
%
% Output:
% - Q: if type=0: test statistic / if type=1: Test R1 (TR1)
%
% See also MOMENT_CONDITIONS.

% By Gabrielle Fack, Julien Grenet and Yinghua He
% 4 Oct. 2018

% NOTE: in Bugni, Canay and Shi (2014)'s pseudo code, the standardization of moment is performed through the S(.) function (line 19)
% In their application, S(.) is the Modified Method of Moments criterion function which involves dividing the moments
% by their standard deviation (i.e. does not involve covariance terms). For simplicity, the moments are here standardized 
% before plugging them into the S(.) function. This explains some differences between this code and the BCS pseudo-code

% ----------------- %
% MOMENT CONDITIONS %
% ----------------- %

% Try the fast version first (not robust to large V_ij)
[m_eq, m_ineq] = Moment_Conditions('fast',is_eq,is_ineq,theta,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,is_WTT);   

% If moments have NaN values (among assigned students only for moment equalities), try the slow version (robust to large V_ij)
if isempty(m_eq) == 0
    if (sum(sum(isnan(m_eq(is_assigned==1,:)))) > 0 || sum(sum(isnan(m_ineq))) > 0)
        [m_eq, m_ineq] = Moment_Conditions('slow',is_eq,is_ineq,theta,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,is_WTT);   
    end
end
if isempty(m_eq) == 1
    if sum(sum(isnan(m_ineq))) > 0
        [m_eq, m_ineq] = Moment_Conditions('slow',is_eq,is_ineq,theta,is_assigned,covariates,nb_pairs,pairwise,feasible,e_eq,elbnd_ineq,eubnd_ineq,is_WTT);   
    end
end

% Expand moments if instruments are used
if isempty(z_ineq)==0 
    ratio_ineq = size(z_ineq,2)/size(m_ineq,2);
    m_ineq = repmat(m_ineq, [1 ratio_ineq]).*z_ineq;
end

% I: # of students (whether assigned or unassigned
I = length(is_assigned);

% N_eq = number of moment equalities
N_eq = size(m_eq,2);

% N_ineq = number of moment inequalities
N_ineq = size(m_ineq,2);

% ----------------------------- %
% STANDARDIZE MOMENT CONDITIONS %
% ----------------------------- %    

% MOMENT EQUALITIES %
% ----------------- %

% Note: for moment equalities, we don't use unassigned students. The corresponding values in m_eq are set to NaN

if is_eq == 1
    m_eq_NaN = m_eq;  
    % replace moments for non feasible schools by NaNs
    if is_ineq == 0
        m_eq_NaN(not(feasible(:,2:end))) = NaN;
    % If moment inequalities are supplied, keep first redundant moment equality
    elseif is_ineq == 1
        m_eq_NaN(not(feasible))=NaN;
    end  
    % s.d. of moment equalities excluding NaN
    m_eq_NaN_std = nanstd(m_eq_NaN,1,1);  
    % Standardized moment equalities
    mbar_NaN_std_eq = nanmean(m_eq_NaN)./m_eq_NaN_std;
    % Number of "relevant" students for each moment equality (row vector)
    I_eq_NaN = sum(not(isnan(m_eq_NaN)),1);         
end

% MOMENT INEQUALITIES %
% ------------------- %

% N.B.: contrary to moment equalities, moment inequalities are defined for all assigned students and never have NaNs

if is_ineq ==1
    % s.d. of moment inequalities
    Sigma_ineq = std(m_ineq,1,1);
    % Replace zeros by very small number
    Sigma_ineq(Sigma_ineq == 0) = 1e-20;     
    % Mean of moment inequalities
    mbar_ineq = mean(m_ineq,1);
    % Standardized moment inequalities
    mbar_std_ineq  = mbar_ineq./Sigma_ineq;
end

% ----------------------------- %
% COMPUTE V(THETA) and L(THETA) %
% ----------------------------- %    

% TYPE = 0 (TEST STATISTIC) %
% ------------------------- %

if type == 0
    % Moment equalities
    if is_eq == 1
        % v(theta)
        v_theta_eq = sqrt(I_eq_NaN).*mbar_NaN_std_eq; % scaled average (line 10 in BCS). N.B.: standardization is done here. In BCS, it's done on line 19.
        % l(theta)
        l_theta_eq = zeros(1,N_eq); % test statistic does not involve l(.) (line 11 in BCS)
    end
    % Moment inequalities
    if is_ineq ==1
        % v(theta)
        v_theta_ineq = sqrt(I).*mbar_std_ineq; % scaled average (line 10 in BCS)
        % l(theta)
        l_theta_ineq = zeros(1,N_ineq); % test statistic does not involve l(.) (line 11 in BCS)
    end          
end

% TYPE = 1 (TEST R1) %
% ------------------ %

if type == 1
    % Moment equalities
    if is_eq == 1       
        % v(theta): use asymptotic approximation
        xmu=nanmean(m_eq_NaN,1);
        zscore_NaN=(m_eq_NaN-repmat(xmu,I,1))./repmat(m_eq_NaN_std,I,1);
        % For cells with NaN values, put to zero to ignore in computation of v_theta_eq           
        zscore_NaN(isnan(zscore_NaN))=0;
        v_theta_eq = (ones(1,N_eq)./sqrt(I_eq_NaN)).*(draws'*zscore_NaN); % define stochastic process (line 13 in BCS)    
        % l(theta)
        l_theta_eq = zeros(1,N_eq); % N.B.: Phi = 0 for moment equalities (see p. 6 of BCS)
    end
    % 2/ Moment inequalities
    if is_ineq == 1
        % v(theta): use asymptotic approximation
        zscore_ineq = zscore(m_ineq,1);
        v_theta_ineq = (1/sqrt(I)).*(draws'*zscore_ineq); % define stochastic process (line 13 in BCS)
        % l(theta)
        B_n = sqrt(0.4*log(I)/log(log(I)));
        l_theta_ineq = ((sqrt(I).*mbar_std_ineq./Kappa) >=1 ).*10^10;%B_n; % BCS use 10^10 in their code. Using B_n instead makes no difference (line 14 in BCS)
    end                     
end

% ------------ %
% RETURN QSTAT %
% ------------ %    

% Q-statistic is based on Modified Methods of Moments criterion
Q = 0;

% Moment equalities
if is_eq == 1
    % Print error if some moments equalities are undefined
    if sum(sum(isnan(l_theta_eq))) > 0 || sum(sum(isnan(v_theta_eq))) > 0
        disp(mat2str(theta))
        error('ERROR: m_eq undefined')          
    else
    % Sum of squares of moment equalities
        Q_eq = sum((l_theta_eq + v_theta_eq).^2,2);
        Q = Q + Q_eq; 
    end
end

% Moment inequalities
if is_ineq == 1
    % Print error if some moments inequalities are undefined
    if sum(sum(isnan(l_theta_ineq)))>0 || sum(sum(isnan(v_theta_ineq)))>0
        error('ERROR: m_ineq undefined')    
    else
    % Sum of squares of moment inequalities (if >0)
        Q_ineq = sum((min(0,(l_theta_ineq + v_theta_ineq))).^2,2);
        Q = Q + Q_ineq;      
    end
end