% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
% ESTIMATION OF STUDENT PREFERENCES UNDER CONSTRAINED DA / DA WITH COST - MONTE CARLO SIMULATIONS %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% By Gabrielle Fack, Julien Grenet and Yinghua He
% 4 Oct. 2018

fprintf('\n-------------------------------------------------------------------------------------------------\n');
fprintf(' ESTIMATION OF STUDENT PREFERENCES UNDER CONSTRAINED DA / DA WITH COST - MONTE CARLO SIMULATIONS ');
fprintf('\n-------------------------------------------------------------------------------------------------\n');
fprintf(' # By Gabrielle Fack, Julien Grenet and Yinghua He\n');

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                     TUNING PARAMETERS                                           %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% Specify here whether coefficients should be normalized as willingness to travel (WTT)
SETTINGS.is_WTT = 0;

% --------------------- %
% True parameter values %
% --------------------- %

SETTINGS.theta_true = (PARAM.FE - PARAM.FE(1))';
SETTINGS.theta_true = [SETTINGS.theta_true(2:end); PARAM.coeff_score; PARAM.coeff_dist];

% # of parameters
SETTINGS.nb_param = size(SETTINGS.theta_true,1);

% ------------------------------- %
% Specify which results to output %
% ------------------------------- %

% Estimation %
% ---------- %

%   'TT_ML'         estimates based on weak truth-telling (MLE)
%   'ST_ML'         estimates based on stability (MLE)
%   'ST_ME'         estimates based on stability (Moment equalities)
%   'ST_MEI'        estimates based on stability + undominated strategies (Moment equalities + inequalities)

SETTINGS.estimates = {'TT_ML' 'ST_ML' 'ST_MEI'};  

% Tests %
% ----- %

%   'Test_TT_ST'    Hausman test for Truth-telling vs. Stability
%   'Test_ST_US'    Hausman test for Stability vs. Undominated Strategies

SETTINGS.tests = {'TT_vs_ST' 'ST_vs_US'};

% 95% coverage probabilities %
% -------------------------- %

SETTINGS.cp = {'ST_MEI_CP'};

% Internal checks %
% --------------- %

if  (sum(strcmp('TT_vs_ST',SETTINGS.tests))==1 )
    % TT estimates required
    if sum(strcmp('TT_ML',SETTINGS.estimates))~=1
        SETTINGS.estimates = [SETTINGS.estimates 'TT_ML'];
    end
    % ST estimates required
    if sum(strcmp('ST_ML',SETTINGS.estimates))~=1
        SETTINGS.estimates = [SETTINGS.estimates 'ST_ML'];
    end
end

% MEI estimates required
if  ( sum(strcmp('ST_vs_US',SETTINGS.tests))==1  ) 
    if sum(strcmp('ST_MEI',SETTINGS.estimates))~=1
        SETTINGS.estimates = [SETTINGS.estimates 'ST_MEI'];
    end
end

% ------------------------------%
% Initialize estimation results %
% ------------------------------%
 
init_zeros = NaN(SETTINGS.nb_param,PARAM.M);
init_zerosH = NaN(SETTINGS.nb_param,SETTINGS.nb_param,PARAM.M);

% True preferences (benchmark)
True_ML_est = init_zeros;   
True_ML_se = init_zeros;  
True_ML_hessian = init_zerosH;
True_ML_cp = init_zeros;   

% Truth-telling
TT_ML_est = init_zeros;   
TT_ML_se = init_zeros;  
TT_ML_hessian = init_zerosH;
TT_ML_cp = init_zeros;   

% Stability (MLE)
ST_ML_est = init_zeros;   
ST_ML_se = init_zeros;
ST_ML_hessian = init_zerosH;
ST_ML_cp = init_zeros;   

% Stability (ME)
ST_ME_est = init_zeros; 
ST_ME_cp = init_zeros;  

% Stability + undrominated strategies (ME+MI)
ST_MEI_est = init_zeros; 
ST_MEI_cp = init_zeros;

% ----------------------------%
% Initialize results of tests %
% ----------------------------%

% Truth-telling vs. Stability (Hausman test)
TT_vs_ST_Hausman_test = NaN(PARAM.M,1);
TT_vs_ST_Hausman_rejectTT = NaN(PARAM.M,1);
TT_vs_ST_pvalue = NaN(PARAM.M,1);

% Stability vs. Undominated Strategies
Stability_MEI = NaN(PARAM.M,1);

clear init_zeros init_zerosH;

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                       PREPARE DATA                                              %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

fprintf('\nESTIMATION\n');

% To monitor progress of parfor loops
parfor_progress(PARAM.M);

parfor mm=1:PARAM.M % loop over MC samples

% ------------ %
% PREPARE DATA %
% ------------ %
    
% Retrieve MC sample from simulated data
sim_data_mm = SIM_DATA(SIM_DATA.Sample_id==mm,:);

% School id
sch_id_mm = sim_data_mm.sch_id;

I = PARAM.I;
J = PARAM.J;
A = PARAM.A;

% SCHOOL FIXED EFFECTS %

% School dummy variables
sch_fe_mm = [];
for jj=1:J
   sch_fe_mm(:,jj) = (sch_id_mm == jj);
end

% School scores
sch_score = repmat(MC.School_mscores',[I 1]);

% Own score x school score
sch_stu_score = sim_data_mm.stu_score .* sch_score;

% COVARIATES %

% 2-D array (I*J,Z)
z_all_2D_mm = [sch_fe_mm, sch_stu_score, sim_data_mm.distance];

% # of covariates
nb_z = size(z_all_2D_mm,2);

% List of assigned students
stu_ass_id_mm = sim_data_mm.stu_id(sim_data_mm.sch_assignment==1,:);

% # of assigned students (scalar)
nb_ass_mm = numel(stu_ass_id_mm);

% Logical array to select assigned students
is_assigned_mm = ismember(sim_data_mm.stu_id, stu_ass_id_mm);

% Matrix of dummies for student's asigned school (IA x 5)
stu_ass_assignment_mm = reshape(sim_data_mm.sch_assignment(is_assigned_mm),J,[])';

% Covariates (3-D array (I,J,Z))
z_all_3D_mm = [];
for zz=1:nb_z
    z_all_3D_mm(:,:,zz)= reshape(z_all_2D_mm(:,zz),J,[])';
end
z_ass_3D_mm = z_all_3D_mm(stu_ass_id_mm,:,:);

% Ex post feasible schools (IA x J matrix of dummy variables)
all_sch_feasible_mm = reshape(sim_data_mm.sch_feasible,J,[])';
ass_sch_feasible_mm = reshape(sim_data_mm.sch_feasible(is_assigned_mm),J,[])';

% Students' true preference ranking of schools (I,5)
all_true_rks_mm = reshape(sim_data_mm.true_rk,J,[])';

% Students' submitted ranking of schools (IA,5). N.B.: 0 if not listed
all_sub_rks_mm = reshape(sim_data_mm.choice_rk,J,[])';
ass_sub_rks_mm = reshape(sim_data_mm.choice_rk(is_assigned_mm),J,[])';

% Student is assigned to a school
all_is_assigned_mm = sum(reshape(sim_data_mm.sch_assignment,J,[])'>0,2);

% INSTRUMENTS FOR MOMENT CONDITIONS
all_stu_score = sim_data_mm.stu_score(sim_data_mm.sch_id == 1);
all_distance_s1 = sim_data_mm.distance(sim_data_mm.sch_id == 1);
all_distance_s2 = sim_data_mm.distance(sim_data_mm.sch_id == 2);

% -------------------------------- %
%  EMPIRICAL MOMENT (IN)EQUALITIES %
% -------------------------------- %

% Empirical moments for equalities %
% -------------------------------- %

% J+C-1 empirical moments: J-1 for school assignment + C for covariates other than school fixed effects
all_assignment_mm =  reshape(sim_data_mm.sch_assignment,J,[])';
all_e_eq_mm = permute(sum(z_all_3D_mm.*repmat(all_assignment_mm,1,1,nb_z),2),[1 3 2]);

% Empirical moments for inequalities %
% ---------------------------------- %

% Matrix of all possible pairs of schools out of J schools
pairs = flipud(combnk(1:J,2));

% # of pairs
nb_pairs=size(pairs,1);

% Both schools are ranked
ass_both_ranked_mm = zeros(nb_ass_mm,nb_pairs);
for pp = 1:nb_pairs
    ass_both_ranked_mm(:,pp) = ( ass_sub_rks_mm(:,pairs(pp,1))>0 & ass_sub_rks_mm(:,pairs(pp,2))>0 );    
end    

% Empirical lower bounds : (N,K) matrix with dummy variables = 1 for students s.t. S_i and S_j are both ranked and S_i<S_j
all_elbnd_ineq_mm = zeros(I,nb_pairs);
% Empirical upper bounds : (N,K) matrix with dummy variables = 1 for students s.t. S_i and S_j are both ranked and either S_i>S_j or S_i and S_j are not both ranked
all_eubnd_ineq_mm = zeros(I,nb_pairs);
for pp = 1:nb_pairs
    % Lower bound on probability of preferring second school j to first school i = both school are ranked and second school preferred to first (lower rank)
    all_elbnd_ineq_mm(:,pp) = (   all_sub_rks_mm(:,pairs(pp,1)) > 0 & all_sub_rks_mm(:,pairs(pp,2)) > 0 ...
                                & all_sub_rks_mm(:,pairs(pp,1)) > all_sub_rks_mm(:,pairs(pp,2)));    
    % Upper bound on probability of preferring second school j to first school i = second school is preferred to first + all those who did not rank both schools
    all_eubnd_ineq_mm(:,pp) = all_elbnd_ineq_mm(:,pp) + ( all_sub_rks_mm(:,pairs(pp,1))==0 | all_sub_rks_mm(:,pairs(pp,2))==0);    
end

% Instruments for moment inequalities %
% ----------------------------------- %

% N.B.: no instruments for moment equalities
% Instruments for moment inequalities: student score, distance to school 1, distance to school 2
all_z_ineq_mm = [ ones(I,1),  all_stu_score, all_distance_s1, all_distance_s2];
% Expand instrument to size of (elbnd_ineq, eubnd_ineq)
all_z_ineq_mm = repmat(all_z_ineq_mm,[size(all_elbnd_ineq_mm,2) + size(all_eubnd_ineq_mm,2) 1]);
all_z_ineq_mm = reshape(all_z_ineq_mm,I,[]);

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                       ESTIMATION                                                %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% ------------------------------------ %
% Parameters of optimization algorithm %
% ------------------------------------ %

MLE_optfminunc = optimoptions(@fminunc,'display','none','Algorithm','quasi-newton','MaxIter',5e+4,'MaxFunEvals',5e+4,'TolFun',1e-6,'TolX',1e-6);     
MLE_optfminsearch = optimset('display','off','MaxIter',5e+3,'MaxFunEvals',5e+3,'TolFun',1e-6,'TolX',1e-6);     

% --------------------------------------------- %
% Tuning parameters from Andrews and Shi (2013) %
% --------------------------------------------- %

BCS_B_n=sqrt(0.4*log(I)/log(log(I))); % Tuning parameter B_n
BCS_eta = 10^(-6); % Tuning parameter eta
BCS_conf_level=95 + BCS_eta*100; % Confidence level adjusted by eta

% --------------- %
% Starting values %
% --------------- %

% (1) True value of theta
% (2) Vector of ones
% (3) Vector of random values with seed determined by sample no.

rng(mm,'twister')
start_val_mm = [SETTINGS.theta_true, zeros(SETTINGS.nb_param,1), rand(SETTINGS.nb_param,1)];

% If distance normalized to -1, cannot have negative sigma: change starting values 1 and 3
if SETTINGS.is_WTT==1
   start_val_mm(end,[1 3]) = 1;
end

% ---------------------------- %
% TRUE PREFERENCES (BENCHMARK) %
% ---------------------------- %

% Method: MLE (rank-ordered/exploded logit)

if sum(strcmp('True_ML',SETTINGS.estimates)) == 1
     
    % Choice set: all schools (exploded)
    True_ML_choice_set_mm = [];
    block= [];
    for rr = 1:J
        if rr==1
            True_ML_choice_set_mm = (all_true_rks_mm >= 0);
            block = True_ML_choice_set_mm;
        else
            block = block - (all_true_rks_mm == rr-1);
            True_ML_choice_set_mm = vertcat(True_ML_choice_set_mm, block);                    
        end
    end  

    % Choices: most preferred school in exploded choice set
    True_ML_choices_mm = [];
    for jj = 1:J
        if jj == 1
            True_ML_choices_mm = (all_true_rks_mm == jj);
        else 
            True_ML_choices_mm = [True_ML_choices_mm; all_true_rks_mm == jj];                    
        end
    end
        
    % Objective function: log-likelihood    
    True_ML_objfunc = @(x)LL(x, True_ML_choices_mm, True_ML_choice_set_mm, z_all_3D_mm, SETTINGS.is_WTT);    
    
    % Point estimates / standard errors
    True_ML_est_= NaN(SETTINGS.nb_param,3);
    True_ML_se_ = NaN(SETTINGS.nb_param,3);
    True_ML_fval_ = NaN(1,3);
    True_ML_hessian_ = NaN(SETTINGS.nb_param,SETTINGS.nb_param,3);
    for ss=1:3 % loop over starting values
        % Estimates
        [True_ML_est_(:,ss), True_ML_fval_(:,ss), ~, ~, ~, True_ML_hessian_(:,:,ss)] = fminunc(True_ML_objfunc,start_val_mm(:,ss),MLE_optfminunc);
        % Standard errors
        True_ML_se_(:,ss) = sqrt(diag(inv(True_ML_hessian_(:,:,ss))));
    end
    
    % Select estimates that minimize the likelihood function
    [~, True_ML_id] = min(True_ML_fval_);
    True_ML_est(:,mm) = True_ML_est_(:,True_ML_id);
    True_ML_se(:,mm) = True_ML_se_(:,True_ML_id);
    True_ML_hessian(:,:,mm) = True_ML_hessian_(:,:,True_ML_id);

    % Coverage probabilities of 95% confidence intervals
    % Test: is the true value of theta(i) in the 95% confidence interval of theta(i)hat?
    True_ML_cp(:,mm) = ( (SETTINGS.theta_true >= True_ML_est(:,mm)-1.96*True_ML_se(:,mm)) & (SETTINGS.theta_true <= True_ML_est(:,mm)+1.96*True_ML_se(:,mm)) );

end

% ------------------ %
% WEAK TRUTH-TELLING %
% ------------------ %

% Method: MLE (rank-ordered/exploded logit)

if sum(strcmp('TT_ML',SETTINGS.estimates)) == 1

    % Choice set: all schools (exploded)
    TT_ML_choice_set_mm = [];
    block= [];
    for rr = 1:A
        if rr == 1
            TT_ML_choice_set_mm = (all_sub_rks_mm >= 0);
            block = TT_ML_choice_set_mm;
        else
            block = block - (all_sub_rks_mm == rr-1);
            TT_ML_choice_set_mm = vertcat(TT_ML_choice_set_mm, block);                    
        end
    end

    % Choices: top-ranked school within exploded choice set (1/0)
    TT_ML_choices_mm = [];
    for rr=1:A
        if rr==1
            TT_ML_choices_mm = (all_sub_rks_mm == rr);
        else 
            TT_ML_choices_mm = [TT_ML_choices_mm; all_sub_rks_mm == rr];                    
        end
    end       
    
    % Objective function: log-likelihood    
    TT_ML_objfunc = @(x)LL(x, TT_ML_choices_mm, TT_ML_choice_set_mm, z_all_3D_mm, SETTINGS.is_WTT);    
         
    % Point estimates / standard errors
    TT_ML_est_ = NaN(SETTINGS.nb_param,3);
    TT_ML_se_ = NaN(SETTINGS.nb_param,3);
    TT_ML_fval_ = NaN(1,3);
    TT_ML_hessian_ = NaN(SETTINGS.nb_param,SETTINGS.nb_param,3);
    
    for ss=1:3 % loop over starting values
        % Estimates
        [TT_ML_est_(:,ss), TT_ML_fval_(:,ss), ~, ~, ~, TT_ML_hessian_(:,:,ss)] = fminunc(TT_ML_objfunc, start_val_mm(:,ss),MLE_optfminunc);
        % Standard errors
        TT_ML_se_(:,ss)=sqrt(diag(inv(TT_ML_hessian_(:,:,ss))));
    end
    
    % Selecte estimates that minimize the likelihood function
    [~, TT_ML_id] = min(TT_ML_fval_);
    TT_ML_est(:,mm) = TT_ML_est_(:,TT_ML_id);
    TT_ML_se(:,mm) = TT_ML_se_(:,TT_ML_id);
    TT_ML_hessian(:,:,mm) = TT_ML_hessian_(:,:,TT_ML_id); 
    
    % Coverage probabilities of 95% confidence intervals 
    TT_ML_cp(:,mm) = ( (SETTINGS.theta_true >= TT_ML_est(:,mm)-1.96*TT_ML_se(:,mm))  & (SETTINGS.theta_true <= TT_ML_est(:,mm)+1.96*TT_ML_se(:,mm)) );
    
end

% --------------- %
% STABILITY (MLE) %
% --------------- %

% Estimation method: MLE (conditional logit)

if sum(strcmp('ST_ML',SETTINGS.estimates)) == 1 
    
    % Choice set: ex post feasible schools
    ST_ML_choice_set_mm = ass_sch_feasible_mm;

    % Choices: assigned school (1/0)
    ST_ML_choices_mm = stu_ass_assignment_mm;

    % Objective function: log-likelihood
    ST_ML_objfunc = @(x)LL(x,ST_ML_choices_mm,ST_ML_choice_set_mm,z_ass_3D_mm,SETTINGS.is_WTT);    
    
    % Point estimates and standard errors
    ST_ML_est_ = NaN(SETTINGS.nb_param,3);
    ST_ML_se_ = NaN(SETTINGS.nb_param,3);
    ST_ML_fval_ = NaN(1,3);
    ST_ML_hessian_ = NaN(SETTINGS.nb_param,SETTINGS.nb_param,3);
     
    for ss=1:3 % loop over starting values
        % Estimates
        [ST_ML_est_(:,ss), ST_ML_fval_(:,ss), ~, ~, ~, ST_ML_hessian_(:,:,ss)] = fminunc(ST_ML_objfunc,start_val_mm(:,ss),MLE_optfminunc);
        % Standard errors
        ST_ML_se_(:,ss)=sqrt(diag(inv(ST_ML_hessian_(:,:,ss))));
    end

    % Select estimates that minimize the likelihood function
    [~, ST_ML_id] = min(ST_ML_fval_);
    ST_ML_est(:,mm) = ST_ML_est_(:,ST_ML_id);
    ST_ML_se(:,mm) = ST_ML_se_(:,ST_ML_id);
    ST_ML_hessian(:,:,mm) = ST_ML_hessian_(:,:,ST_ML_id);     
    
    % Coverage probabilities of 95% confidence intervals
    ST_ML_cp(:,mm) = ( (SETTINGS.theta_true >= ST_ML_est(:,mm)-1.96*ST_ML_se(:,mm))  & (SETTINGS.theta_true <= ST_ML_est(:,mm)+1.96*ST_ML_se(:,mm)) );  

end

% -------------- %
% STABILITY (ME) %
% -------------- %

% Estimation method: moment equalities

if sum(strcmp('ST_ME',SETTINGS.estimates))==1
    
    % Objective function: Cramer-von Mises (CvM)-type criterion
    ST_ME_args = {all_is_assigned_mm, z_all_3D_mm, [], [], all_sch_feasible_mm, all_e_eq_mm, [],[],[],[],[],[], SETTINGS.is_WTT};
    ST_ME_objfunc = @(x)QStat(0,1,0,x,ST_ME_args{:});

    % Point estimates
    ST_ME_est_ = NaN(SETTINGS.nb_param,3);
    ST_ME_se_ = NaN(SETTINGS.nb_param,3);
    ST_ME_fval_ = NaN(1,3);
    ST_ME_hessian_ = NaN(SETTINGS.nb_param,SETTINGS.nb_param,3);
    
    for ss=1:3 % loop over starting values
        % Estimates
        [ST_ME_est_(:,ss), ST_ME_fval_(:,ss), ~, ~, ~, ST_ME_hessian_(:,:,ss)] = fminunc(ST_ME_objfunc,start_val_mm(:,ss),MLE_optfminunc);
         % Standard errors (not efficient)
        ST_ME_se_(:,ss)=sqrt(diag(inv(ST_ME_hessian_(:,:,ss))));
    end
    
    % Select estimates that minimize the likelihood function
    [~, ST_ME_id] = min(ST_ME_fval_);
    ST_ME_est(:,mm) = ST_ME_est_(:,ST_ME_id);
    ST_ME_se(:,mm) = ST_ME_se_(:,ST_ME_id);
    ST_ME_hessian(:,:,mm) = ST_ME_hessian_(:,:,ST_ME_id);  

end

% -------------------------------------------- %
% STABILITY AND UNDOMINATED STRATEGIES (ME+MI) %
% -------------------------------------------- %

% Estimation method: moment equalities + moment inequalities

if sum(strcmp('ST_MEI',SETTINGS.estimates)) == 1
    
    % Objective function: CvM-type criterion
    ST_MEI_args = {all_is_assigned_mm, z_all_3D_mm, nb_pairs, pairs, all_sch_feasible_mm, ...
                   all_e_eq_mm, all_elbnd_ineq_mm, all_eubnd_ineq_mm, [], all_z_ineq_mm,[],[],SETTINGS.is_WTT};
    ST_MEI_objfunc = @(x)QStat(0,1,1,x,ST_MEI_args{:});   
    
    % Point estimates
    ST_MEI_est_ = NaN(SETTINGS.nb_param,3);
    ST_MEI_fval_ = NaN(1,3);
    
    for ss=1:3 % loop over starting values
        % Estimates
        [ST_MEI_est_(:,ss), ST_MEI_fval_(:,ss)] = fminunc(ST_MEI_objfunc, start_val_mm(:,ss),MLE_optfminunc);
    end
    
    % Select estimates that minimize the likelihood function
    [~, ST_MEI_id] = min(ST_MEI_fval_);
    ST_MEI_est(:,mm) = ST_MEI_est_(:,ST_MEI_id);      
    
end

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                           TESTS                   	                          %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% ------------------------------------------ %
% TRUTH-TELLING VS. STABILITY (HAUSMAN TEST) %
% ------------------------------------------ %

% Perform Hausman type to test:
% H0: stability and truth-telling are satisfied against Ha: only stability is satisfied
% Under H0, TT_ML and ST_ML are both consistent but TT_ML is efficient
% Under Ha, only ST_ML is consistent

if sum(strcmp('TT_vs_ST',SETTINGS.tests))==1
   
    % Covariance matrix of estimaties
    V_TT = inv(TT_ML_hessian(:,:,mm));
    V_ST = inv(ST_ML_hessian(:,:,mm));
    % Hausman test statistic
    % N.B.: pinv computes Moore-Penrose pseudoinverse of matrix
    TT_vs_ST_Hausman_test(mm,1) = (ST_ML_est(:,mm) - TT_ML_est(:,mm))'*pinv(V_ST - V_TT)*(ST_ML_est(:,mm) - TT_ML_est(:,mm));
    % d.f.: rank of differences in covariance matrices
    TT_vs_ST_v_rank = rank(V_ST - V_TT);
    % Critical value from Chi-2 distribution with d.f. = v_rank
    TT_vs_ST_crit_Hausman = chi2inv(0.95,TT_vs_ST_v_rank);
    % Is TT rejected in sample?
    TT_vs_ST_Hausman_rejectTT(mm,1) = (abs(TT_vs_ST_Hausman_test(mm,1)) > TT_vs_ST_crit_Hausman);
    TT_vs_ST_pvalue(mm,1) = chi2cdf(abs(TT_vs_ST_Hausman_test(mm,1)),TT_vs_ST_v_rank,'upper');

end

% ------------------------------ %
% STABILITY TEST (NON-EMPTY SET) %
% ------------------------------ %

% Tests that theta_hat from ME+MI is in the CS of theta
% Based on Bugni, Canay & Shi (JoE 2015)

if sum(strcmp('ST_vs_US',SETTINGS.tests))==1

    % TUNING PARAMETERS %
    % ----------------- %

    % (I,1000) Matrix of 1000 random draws for each student (used for bootstrap approximation in BCS - stability test)
    rng(4,'twister')
    draws = randn(I,1000);

    % Kappa: same as in BCS
    Kappa = sqrt(log(I));
    
    % STABILITY TEST FOR ME+MI
    % Starting value: estimates from ME+MI
    [Stability_MEI(mm,1), ~, ~] = Stability_test(1,1,ST_MEI_est(:,mm),all_is_assigned_mm,z_all_3D_mm,nb_pairs,pairs, ...
                                                 all_sch_feasible_mm,all_e_eq_mm,all_elbnd_ineq_mm,all_eubnd_ineq_mm,[],all_z_ineq_mm, ...
                                                 Kappa,BCS_conf_level,BCS_eta,draws,SETTINGS.is_WTT);   
    
end

% -------------------------------------------------------------------------- %
% COVERAGE PROABILITIES OF 95% CONFIDENCE INTERVALS FOR ESTIMATES FROM ME+MI %
% -------------------------------------------------------------------------- %

if sum(strcmp('ST_MEI_CP',SETTINGS.cp)) == 1 
    
    % (I,1000) Matrix of 1000 random draws for each student (used for bootstrap approximation in BCS - stability test)
    rng(4,'twister')
    draws = randn(I,1000);
    % Kappa: same as in BCS
    Kappa = sqrt(log(I));
    % Is each point estimate in 95% confidence interval?
    ST_MEI_cp_mm = NaN(SETTINGS.nb_param,1);
    for pp = 1:SETTINGS.nb_param
        ST_MEI_cp_mm(pp,1) = TestMR('cp',1,1,pp,ST_MEI_est(:,mm),SETTINGS.theta_true(pp),all_is_assigned_mm,z_all_3D_mm,nb_pairs,pairs, ...
                                                           all_sch_feasible_mm,all_e_eq_mm,all_elbnd_ineq_mm,all_eubnd_ineq_mm,[],all_z_ineq_mm, ...
                                                           Kappa,BCS_conf_level,BCS_eta,draws,SETTINGS.is_WTT);
    end
    ST_MEI_cp(:,mm) = ST_MEI_cp_mm;
end

% Update parfor monitor
parfor_progress;

end

% Close parfor monitor
parfor_progress(0)

% STORE RESULTS IN STRUCTURE ARRAYS
if sum(strcmp('True_ML',SETTINGS.estimates)) == 1
    RESULTS.True_ML = struct('est',True_ML_est,'se',True_ML_se,'hessian',True_ML_hessian,'cp',True_ML_cp);
end
if sum(strcmp('TT_ML',SETTINGS.estimates)) == 1
    RESULTS.TT_ML = struct('est',TT_ML_est,'se',TT_ML_se,'hessian',TT_ML_hessian,'cp',TT_ML_cp);
end
if sum(strcmp('ST_ML',SETTINGS.estimates)) == 1
    RESULTS.ST_ML = struct('est',ST_ML_est,'se',ST_ML_se,'hessian',ST_ML_hessian,'cp',ST_ML_cp);
end
if sum(strcmp('ST_ME',SETTINGS.estimates)) == 1
    RESULTS.ST_ME = struct('est',ST_ME_est,'se',ST_ME_se,'hessian',ST_ME_hessian,'cp',ST_ME_cp);
end
if sum(strcmp('ST_MEI',SETTINGS.estimates)) == 1
    RESULTS.ST_MEI = struct('est',ST_MEI_est,'cp',ST_MEI_cp);
end
if sum(strcmp('TT_vs_ST',SETTINGS.tests)) == 1
    RESULTS.TT_vs_ST = struct('TT_vs_ST_Hausman_test',TT_vs_ST_Hausman_test,'TT_vs_ST_Hausman_rejectTT',TT_vs_ST_Hausman_rejectTT,'TT_vs_ST_pvalue',TT_vs_ST_pvalue);
end
if sum(strcmp('ST_vs_US',SETTINGS.tests)) == 1
    RESULTS.ST_vs_US = struct('Stability_MEI',Stability_MEI);
end

clear   True_ML_est True_ML_se True_ML_hessian True_ML_cp ...
        TT_ML_est TT_ML_se TT_ML_hessian TT_ML_cp ...
        ST_ML_est ST_ML_se ST_ML_hessian ST_ML_cp ...
        ST_ME_est ST_ME_se ST_ME_hessian ST_ME_cp ...
        ST_MEI_est ST_MEI_cp ...
        TT_vs_ST_Hausman_test TT_vs_ST_Hausman_rejectTT TT_vs_ST_pvalue...
        Stability_MEI;

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                    SUMMARY OF RESULTS                                           %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

fprintf('\nEXPORT RESULTS\n');

% ----------------------- %
% MEAN OF POINT ESTIMATES %
% ----------------------- %

RESULTS.TT_ML.m_est = mean(RESULTS.TT_ML.est,2);
RESULTS.ST_ML.m_est = mean(RESULTS.ST_ML.est,2);
RESULTS.ST_MEI.m_est = mean(RESULTS.ST_MEI.est,2);

% ----------------------- %
% S.D. OF POINT ESTIMATES %
% ----------------------- %

RESULTS.TT_ML.sd_est = std(RESULTS.TT_ML.est,0,2);
RESULTS.ST_ML.sd_est = std(RESULTS.ST_ML.est,0,2);
RESULTS.ST_MEI.sd_est = std(RESULTS.ST_MEI.est,0,2);

% ---------------------------------- %
% COVERAGE PROBABILITIES OF 95% C.I. %
% ---------------------------------- %

RESULTS.TT_ML.m_cp = mean(RESULTS.TT_ML.cp,2);
RESULTS.ST_ML.m_cp = mean(RESULTS.ST_ML.cp,2);
RESULTS.ST_MEI.m_cp = mean(RESULTS.ST_MEI.cp,2);

% ----- %
% TESTS %
% ----- %

RESULTS.TT_vs_ST.m_TT_vs_ST_Hausman_rejectTT = mean(RESULTS.TT_vs_ST.TT_vs_ST_Hausman_rejectTT,1);
RESULTS.ST_vs_US.m_reject_Stability = 1-mean(RESULTS.ST_vs_US.Stability_MEI,1);