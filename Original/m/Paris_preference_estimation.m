% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%     ESTIMATION OF STUDENT PREFERENCES UNDER CONSTRAINED DA - SOUTERN DISTRICT OF PARIS          %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% By Gabrielle Fack, Julien Grenet and Yinghua He
% 4 Oct. 2018

fprintf('\n-------------------------------------------------------------------------------------\n');
fprintf(' ESTIMATION OF STUDENT PREFERENCES UNDER CONSTRAINED DA - SOUTHERN DISTRICT OF PARIS       ');
fprintf('\n-------------------------------------------------------------------------------------\n');
fprintf(' # By Gabrielle Fack, Julien Grenet and Yinghua He\n');

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                     TUNING PARAMETERS                                           %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% Specify whether coefficients should be normalized as WTT
SETTINGS.is_WTT = 1;

% ------------------------------- %
% Specify which results to output %
% ------------------------------- %

% Estimation %
% ---------- %

%   'TT_ML'         estimates based on weak truth-telling (MLE)
%   'ST_ML'         estimates based on stability (MLE)
%   'ST_ME'         estimates based on stability (Moment equalities)
%   'ST_MEI'        estimates based on stability + undominated strategies (Moment equalities + inequalities)

% N.B.: Settings are set outside of program
% SETTINGS.estimates = {'TT_ML' 'ST_ML' 'ST_ME' 'ST_MEI'};

% Tests %
% ----- %

%   'Test_TT_ST'    Hausman test for Truth-telling vs. Stability
%   'Test_ST_US'    Test for Stability vs. Undominated Strategies

% Settings are set outside of program
% SETTINGS.tests = {'TT_vs_ST' 'ST_vs_US'};

% Marginal confidence intervals %
% ----------------------------- %

%   'MCI_ST_MEI'    Marginal confidence intervals for estimates based on stability + undominated strategies

% Settins are set outside of program
% SETTINGS.mci = {'MCI_ST_MEI'};

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
if  ( sum(strcmp('ST_vs_US',SETTINGS.tests))==1 || sum(strcmp('MCI_ST_MEI',SETTINGS.mci))==1  ) 
    if sum(strcmp('ST_MEI',SETTINGS.estimates))~=1
        SETTINGS.estimates = [SETTINGS.estimates 'ST_MEI'];
    end
end

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                       PREPARE DATA                                              %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

fprintf('\nDATA\n');
fprintf('  - Import data...');
% import school choice data for Southern District of Paris
DATA = importdata(fullfile(dir_data,'district_sud.xlsx'));
% Assign data to individual variables
for nn=1:size(DATA.data,2)
    VAR.(DATA.colheaders{nn})= DATA.data(:,nn);
end
clear nn;
fprintf(' done\n');

fprintf('  - Prepare data...');

% # of allowed choices
NB.rol = 8; 
% # of students
NB.stu = numel(unique(VAR.stu_id));
% # of schools
NB.sch = numel(unique(VAR.sch_id));
% # of observations
NB.obs = size(VAR.stu_id,1);

% SCHOOL CAPACITIES %
SCH.capacities = VAR.sch_capacity(1:NB.sch);

% OBSERVED CUTOFFS %
SCH.cutoffs = VAR.sch_cutoff(1:NB.sch);

% SCHOOL FIXED EFFECTS %

% School dummy variables
STU_SCH.sch_fe = NaN(NB.obs,NB.sch);
LABELS.sch_fe = cell(NB.sch,1);
for jj=1:NB.sch
   STU_SCH.sch_fe(:,jj) = (VAR.sch_id == jj);
   LABELS.sch_fe(jj,:) = {sprintf('school %g',jj)};
end
clear jj;

% INTERACTION TERMS %

% Student dnb x school dnb (French)
VAR.stu_dnbf_norm =  VAR.stu_dnbf_pct - median(VAR.stu_dnbf_pct);
VAR.stu_sch_dnbf = VAR.stu_dnbf_norm.*VAR.sch_dnbf.*100; 
% Student dnb x school dnb (Math)
VAR.stu_dnbm_norm =  VAR.stu_dnbm_pct - median(VAR.stu_dnbm_pct);
VAR.stu_sch_dnbm = VAR.stu_dnbm_norm.*VAR.sch_dnbm.*100; 
% Student SES x school with high SES (ref cat.: low SES students)
VAR.stu_sch_highses = VAR.stu_highses.*VAR.sch_highses;

% ALL COVARIATES %

% 2-D array (I*J,Z)

Z.ALL.Z_2D = [STU_SCH.sch_fe, VAR.closest_sch, VAR.collocated_sch, VAR.stu_sch_dnbf, VAR.stu_sch_dnbm, VAR.stu_sch_highses, VAR.distance];

LABELS.Z =  vertcat(LABELS.sch_fe, ...
                 'closest school', ...
                 'collocated schools', ...
                 'own score x school dnb (French)', ...
                 'own score x school dnb (Math)', ...
                 'high SES x % high SES in school', ...
                 'distance home-school');

% If coefficient on distance normalized to -1, last coefficient is scaling parameter sigma                    
if SETTINGS.is_WTT == 1
    LABELS.Z(end,1) =  {'scale parameter (sigma)'};
end
                  
% # of covariates
NB.Z = size(Z.ALL.Z_2D,2);

% Determine whether each student was assigned to a school (N x 2 matrix)
% list of assigned students
STU.ASS.stu_id = VAR.stu_id(VAR.sch_assignment==1,:);
% Number of assigned students (scalar)
NB.assigned = numel(STU.ASS.stu_id);
% logical array that will selected assigned students from initial data
VAR.is_assigned = ismember(VAR.stu_id,STU.ASS.stu_id);

% Matrix of dummies for student's asigned school (IA x 5)
STU_SCH.ASS.assignment = reshape(VAR.sch_assignment(VAR.is_assigned),NB.sch,[])';

% Covariates (3-D array IxJxZ)
Z.ALL.Z_3D = [];
Z.ASS.Z_3D = [];

for zz=1:NB.Z
    Z.ALL.Z_3D(:,:,zz) = reshape(Z.ALL.Z_2D(:,zz),NB.sch,[])';
end
clear zz;
Z.ASS.Z_3D = Z.ALL.Z_3D(STU.ASS.stu_id,:,:);

% Ex post feasible schools ((I,J) matrix of dummy variables)
STU_SCH.ALL.sch_feasible = reshape(VAR.sch_feasible,NB.sch,[])';
STU_SCH.ASS.sch_feasible = reshape(VAR.sch_feasible(VAR.is_assigned),NB.sch,[])';

% Students' submitted ranking of schools (I,5)
STU_SCH.ALL.sub_rks = reshape(VAR.choice_rk,NB.sch,[])';
STU_SCH.ASS.sub_rks = reshape(VAR.choice_rk(VAR.is_assigned),NB.sch,[])';

% SELECTORS *

% Student is assigned to a school
STU.ALL.is_assigned = sum(reshape(VAR.sch_assignment,NB.sch,[])' > 0,2);

% INSTRUMENTS *

STU.ALL.stu_dnbm_pct = VAR.stu_dnbm_pct(VAR.sch_id == 1);
STU.ALL.stu_dnbf_pct = VAR.stu_dnbf_pct(VAR.sch_id == 1);
STU.ALL.distance_s1 = VAR.distance(VAR.sch_id == 1);
STU.ALL.distance_s2 = VAR.distance(VAR.sch_id == 2);

fprintf(' done\n');

% STATS *

STU_SCH.ALL.sch_infeasible_and_ranked_top1 = sum((STU_SCH.ALL.sch_feasible==0).*(STU_SCH.ALL.sub_rks==1),2);
STU_SCH.ALL.sch_infeasible_and_ranked_top2 = sum((STU_SCH.ALL.sch_feasible==0).*(STU_SCH.ALL.sub_rks==2),2);
STU_SCH.ALL.sch_infeasible_and_ranked = sum((STU_SCH.ALL.sch_feasible==0).*(STU_SCH.ALL.sub_rks>0),2);
STU_SCH.ALL.nb_ranked = sum(STU_SCH.ALL.sub_rks>=1,2);

disp(['  Fraction ranked at least one infeasible school among top 2 choices: ', ...
    sprintf('%.2f',mean(STU_SCH.ALL.sch_infeasible_and_ranked_top1==1 | STU_SCH.ALL.sch_infeasible_and_ranked_top2 ==1))]);
disp(['  Fraction ranked exactly one infeasible school among top 2 choices: ', ...
    sprintf('%.2f',mean((STU_SCH.ALL.sch_infeasible_and_ranked_top1==1 & STU_SCH.ALL.sch_infeasible_and_ranked_top2 ==0) | ...
    (STU_SCH.ALL.sch_infeasible_and_ranked_top1==0 & STU_SCH.ALL.sch_infeasible_and_ranked_top2 ==1)))]);
disp(['  Fraction ranked exactly two infeasible school among top 2 choices: ', ...
    sprintf('%.2f',mean(STU_SCH.ALL.sch_infeasible_and_ranked_top1==1 & STU_SCH.ALL.sch_infeasible_and_ranked_top2 ==1))]);
disp(['  Average fraction of non feasible choices among all ranked choices: ', ...
    sprintf('%.2f',mean(STU_SCH.ALL.sch_infeasible_and_ranked./STU_SCH.ALL.nb_ranked))]);

% ----------------------------------------------------------------------------------------------- %
%                                EMPIRICAL MOMENT (IN)EQUALITIES                                  %
% ----------------------------------------------------------------------------------------------- % 

fprintf('EMPIRICAL MOMENT (IN)EQUALITIES\n');

% Empirical moments for equalities %
% -------------------------------- %
fprintf('  - Moment equalities...');

% J+N_Z-1 empirical moments: J-1 for school assignment + N_Z for covariates other than school f.e.
STU_SCH.ALL.assignment = reshape(VAR.sch_assignment,NB.sch,[])';
MOMENTS.ALL.e_eq = permute(sum(Z.ALL.Z_3D.*repmat(STU_SCH.ALL.assignment,1,1,NB.Z),2),[1 3 2]);

% Empirical moments for inequalities %
% ---------------------------------- %
fprintf('  - Moment inequalities...');

% -------------------------- %
% Matrix of all school pairs %
% -------------------------- %

% K x 2 matrix indicating the school ids of all possible pairs of schools (used in construction of moment inequalities)

% Matrix of all possible pairs of schools out of J schools
MOMENTS.pairs = flipud(combnk(1:NB.sch,2));

% # of pairs
NB.pairs=size(MOMENTS.pairs,1);
                                   
% Both schools are ranked
MOMENTS.ASS.both_ranked = zeros(NB.assigned,NB.pairs);
for pp=1:NB.pairs
    MOMENTS.ASS.both_ranked(:,pp) = ( STU_SCH.ASS.sub_rks(:,MOMENTS.pairs(pp,1))>0 & STU_SCH.ASS.sub_rks(:,MOMENTS.pairs(pp,2))>0 );    
end    

STATS.m_both_ranked = mean(MOMENTS.ASS.both_ranked,1);

% Empirical lower bounds : (N,K) matrix with dummy variables = 1 for students s.t. S_i and S_j are both ranked and S_i<Sçj
MOMENTS.ALL.elbnd_ineq = zeros(NB.stu,NB.pairs);
% Empirical upper bounds : (N,K) matrix with dummy variables = 1 for students s.t. S_i and S_j are both ranked and S_i>S_j or S_i and S_j are not both ranked
MOMENTS.ALL.eubnd_ineq = zeros(NB.stu,NB.pairs);
for pp=1:NB.pairs
    % lower bound on probability of preferring second school j to first school i = both school are ranked and second school preferred to first (lower rank)
    MOMENTS.ALL.elbnd_ineq(:,pp) = (  STU_SCH.ALL.sub_rks(:,MOMENTS.pairs(pp,1)) > 0 & STU_SCH.ALL.sub_rks(:,MOMENTS.pairs(pp,2)) > 0 ...
                                    & STU_SCH.ALL.sub_rks(:,MOMENTS.pairs(pp,1)) > STU_SCH.ALL.sub_rks(:,MOMENTS.pairs(pp,2)));    
    % upper bound on probability of preferring second school j to first school i = second school is preferred to first + all those who did not rank both schools
    MOMENTS.ALL.eubnd_ineq(:,pp) = MOMENTS.ALL.elbnd_ineq(:,pp) + ( STU_SCH.ALL.sub_rks(:,MOMENTS.pairs(pp,1))==0 | STU_SCH.ALL.sub_rks(:,MOMENTS.pairs(pp,2))==0);    
end
clear pp;

% Statistics on empirical bounds for moment inequalities
STATS.m_elbnd_ineq = mean(MOMENTS.ALL.elbnd_ineq,1);
STATS.m_eubnd_ineq = mean(MOMENTS.ALL.eubnd_ineq,1);
STATS.m_diff_ineq = STATS.m_eubnd_ineq - STATS.m_elbnd_ineq;
STATS.empirical_bounds = [STATS.m_elbnd_ineq', STATS.m_eubnd_ineq'];

% Instruments for moment inequalities %
% ----------------------------------- %

% N.B.: no instruments for moment equalities
% Instruments for moment inequalities (all students, whether assigned or unassigned): student score in French, Math, distance to school 1, distance to school 2
MOMENTS.ALL.z_ineq = [ ones(NB.stu,1),  STU.ALL.stu_dnbf_pct, STU.ALL.stu_dnbm_pct, STU.ALL.distance_s1, STU.ALL.distance_s2];
% Expand instrument to size of (elbnd_ineq, eubnd_ineq)
MOMENTS.ALL.z_ineq = repmat(MOMENTS.ALL.z_ineq,[size(MOMENTS.ALL.elbnd_ineq,2) + size(MOMENTS.ALL.eubnd_ineq,2) 1]);
MOMENTS.ALL.z_ineq = reshape(MOMENTS.ALL.z_ineq,NB.stu,[]);
                            
% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                        ESTIMATION                                               %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

fprintf('\nESTIMATION\n');
  
% ------------------------------------ %
% Parameters of optimization algorithm %
% ------------------------------------ %

SETTINGS.MLE.optfminunc = optimoptions(@fminunc, 'display','none','Algorithm','quasi-newton','MaxIter',5e+4,'MaxFunEvals',5e+4,'TolFun',1e-5,'TolX',1e-5);     
SETTINGS.MLE.optfminsearch = optimset('display','off','MaxIter',5e+3,'MaxFunEvals',5e+3,'TolFun',1e-5,'TolX',1e-5);     

% --------------------------------------------- %
% Tuning parameters from Andrews and Shi (2013) %
% --------------------------------------------- %

SETTINGS.BCS.B_n=sqrt(0.4*log(NB.stu)/log(log(NB.stu))); % Tuning parameter B_n
SETTINGS.BCS.eta = 10^(-6); % Tuning parameter eta
SETTINGS.BCS.conf_level=95 + SETTINGS.BCS.eta*100; % Confidence level adjusted by eta

% Initialize structure array where estimates will be stored
ESTIMATES = struct('TT_ML',{[]},'ST_ML',{[]},'ST_ME',{[]},'ST_MEI',{[]});

% ------------ %
% OUTPUT TABLE %
% ------------ %

% Basic structure of table
empty_tab = cell(NB.Z +1,6);
empty_tab(1,:) = [{'variable'},{'estimate'},{'s.e.'},{'t-stat'},{'LB'},{'UB'}];
empty_tab(2:end,1) = LABELS.Z;

% Specific results tables
RESULTS = struct('TT_ML',{empty_tab},'ST_ML',{empty_tab},'ST_ME',{empty_tab},'ST_MEI',{empty_tab});

clear empty_tab;

% ------------------------------%
% Initialize estimation results %
% ------------------------------%

% Number of parameters
NB.param = NB.Z - 1;

% --------------- %
% STARTING VALUES %
% --------------- %

% (1) Vector of zeros
% (2) Vector of ones
% (3) Vector of random values

rng(4)
OPTIM.start_val = [zeros(NB.param,1), ones(NB.param,1), rand(NB.param,1)];

% If distance normalized to -1, cannot have negative sigma: change starting values 1 and 3
if SETTINGS.is_WTT==1
   OPTIM.start_val(end,[1 3]) = 3;
end

% ------------------ %
% WEAK TRUTH-TELLING %
% ------------------ %

% Method: MLE (rank-ordered/exploded logit)

if sum(strcmp('TT_ML',SETTINGS.estimates))==1 
   
    fprintf('  - Truth-telling (Maximum likelihood) ...\n');   

    % Choice set: all schools (exploded)
    TT_ML.choice_set = [];
    block= [];
    for rr = 1:NB.rol
        if rr==1
            TT_ML.choice_set = (STU_SCH.ALL.sub_rks >= 0);
            block = TT_ML.choice_set;
        else
            block = block - (STU_SCH.ALL.sub_rks == rr-1);
            TT_ML.choice_set = vertcat(TT_ML.choice_set, block);                    
        end
    end
    clear rr block;

    % Choices: top-ranked school within exploded choice set (1/0)
    TT_ML.choices = [];
    for rr=1:NB.rol
        if rr==1
            TT_ML.choices = (STU_SCH.ALL.sub_rks == rr);
        else 
            TT_ML.choices = [TT_ML.choices; STU_SCH.ALL.sub_rks == rr];                    
        end
    end
    clear rr;

    % Objective function: log-likelihood
    TT_ML.objfunc = @(x)LL(x,TT_ML.choices,TT_ML.choice_set,Z.ALL.Z_3D,SETTINGS.is_WTT);    
 
    % Point estimates and standard errors
    TT_ML_est_ = NaN(NB.param,3);
    TT_ML_se_ = NaN(NB.param,3);
    TT_ML_fval_ = NaN(1,3);
    TT_ML_hessian_ = NaN(NB.param,NB.param,3);

    parfor ss=1:3 % loop over starting values
        % Estimates
        [TT_ML_est_(:,ss), TT_ML_fval_(:,ss), ~, ~, ~, TT_ML_hessian_(:,:,ss)] = fminunc(TT_ML.objfunc,OPTIM.start_val(:,ss),SETTINGS.MLE.optfminunc);
        % Standard errors
        TT_ML_se_(:,ss)=sqrt(diag(inv(TT_ML_hessian_(:,:,ss))));
    end
    % Back to structure arrays
    TT_ML.est_ = TT_ML_est_;
    TT_ML.fval_ = TT_ML_fval_;
    TT_ML.se_ = TT_ML_se_;
    TT_ML.hessian_ = TT_ML_hessian_;
    clear ss TT_ML_est_ TT_ML_fval_ TT_ML_se_ TT_ML_hessian_;

    % Select estimates that minimize the likelihood function
    [~, TT_ML.id] = min(TT_ML.fval_);
    TT_ML.est = TT_ML.est_(:,TT_ML.id);
    TT_ML.se = TT_ML.se_(:,TT_ML.id);
    TT_ML.hessian = TT_ML.hessian_(:,:,TT_ML.id); 
    
    % Save results
    RESULTS.TT_ML(3:end, 2:4) = [num2cell(TT_ML.est), num2cell(TT_ML.se), num2cell(abs(TT_ML.est./TT_ML.se))];
    RESULTS.TT_ML(3:end, 5:6) = [num2cell(TT_ML.est-1.96.*TT_ML.se), num2cell(TT_ML.est+1.96.*TT_ML.se)];
    
    % if is_WTT = 1, compute non-normalized estimates to perform Hausman test
    if SETTINGS.is_WTT == 1
        TT_ML.objfunc2 = @(x)LL(x,TT_ML.choices,TT_ML.choice_set,Z.ALL.Z_3D,0);  
        [TT_ML.est2, ~, ~, ~, ~, TT_ML.hessian2] = fminunc(TT_ML.objfunc2,OPTIM.start_val(:,TT_ML.id),SETTINGS.MLE.optfminunc);        
    end
 
    ESTIMATES.TT_ML = TT_ML;
    clear TT_ML;
    
    fprintf(' done\n');
    
end

% --------------- %
% STABILITY (MLE) %
% --------------- %

% Estimation method: MLE (conditional logit)

if sum(strcmp('ST_ML',SETTINGS.estimates)) == 1 
    
    fprintf('  - Stability (Maximum likelihood)...');   

    % Choice set: ex post feasible schools
    ST_ML.choice_set = STU_SCH.ASS.sch_feasible;

    % Choices: assigned school (1/0)
    ST_ML.choices = STU_SCH.ASS.assignment;

    % Objective function: log-likelihood
    ST_ML.objfunc = @(x)LL(x,ST_ML.choices,ST_ML.choice_set,Z.ASS.Z_3D,SETTINGS.is_WTT);    
    
    % Point estimates and standard errors
    ST_ML_est_ = NaN(NB.param,3);
    ST_ML_se_ = NaN(NB.param,3);
    ST_ML_fval_ = NaN(1,3);
    ST_ML_hessian_ = NaN(NB.param,NB.param,3);
     
    parfor ss=1:3 % loop over starting values
        % Estimates
        [ST_ML_est_(:,ss), ST_ML_fval_(:,ss), ~, ~, ~, ST_ML_hessian_(:,:,ss)] = fminunc(ST_ML.objfunc,OPTIM.start_val(:,ss),SETTINGS.MLE.optfminunc);
        % Standard errors
        ST_ML_se_(:,ss)=sqrt(diag(inv(ST_ML_hessian_(:,:,ss))));
    end
    
    % Back to structure arrays
    ST_ML.est_ = ST_ML_est_;
    ST_ML.fval_ = ST_ML_fval_;
    ST_ML.se_ = ST_ML_se_;
    ST_ML.hessian_ = ST_ML_hessian_;
    clear ss ST_ML_est_ ST_ML_fval_ ST_ML_se_ ST_ML_hessian_;    

    % Select estimates that minimize the likelihood function
    [~, ST_ML.id] = min(ST_ML.fval_);
    ST_ML.est = ST_ML.est_(:,ST_ML.id);
    ST_ML.se = ST_ML.se_(:,ST_ML.id);
    ST_ML.hessian = ST_ML.hessian_(:,:,ST_ML.id);     
    
    % Save results
    RESULTS.ST_ML(3:end, 2:4) = [num2cell(ST_ML.est), num2cell(ST_ML.se), num2cell(abs(ST_ML.est./ST_ML.se))];
    RESULTS.ST_ML(3:end, 5:6) = [num2cell(ST_ML.est-1.96.*ST_ML.se), num2cell(ST_ML.est+1.96.*ST_ML.se)];
    
    % if is_WTT = 1, compute non-normalized estimates to perform Hausman test
    if SETTINGS.is_WTT == 1
        ST_ML.objfunc2 = @(x)LL(x,ST_ML.choices,ST_ML.choice_set,Z.ASS.Z_3D,0);  
        [ST_ML.est2, ~, ~, ~, ~, ST_ML.hessian2] = fminunc(ST_ML.objfunc2,OPTIM.start_val(:,ST_ML.id),SETTINGS.MLE.optfminunc);        
    end
  
    ESTIMATES.ST_ML = ST_ML;
    clear ST_ML;
    
    fprintf(' done\n');   

end

% -------------- %
% STABILITY (ME) %
% -------------- %

% Estimation method: moment equalities

if sum(strcmp('ST_ME',SETTINGS.estimates))==1
    
    fprintf('  - Stability (Moment equalities)...');   

    % Objective function: Cramer-von Mises (CvM)-type criterion
    ST_ME.args = {STU.ALL.is_assigned, Z.ALL.Z_3D, [], [], STU_SCH.ALL.sch_feasible, MOMENTS.ALL.e_eq, [],[],[],[],[],[], SETTINGS.is_WTT};    
    ST_ME.objfunc = @(x)QStat(0,1,0,x,ST_ME.args{:});

    % Point estimates
    ST_ME_est_ = NaN(NB.param,3);
    ST_ME_se_ = NaN(NB.param,3);
    ST_ME_fval_ = NaN(1,3);
    ST_ME_hessian_ = NaN(NB.param,NB.param,3);
    
    parfor ss=1:3 % loop over starting values
        % Estimates
        [ST_ME_est_(:,ss), ST_ME_fval_(:,ss), ~, ~, ~, ST_ME_hessian_(:,:,ss)] = fminunc(ST_ME.objfunc,OPTIM.start_val(:,ss),SETTINGS.MLE.optfminunc);
         % Standard errors (not efficient)
        ST_ME_se_(:,ss)=sqrt(diag(inv(ST_ME_hessian_(:,:,ss))));
    end

    % Back to structure arrays
    ST_ME.est_ = ST_ME_est_;
    ST_ME.fval_ = ST_ME_fval_;
    ST_ME.se_ = ST_ME_se_;
    ST_ME.hessian_ = ST_ME_hessian_;
    clear ss ST_ME_est_ ST_ME_fval_ ST_ME_se_ ST_ME_hessian_;      
    
    % Keep the estimate that minimizes the likelihood function
    [~, ST_ME.id] = min(ST_ME.fval_);
    ST_ME.est = ST_ME.est_(:,ST_ME.id);
    ST_ME.se = ST_ME.se_(:,ST_ME.id);
    ST_ME.hessian = ST_ME.hessian_(:,:,ST_ME.id); 
    
    % Save results
    RESULTS.ST_ME(3:end, 2) = num2cell(ST_ME.est);
    RESULTS.ST_ME(3:end, 2:4) = [num2cell(ST_ME.est), num2cell(ST_ME.se), num2cell(abs(ST_ME.est./ST_ME.se))];
    RESULTS.ST_ME(3:end, 5:6) = [num2cell(ST_ME.est-1.96.*ST_ME.se), num2cell(ST_ME.est+1.96.*ST_ME.se)];

    ESTIMATES.ST_ME = ST_ME;
    clear ST_ME;
  
    fprintf(' done\n');   

end
    
% -------------------------------------------- %
% STABILITY AND UNDOMINATED STRATEGIES (ME+MI) %
% -------------------------------------------- %

% Estimation method: moment equalities + moment inequalities

if sum(strcmp('ST_MEI',SETTINGS.estimates))==1
    
    fprintf('  - Stability (Moment equalities + inequalities)...');   

    % Objective function: CvM-type criterion
    ST_MEI.args = {STU.ALL.is_assigned, Z.ALL.Z_3D, NB.pairs, MOMENTS.pairs, STU_SCH.ALL.sch_feasible, ...
                   MOMENTS.ALL.e_eq, MOMENTS.ALL.elbnd_ineq, MOMENTS.ALL.eubnd_ineq, [], MOMENTS.ALL.z_ineq,[],[],SETTINGS.is_WTT};
    ST_MEI.objfunc = @(x)QStat(0,1,1,x,ST_MEI.args{:});   
    
    % Point estimates
    ST_MEI_est_ = NaN(NB.param,3);
    ST_MEI_fval_ = NaN(1,3);
    
    parfor ss=1:3 % loop over starting values
        % Estimates
        [ST_MEI_est_(:,ss), ST_MEI_fval_(:,ss)] = fminunc(ST_MEI.objfunc, OPTIM.start_val(:,ss),SETTINGS.MLE.optfminunc);
    end

    % Back to structure arrays
    ST_MEI.est_ = ST_MEI_est_;
    ST_MEI.fval_ = ST_MEI_fval_;
    clear ss ST_MEI_est_ ST_MEI_fval_ ;
    
    % Select estimates that minimize the likelihood function
    [~, ST_MEI.id] = min(ST_MEI.fval_);
    ST_MEI.est = ST_MEI.est_(:,ST_MEI.id);      

    % Save results
    RESULTS.ST_MEI(3:end, 2) = [num2cell(ST_MEI.est)];
  
    ESTIMATES.ST_MEI = ST_MEI;
    clear ST_MEI;

    fprintf(' done\n');
    
end

% Compute marginal confidence intervals for parameters estimated using ME+MI

if sum(strcmp('MCI_ST_MEI',SETTINGS.mci))==1

    fprintf('  - Marginal confidence intervals...\n')
   
    % ----------------- %
    % TUNING PARAMETERS %
    % ----------------- %

    % (I,300) Matrix of 300 random draws for each student (used for bootstrap approximation in BCS - marginal confidence intervals)
    rng(5,'twister')
    draws = randn(NB.stu, 300);   
    
    % Kappa: as in BCS
    Kappa = sqrt(log(NB.stu));
    
    MCI = struct('draws',draws,'Kappa',Kappa);

    % Optimization algorithm options

    % Version 1: larger tolerance to speed-up estimation
    optfminsearch2 = optimset('display','off','MaxIter',1e+3,'MaxFunEvals',1e+3,'TolFun',1e-2,'TolX',1e-2);   

    % Version 2: smaller tolerance if previous tolerance is too high (slower than v1)
    optfminsearch3 = optimset('display','off','MaxIter',1e+3,'MaxFunEvals',1e+3,'TolFun',1e-3,'TolX',1e-3);           
   
    % Starting value for theta_1 bounds: point estimates from ME+MI
    theta_0 = ESTIMATES.ST_MEI.est;
    
    parfor position = 1:NB.param % loop over parameters

        optfminsearch = optfminsearch2;
        if position == 11 || position == 13 || position == 14
            optfminsearch = optfminsearch3; 
        end   

        % Objective function to find lower and upper bound for theta(i)
        % Find the smallest/largest value of theta(i)0 that belongs to the marginal confidence interval of theta(i) using Test MR in Bugni, Canay and Shi (QE 2017)
        find_lb = @(x)TestMR('lb',1,1,position,theta_0,x,STU.ALL.is_assigned,Z.ALL.Z_3D,NB.pairs,MOMENTS.pairs, ...
                             STU_SCH.ALL.sch_feasible,MOMENTS.ALL.e_eq,MOMENTS.ALL.elbnd_ineq,MOMENTS.ALL.eubnd_ineq,[],MOMENTS.ALL.z_ineq, ...
                             Kappa,SETTINGS.BCS.conf_level,SETTINGS.BCS.eta,draws,SETTINGS.is_WTT);                                   
        find_ub = @(x)TestMR('ub',1,1,position,theta_0,x,STU.ALL.is_assigned,Z.ALL.Z_3D,NB.pairs,MOMENTS.pairs, ...
                             STU_SCH.ALL.sch_feasible,MOMENTS.ALL.e_eq,MOMENTS.ALL.elbnd_ineq,MOMENTS.ALL.eubnd_ineq,[],MOMENTS.ALL.z_ineq, ...
                             Kappa, SETTINGS.BCS.conf_level,SETTINGS.BCS.eta,draws,SETTINGS.is_WTT);                                     

        % Search for lower bound of MCI of theta(i)
        [theta_lb(position,1), ~, ~, ~] = fminsearch(find_lb,theta_0(position),optfminsearch);
 
        % Search for upper bound of MCI of theta(i)
        [theta_ub(position,1), ~, ~, ~] = fminsearch(find_ub,theta_0(position),optfminsearch);
    end

    RESULTS.ST_MEI(3:end, 5:6) = [num2cell(theta_lb), num2cell(theta_ub)];

    clear optfminsearch optfminsearch2 optfminsearch3 draws theta_lb theta_ub output_lb output_ub position;
        
    fprintf(' done\n');
    
end

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                          TESTS                    	                          %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

fprintf('\nTESTS\n');

% ------------------------------------------ %
% TRUTH-TELLING VS. STABILITY (HAUSMAN TEST) %
% ------------------------------------------ %

% Perform Hausman type to test:
% H0: stability and truth-telling are satisfied against Ha: only stability is satisfied
% Under H0, TT_ML and ST_ML are both consistent but TT_ML is efficient
% Under Ha, only ST_ML is consistent

if sum(strcmp('TT_vs_ST',SETTINGS.tests))==1

    fprintf('  - Truth-telling vs. Stability (Hausman test)...')
    
    % Covariance matrix of estimaties
    V_TT = inv(ESTIMATES.TT_ML.hessian);
    V_ST = inv(ESTIMATES.ST_ML.hessian);
    % Hausman test statistic
    % N.B.: pinv computes Moore-Penrose pseudoinverse of matrix
    TT_vs_ST.Hausman_test = (ESTIMATES.ST_ML.est - ESTIMATES.TT_ML.est)'*pinv(V_ST - V_TT)*(ESTIMATES.ST_ML.est - ESTIMATES.TT_ML.est);
    % d.f.: rank of differences in covariance matrices
    TT_vs_ST.v_rank = rank(V_ST - V_TT);
    % Critical value from Chi-2 distribution with d.f. = v_rank
    TT_vs_ST.crit_Hausman = chi2inv(0.95,TT_vs_ST.v_rank);
    % Is TT rejected in sample?
    TT_vs_ST.Hausman_rejectTT = (abs(TT_vs_ST.Hausman_test) > TT_vs_ST.crit_Hausman);
    TT_vs_ST.pvalue = chi2cdf(abs(TT_vs_ST.Hausman_test),TT_vs_ST.v_rank,'upper');
    clear V_TT V_ST;
    
    % if is_WTT = 1, perform Hausman test on non-normalized estimates
    % Covariance matrix of estimaties
    if SETTINGS.is_WTT == 1
        V_TT2 = inv(ESTIMATES.TT_ML.hessian2);
        V_ST2 = inv(ESTIMATES.ST_ML.hessian2);    
        TT_vs_ST.Hausman_test2 = (ESTIMATES.ST_ML.est2 - ESTIMATES.TT_ML.est2)'*pinv(V_ST2 - V_TT2)*(ESTIMATES.ST_ML.est2 - ESTIMATES.TT_ML.est2);
        TT_vs_ST.v_rank2 = rank(V_ST2 - V_TT2);
        TT_vs_ST.crit_Hausman2 = chi2inv(0.95,TT_vs_ST.v_rank2);
        TT_vs_ST.Hausman_rejectTT2 = (abs(TT_vs_ST.Hausman_test2) > TT_vs_ST.crit_Hausman2);
        TT_vs_ST.pvalue2 = chi2cdf(TT_vs_ST.Hausman_test2,TT_vs_ST.v_rank2,'upper');
    end
    
    TESTS.TT_vs_ST = TT_vs_ST;
    clear V_TT2 V_ST2 TT_vs_ST;
 
    fprintf('done\n')

end
   
% ------------------------------ %
% STABILITY TEST (NON-EMPTY SET) %
% ------------------------------ %

% Tests that theta_hat from ME+MI is in the CS of theta
% Based on Bugni, Canay & Shi (JoE 2015)

if ( sum(strcmp('ST_vs_US',SETTINGS.tests))==1 )

    fprintf('  - Stability vs. Undominated strategies (Bugni, Canay & Shi, 2015)...')

    % TUNING PARAMETERS %
    % ----------------- %

    % (I,1000) Matrix of 1000 random draws for each student (used for bootstrap approximation in BCS - stability test)
    rng(4,'twister')
    draws = randn(NB.stu,1000);

    % Kappa: same as in CS
    Kappa = sqrt(log(NB.stu));
    
    TESTS.ST_vs_US = struct('draws',draws,'Kappa',Kappa);
    
    % STABILITY TEST FOR ME+MI
    % Starting value: estimates from ME+MI
    [BCS_MEI_Stability_accept, BCS_MEI_Tn, BCS_MEI_crit_val] = Stability_test(1,1,ESTIMATES.ST_MEI.est,STU.ALL.is_assigned,Z.ALL.Z_3D,NB.pairs,MOMENTS.pairs, ...
                                                                       STU_SCH.ALL.sch_feasible,MOMENTS.ALL.e_eq,MOMENTS.ALL.elbnd_ineq,MOMENTS.ALL.eubnd_ineq,[],MOMENTS.ALL.z_ineq, ...
                                                                       Kappa,SETTINGS.BCS.conf_level,SETTINGS.BCS.eta,draws,SETTINGS.is_WTT);   
    TESTS.ST_vs_US = struct('BCS_MEI_Stability_reject',1-BCS_MEI_Stability_accept,'BCS_MEI_Tn',BCS_MEI_Tn,'BCS_MEI_crit_val',BCS_MEI_crit_val);
    clear draws Kappa BCS_MEI_Stability_accept BCS_MEI_Tn BCS_MEI_crit_val;
 
    fprintf('done\n')
    
end

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                    EXPORT RESULTS                    	                          %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %
   
fprintf('\nEXPORT RESULTS\n');

% WRITE RESULTS IN PRINT-FRIENDLY MODE

if sum(strcmp('TT_ML',SETTINGS.estimates))==1 && sum(strcmp('ST_ML',SETTINGS.estimates))==1 && sum(strcmp('ST_MEI',SETTINGS.estimates)) == 1 ...
&& sum(strcmp('TT_vs_ST',SETTINGS.tests))==1 && sum(strcmp('ST_vs_US',SETTINGS.tests))==1 ...
&& sum(strcmp('MCI_ST_MEI',SETTINGS.mci ))==1
 
table = cell(2.*(NB.Z)+20,4);
emptyCells = cellfun(@isempty,table);
table(emptyCells) = {''};
clear emptyCells;

% Estimates
table(1,:) = [{''},{'TT_ML'},{'ST_ML'},{'ST_MEI'}];
for ll=3:NB.Z+1   
    table(2.*(ll-2),1) = LABELS.Z(ll-1,1);
    table(2.*(ll-2),2) = {sprintf('%.2f',cell2mat(RESULTS.TT_ML(ll,2)))};
    table(2.*(ll-2)+1,2) = {sprintf('[%.2f;%.2f]',cell2mat(RESULTS.TT_ML(ll,5)),cell2mat(RESULTS.TT_ML(ll,6)))};        
    table(2.*(ll-2),3) = {sprintf('%.2f',cell2mat(RESULTS.ST_ML(ll,2)))};
    table(2.*(ll-2)+1,3) = {sprintf('[%.2f;%.2f]',cell2mat(RESULTS.ST_ML(ll,5)),cell2mat(RESULTS.ST_ML(ll,6)))};           
    table(2.*(ll-2),4) = {sprintf('%.2f',cell2mat(RESULTS.ST_MEI(ll,2)))};
    table(2.*(ll-2)+1,4) = {sprintf('[%.2f;%.2f]',cell2mat(RESULTS.ST_MEI(ll,5)),cell2mat(RESULTS.ST_MEI(ll,6)))};                  
end
clear ll;

% Tests
table(2.*(NB.Z)+1,1) = {'Truth-Telling vs. Stability'};
table(2.*(NB.Z)+2,1:2) = [{'Hausman test statistic (ST)'},{sprintf('%.1f',TESTS.TT_vs_ST.Hausman_test2)}];
table(2.*(NB.Z)+3,1:2) = [{'Critical value (ST)'},{sprintf('%.1f',TESTS.TT_vs_ST.crit_Hausman2)}];
table(2.*(NB.Z)+4,1:2) = [{'p-value (ST)'},{sprintf('%.3f',TESTS.TT_vs_ST.pvalue2)}];
table(2.*(NB.Z)+5,1) = {'Stability vs. Undominated Strategies'};
table(2.*(NB.Z)+6,1:2) = [{'BCS test statistic (MEI)'},{sprintf('%.1f',TESTS.ST_vs_US.BCS_MEI_Tn)}];
table(2.*(NB.Z)+7,1:2) = [{'BCS critical value (MEI)'},{sprintf('%.1f',TESTS.ST_vs_US.BCS_MEI_crit_val)}];

RESULTS.Table3_Estimates = table;
clear table;

end