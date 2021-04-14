% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                   GOODNESS OF FIT STATISTICS - SOUTERN DISTRICT OF PARIS                        %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% By Gabrielle Fack, Julien Grenet and Yinghua He
% 4 Oct. 2018

fprintf('\n---------------------------------------------------------\n');
fprintf(' GOODNESS OF FIT STATISTICS - SOUTHERN DISTRICT OF PARIS       ');
fprintf('\n---------------------------------------------------------\n');
fprintf(' # By Gabrielle Fack, Julien Grenet and Yinghua He\n');  
    
% Which estimates to use?
GOF.estimates_to_use = [{'TT_ML'},{'ST_ML'},{'ST_MEI'}];

% ESTIMATES
GOF.ESTIMATES.TT_ML = [0;ESTIMATES.TT_ML.est];
GOF.ESTIMATES.ST_ML = [0;ESTIMATES.ST_ML.est];
GOF.ESTIMATES.ST_MEI = [0;ESTIMATES.ST_MEI.est];

% PRIORITIES
GOF.priorities = reshape(VAR.stu_priority,NB.sch,[]);
% Replace priorities by ranks (to avoid ties)
[~, GOF.priorities_ranks] = sort(GOF.priorities,2);
[A1, A2]=ndgrid( 1:size(GOF.priorities_ranks,1), 1:size(GOF.priorities_ranks,2));
GOF.priorities_ranks(sub2ind(size(GOF.priorities_ranks), A1, GOF.priorities_ranks)) = A2;
clear A1 A2;

% STUDENTS' OBSERVED ASSIGNMENT
GOF.obs_assignment = sum(reshape(VAR.sch_assignment,NB.sch,[])'.*repmat((1:NB.sch), [NB.stu 1]),2);

% ------------------------------------------------------------------------------- %
% SIMULATE PREFERENCES USING ESTIMATES BASED ON DIFFERENT IDENTIFYING ASSUMPTIONS %
% ------------------------------------------------------------------------------- %

for type_estimate = GOF.estimates_to_use

disp(['Estimate: ', char(type_estimate)]);

% Select WTT estimate
if strcmp('TT_ML',type_estimate)==1
   estimates_sim = GOF.ESTIMATES.TT_ML;
elseif strcmp('ST_ML',type_estimate)==1   
   estimates_sim = GOF.ESTIMATES.ST_ML;
elseif strcmp('ST_MEI',type_estimate)==1   
   estimates_sim = GOF.ESTIMATES.ST_MEI;
end

% Reconstruct un-normalized theta from estimates
theta_sim = [estimates_sim(1:end-1)./estimates_sim(end); -1./estimates_sim(end)]; 
clear estimates_sim;

% Fixed component of utility
V_ij = repmat(theta_sim(1:NB.sch,1)',[NB.stu,1]);
    for pp=NB.sch+1:NB.Z
        V_ij = V_ij + theta_sim(pp).*Z.ALL.Z_3D(:,:,pp);
    end
clear pp theta_sim;

% Idiosyncratic component of utility
rng(25)
E_ij_good = -log(-log(rand(NB.stu,NB.sch,300)));

% U_ij: Utility of each school
U_ij = repmat(V_ij,[1 1 300]) + E_ij_good;
clear E_ij_good;

% ROLs
[~, ROL] = sort(U_ij,2,'descend');

% If TT estimates, keep only top-k choices
% Create selector based on number of submitted choices
nb_choices = sum(reshape(VAR.choice_rk,NB.sch,[])' > 0,2);
choice_selector = (repmat((1:NB.sch), [NB.stu 1]) <= repmat(nb_choices,[1 NB.sch]));
if strcmp('TT_ML',type_estimate)==1
    ROL = ROL .* repmat(choice_selector,[1 1 300]);
end
clear choice_selector nb_choices;

%Student preferences over schools
% (J,I,M) matrix where, for a given MC sample, (j,i) is the preference of student i for school j, in ranks (from 0 if unranked to J for most preferred)
sim_pref = zeros(NB.stu, NB.sch, 300);
temp = flip(ROL,2);
[A1, A2, A3]=ndgrid(1:NB.stu, 1:NB.sch, 1:300);
A1(temp==0)=0;
A2(temp==0)=0;
A3(temp==0)=0;
sim_pref(sub2ind(size(temp), nonzeros(A1), nonzeros(temp), nonzeros(A3))) = nonzeros(A2);
clear A1 A2 A3 temp;

GOF.(char(type_estimate)).Sim_pref = sim_pref;

% ------------------------ %
% RUN STUDENT-PROPOSING DA %
% ------------------------ %

assigned = NaN(NB.stu,300);
cutoffs = NaN(NB.sch,300);

parfor mm = 1:300 % loop over 300 MC samples

    % Students' preference ranks
    sim_pref_mm = sim_pref(:,:,mm);

    % Assigned: each student's assigned school / Matched: each school's assigned students
    [assigned_mm, matched_mm] = DA('stu', sim_pref_mm', GOF.priorities_ranks, SCH.capacities); 

    % Find school cutoffs
    cutoffs_mm=NaN(NB.sch,1); % loop over schools
    for jj = 1:NB.sch 
        cutoffs_mm(jj) = min( GOF.priorities(jj, matched_mm(jj,:) == 1) );
        % If school capacity is not filled, then cutoff is 0
        if sum(matched_mm(jj,:)) < SCH.capacities(jj)
            cutoffs_mm(jj)=0;
        end
        % If the school cutoff is the student with lowest priority, then cutoff is zero
        if cutoffs_mm(jj) == min(GOF.priorities(jj,:))
            cutoffs_mm(jj)=0;
        end
    end

    assigned(:,mm) = assigned_mm;
    cutoffs(:,mm) = cutoffs_mm;
end

clear sim_pref;

% ----------------------------------------------------------------------------------------------- %
%                                  OBSERVED VS. SIMULATED CUTOFFS                                 %
% ----------------------------------------------------------------------------------------------- %

GOF.(char(type_estimate)).Assigned = assigned;
GOF.(char(type_estimate)).Cutoffs = cutoffs;
GOF.(char(type_estimate)).m_cutoffs = mean(cutoffs,2);
GOF.(char(type_estimate)).sd_cutoffs = std(cutoffs,0,2);

clear U_ij cutoffs;

% ----------------------------------------------------------------------------------------------- %
%                            PREDICTED VS. OBSERVED PARTIAL PREFERENCE ORDER                      %
% ----------------------------------------------------------------------------------------------- %

% 1. MEAN PREDICTTED FRACTION OF STUDENTS ASSIGNED TO OBSERVED ASSIGNMENT
GOF.(char(type_estimate)).prob_assignment = mean(mean(assigned == repmat(GOF.obs_assignment,[1 300]),2),1);
% s.d. across samples
GOF.(char(type_estimate)).prob_assignment_sd = std(mean(assigned == repmat(GOF.obs_assignment,[1 300]),1),0,2);

clear assigned;

% Keep only students who ranked at least two schools (1,587)
stu_choices = reshape(VAR.choice_rk,NB.sch,[])';
sel2 = (sum(stu_choices>0,2)>=2); % selector of students with at least two choices
nb2 = sum(sel2); % nb of students with at least two choices
V_ij2 = V_ij(sel2,:); % Fixed component of utility of students with at least two choices
stu_choices2 = stu_choices(sel2,:); % choices of students with at least two choices
exp_V2 = exp(V_ij2);

clear stu_choices sel2 V_ij;

% 2. PREDICTING OBSERVED ORDERING OF TOP TWO CHOICES

V_ij2_choice1 = sum(V_ij2.*(stu_choices2==1),2); % First choice
V_ij2_choice2 = sum(V_ij2.*(stu_choices2==2),2); % Second choice
prob_ordering = exp(V_ij2_choice1)./(exp(V_ij2_choice1)+exp(V_ij2_choice2));
GOF.(char(type_estimate)).pred_ordering_top2_choices = mean(prob_ordering);
clear V_ij2 V_ij2_choice1 V_ij2_choice2 prob_ordering;

% 3. PREDICTING OBSERVED ORDERING OF ALL SUBMITTED CHOICES

for rr=1:NB.rol
    if rr==1
        choices = (stu_choices2 == rr);
        choice_set = (stu_choices2>0);
        block = choice_set;
    else
        choices = vertcat(choices, stu_choices2 == rr); 
        block = block-(stu_choices2==rr-1);      
        choice_set = vertcat(choice_set, block);                    
    end
end
clear rr block;

prob = repmat(exp_V2, [NB.rol 1]).*choices./repmat(sum(repmat(exp_V2, [NB.rol 1]).*choice_set,2),[1 NB.sch]); % choice probabilities

% transform 2D choice matrix into 3D matrix (A choice situations, I applicants, J schools) 
GOF.tsfm_2D_3D = @(M) permute(reshape(permute(M,[2 1]),NB.sch,nb2,NB.rol),[3 2 1]);

% Transform choices from 2D to 3D
choices_3D = GOF.tsfm_2D_3D(choices);

% Transform choice probabilities from 2D to 3D
prob_3D = GOF.tsfm_2D_3D(prob);

% Only keep the probability of the chosen school in each choice situation
prob_chosen_3D = prob_3D.*choices_3D;

% Keep only one probability per choice situation
prob_chosen = sum(prob_chosen_3D,3);

% replace NaNs by 1
prob_chosen(isnan(prob_chosen))=1;

% Compute product of individual probabilities
prob_chosen = prod(prob_chosen,1);
GOF.(char(type_estimate)).pred_ordering_all_choices = mean(prob_chosen);
   
clear type_estimate nb2 stu_choices2 exp_V2 choices choice_set prob prob choices_3D prob_3D prob_chosen_3D prob_chosen;

end

% ----------------------------------------------------------------------------------------------- %
%                                       EXPORT RESULTS                                            %
% ----------------------------------------------------------------------------------------------- %

% WRITE RESULTS IN PRINT-FRIENDLY MODE

% OBSERVED VS. SIMULATED CUTOFFS

table = cell(50,5);
emptyCells = cellfun(@isempty,table);
table(emptyCells) = {''};
clear emptyCells;

table(1,:) = [ {''},{'Observed'},{'TT_ML'},{'ST_ML'},{'ST_MEI'}];
table(3,1) = {'School cutoffs'};
for ss=1:NB.sch
    table(3+2.*ss-1,1) = {sprintf('s%.0f',ss)};                                
end

% Observed cutoffs
for ss=1:NB.sch
    table(3+2.*ss-1,2) = {sprintf('%.3f',SCH.cutoffs(ss))};
end

for type_estimate = GOF.estimates_to_use

    % Cutoffs
    if strcmp('TT_ML',type_estimate)==1
        col = 3;
    end
    if strcmp('ST_ML',type_estimate)==1
        col = 4;
    end
    if strcmp('ST_MEI',type_estimate)==1
        col = 5;
    end    
    
    for ss=1:NB.sch
    	table(3+2.*ss-1,col) = {sprintf('%.3f',GOF.(char(type_estimate)).m_cutoffs(ss))};
    	table(3+2.*ss,col) = {sprintf('(%.3f)',GOF.(char(type_estimate)).sd_cutoffs(ss))};
    end

end
clear col ss type_estimate

RESULTS.TableD3_Cutoffs = table;
clear table;

% SIMULATED VS. OBSERVED ASSIGNMENT/PREFERENCE ORDER

table = cell(50,4);
emptyCells = cellfun(@isempty,table);
table(emptyCells) = {''};
clear emptyCells;

table(1,:) = [ {''},{'TT_ML'},{'ST_ML'},{'ST_MEI'}];
table(3,1) = {'A. Simulated vs. observed assignment (300 simulated samples)'};
table(4,1) = {'Mean predicted fraction of students assigned to observed assignment'};
table(6,1) = {'B. Predicted vs. observed partial preference order'};
table(7,1) = {'Predicted vs. observed partial preference order of top two choices'};
table(8,1) = {'Predicted vs. observed partial preference order of all submitted choices'};
    
for type_estimate = GOF.estimates_to_use

    % Cutoffs
    if strcmp('TT_ML',type_estimate)==1
        col = 2;
    end
    if strcmp('ST_ML',type_estimate)==1
        col = 3;
    end
    if strcmp('ST_MEI',type_estimate)==1
        col = 4;
    end    
    
  	table(4,col) = {sprintf('%.3f',GOF.(char(type_estimate)).prob_assignment)};
    table(5,col) = {sprintf('(%.3f)',GOF.(char(type_estimate)).prob_assignment_sd)};

 	table(7,col) = {sprintf('%.3f',GOF.(char(type_estimate)).pred_ordering_top2_choices)};
    table(8,col) = {sprintf('%.3f',GOF.(char(type_estimate)).pred_ordering_all_choices)};

end

clear col type_estimate

RESULTS.Table_6_GOF_Measures = table;
clear table;