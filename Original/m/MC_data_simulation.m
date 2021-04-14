% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%        MONTE CARLO SIMULATION OF SCHOOL CHOICE DATA UNDER CONSTRAINED DA / DA WITH COST         %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% By Gabrielle Fack, Julien Grenet and Yinghua He
% 4 Oct. 2018

fprintf('\n----------------------------------------------------------------------------------\n');
fprintf(' MONTE CARLO SIMULATION OF SCHOOL CHOICE DATA UNDER CONSTRAINED DA / DA WITH COST ');
fprintf('\n----------------------------------------------------------------------------------\n');
fprintf(' # By Gabrielle Fack, Julien Grenet and Yinghua He\n');

% ----------------------------------------------------------------------------------------------- %
%                                       OPTIONS                                                   %
% ----------------------------------------------------------------------------------------------- %

% Specify colors for figure showing distribution of cutoffs (up to 10 schools)
FIG.Colors=[    0, 0, 1.0; ...          % School 1
                0, 0.5, 0; ...          % School 2
                1.0, 0, 0; ...          % School 3
                0, 0.75, 0.75; ...      % School 4
                0.75, 0, 0.75; ...      % School 5
                0.75, 0.75, 0; ...      % School 6
                0.25, 0.25, 0.25; ...   % School 7
                0.5, 0.5, 0.5; ...      % School 8
                0.75, 0.75, 0.75; ...   % School 9
                1, 0.75, 1];            % School 10

% Set global timer
TIMER.step_total_time=tic;

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                       MODEL SETUP                                               %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% Economy size: I students and J schools
% M Monte Carlo sample replications to simulate unonstrained/constrained choice

% Set timer Part 1
TIMER.step_part1 = tic;

% --------------- %
% MAIN PARAMETERS %
% --------------- %

% Parameters which are specified outside the program

% propor: I/100
% M: # of MC samples
% J: # of schools
% A: maximum size of ROL (A<=J)
% unit_cost: marginal cost of submitting one extra school when |L|>1

% Fixed parameters

% I: # of students
I = propor*100;
% Capacities: school capacities (normalized to total capacity of 95)
Capacities = propor*[10; 10; 5; 10; 30; 30];
% rho: correlation of student priorities across schools
rho = 0.7;
% Z: maximum # of iterations for convergence of distribution of cutoffs under constrained choice
Z = 20;
% Three grops of coefficients: (i) school fixed effects; (ii) coefficient on distance; (iii) own score x school score
% FE: school fixed effects
FE = [10, 10.5, 11, 11.5, 12, 12.5];
% coeff_score: coefficient on own score x school score
coeff_score = 3;
% coeff_dist: coefficient on distance
coeff_dist = -1;

PARAM = struct('propor',propor,'I',I,'J',J,'Capacities',Capacities,'A',A,'M',M,'unit_cost',unit_cost,'Z',Z,'FE',FE,'coeff_score',coeff_score,'coeff_dist',coeff_dist,'rho',rho);

clear Capacities Z FE power coeff_score coeff_dist rho;
        
% Check internal consistency of supplied parameters

% Check: I/J must be an integer. If not, exit program
if mod(I,J) ~=0 && isequal(PARAM.Capacities,(I/J).*ones(J,1))
    disp('Execution stopped: I/J must be an integer.')
    % To exit keyboard, type -dbquit-
    keyboard
end

% Check : A<=J
if (A>J)
    disp('Error: A > J')
    keyboard % To exit keyboard, type -dbquit-
end

% ------------------------------------ %
%  STUDENT SCORES AND PRIORITY INDICES %
% ------------------------------------ %

% M random draws of student scores and priorities between 0 and 1
% Priorities and scores are uniform draws 
% All are correlated with correlation coefficient rho

% MC samples #0;: M samples to compute school quality
MC0.Stu_score = NaN(I,M);
MC0.Priorities = NaN(J,I,M);

% MC samples #1: M new samples to compute equilibrium distribution of cutoff
MC1.Stu_score = NaN(I,M);
MC1.Priorities = NaN(J,I,M);

% MC samples #2: M new samples to generate simulated school choice data
MC2.Stu_score = NaN(I,M);
MC2.Priorities = NaN(J,I,M);

for mm = 1:3*M % Loop over MC set of samples

Cov_Matrix = eye(J+1) + PARAM.rho .* (ones(J+1,J+1) - eye(J+1));
rng(mm) % for reproducibility
U = copularnd('Gaussian',Cov_Matrix,I);
Correlated_draws = unifinv(U,0,1);

% MC Samples #0
if mm <= M
    % Student score is first column
    MC0.Stu_score(:,mm) = Correlated_draws(:,1);
    % Student priorities are columns 2 to J+1 (J in total)
    MC0.Priorities(:,:,mm) = Correlated_draws(:,2:end)';
end

% MC Samples #1
if mm > M && mm <= 2*M
    % Student score is first column
    MC1.Stu_score(:,mm-M) = Correlated_draws(:,1);
    % Student priorities are columns 2 to J+1 (J in total)
    MC1.Priorities(:,:,mm-M) = Correlated_draws(:,2:end)';
end

% MC Samples #2
if mm > 2*M
    % Student score is first column
    MC2.Stu_score(:,mm-2*M) = Correlated_draws(:,1);
    % Student priorities are columns 2 to J+1 (J in total)
    MC2.Priorities(:,:,mm-2*M) = Correlated_draws(:,2:end)';
end

end

clear mm Cov_Matrix U Correlated_draws;

% ------------------ %
% DISTANCE TO SCHOOL %
% ------------------ %

% Geographic coordinates of schools %
% --------------------------------- %

% Schools are located on a circle of radius 1/2 and are equally spaced on that circle

MC.school_x = zeros(J,1);
MC.school_y = zeros(J,1);

for jj=1:J
    MC.school_x(jj) = cosd(90+((jj-1)/J)*360)./2;
    MC.school_y(jj) = sind(90+((jj-1)/J)*360)./2;
end
clear jj;

% Geographic coordinates of students %
% -----------------------((--------- %

% Students are randomly distributed on a disc of radius 1
% New coordinates are randomly drawn for each set of M Monte Carlo samples

coord_x = @(radius,angle,dim1,dim2) (reshape(radius .* cos(angle),dim1,dim2)./2); % X Coordinates
coord_y = @(radius,angle,dim1,dim2) (reshape(radius .* sin(angle),dim1,dim2)./2); % Y Coordinates

% (i) Student coordinates in MC samples #0
rng(12) % to control random number generator
MC0.theta = rand(1,I*M)*(2*pi); % theta: random angle
rng(13) % to control random number generator
MC0.r = sqrt(rand(1,I*M))*2; % r: random radius
MC0.stu_x = coord_x(MC0.r,MC0.theta,I,M); % X Coordinates
MC0.stu_y = coord_y(MC0.r,MC0.theta,I,M); % Y Coordinates

% (ii) Student coordinates in MC samples #1
rng(8) % to control random number generator
MC1.theta = rand(1,I*M)*(2*pi); % theta: random angle
rng(9) % to control random number generator
MC1.r = sqrt(rand(1,I*M))*2; % r: random radius
MC1.stu_x = coord_x(MC1.r,MC1.theta,I,M); % X Coordinates
MC1.stu_y = coord_y(MC1.r,MC1.theta,I,M); % Y Coordinates

% (iii) Student coordinates in MC samples #2
rng(10) % to control random number generator
MC2.theta = rand(1,I*M)*(2*pi); % theta: random angle
rng(11) % to control random number generator
MC2.r = sqrt(rand(1,I*M))*2; % r: random radius
MC2.stu_x = coord_x(MC2.r,MC2.theta,I,M); % X Coordinates
MC2.stu_y = coord_y(MC2.r,MC2.theta,I,M); % Y Coordinates

% Check
MC2.check_x = (MC2.stu_x(:,1));
MC2.check_y = (MC2.stu_y(:,1));

% Circle
MC.circle_x = cos(0:0.01:2*pi);
MC.circle_y = sin(0:0.01:2*pi);

if exist('short_sim','var') == 0 % To perform short version of the simulations (distribution of cutoffs)
    
% Distance to schools %
% ------------------- %

distance_fun = @(x1,x2,y1,y2) (sqrt((x1-x2).^2 + (y1-y2).^2)); 

MC0.distance_school = NaN(I,J,M);
MC1.distance_school = NaN(I,J,M);
MC2.distance_school = NaN(I,J,M);
for ii=1:I % loop over students
    for jj=1:J % loop over schools
        for mm=1:M % loop over MC samples
            % (i) Distance in MC samples #0
            MC0.distance_school(ii,jj,mm) = distance_fun(MC0.stu_x(ii,mm), MC.school_x(jj), MC0.stu_y(ii,mm),  MC.school_y(jj));
            % (ii) Distance in MC samples #1
            MC1.distance_school(ii,jj,mm) = distance_fun(MC1.stu_x(ii,mm), MC.school_x(jj), MC1.stu_y(ii,mm),  MC.school_y(jj));
            % (iii) Distance in MC samples #2
            MC2.distance_school(ii,jj,mm) = distance_fun(MC2.stu_x(ii,mm), MC.school_x(jj), MC2.stu_y(ii,mm),  MC.school_y(jj));           
        end
    end
end

clear ii jj mm;

% -------------- %
% SCHOOL QUALITY %
% -------------- %

% Run unconstrained DA using 100 MC preliminary samples to set the score of each school (Note: can be simplified)

% Random utility model : U_ij = V_ij + E_ij
% Vij : Deterministic component
% E_ij: Random component (Type-I extreme value)

% Deterministic component of utility (without school score)
MC0.V_ij_a = repmat(PARAM.FE,[I,1,M]) + PARAM.coeff_dist .* MC0.distance_school;

% Idiosyncratic component of utility
rng(25)
MC0.E_ij = -log(-log(rand(I,J,M)));

% Replace priorities by ranks (to avoid ties): use priority indices from MC samples #0
[~, Ranks] = sort(MC0.Priorities,2);
[A1, A2, A3]=ndgrid( 1:size(Ranks,1), 1:size(Ranks,2), 1:size(Ranks,3) );
Ranks(sub2ind(size(Ranks), A1, Ranks, A3)) = A2;
MC0.Ranks = Ranks;
clear A1 A2 A3 Ranks;

MC.School_scores=NaN(J,100);

for mm = 1:100 % Loop over 100 MC samples
    for pp = 1:100 % Loop over iterations until convergence of cutoffs is found

        % Part B of deterministic component of utility (school score)
        if pp == 1
            % First iteration: all school scores are set to zero
             School_score_old2 = zeros(1,J);
             School_score_old = zeros(1,J);
        end
        if pp == 2
            % Subsequent iterations: update average school score
            School_score_old2 = zeros(1,J);
            School_score_old = School_score_new;
        end
        if pp > 2
            % Subsequent iterations: update average school score
            School_score_old2 = School_score_old;
            School_score_old = School_score_new;
        end        

        % Deterministic component of utility
        V_ij = MC0.V_ij_a(:,:,mm) + PARAM.coeff_score .*repmat(School_score_old,I,1) .* repmat(MC0.Stu_score(:,mm),1,J);

        % Utility for each school
        U_ij = V_ij + MC0.E_ij(:,:,mm);

        % ROLs: (IxJ) matrix where (i,j) = id of school ranked j-th in student i's preferences (from top = col 1 to bottom = col J)
        [~, Submitted_ROL] = sort(U_ij,2,'descend');

        % Students' preferences over schools: (I,J) matrix where (i,j) is the preference of student i for school j, in ranks (from 1 for least preferred to J for most preferred)
        Stu_rk_pref = zeros(I,J);
        temp = flip(Submitted_ROL,2);
        [A1, A2] = ndgrid(1:I, 1:J);
        A1(temp == 0) = 0;
        A2(temp == 0) = 0;
        Stu_rk_pref(sub2ind(size(temp), nonzeros(A1), nonzeros(temp))) = nonzeros(A2);
       
        % Run student-proposing DA
        [Assigned, ~] = DA('stu', Stu_rk_pref', MC0.Ranks(:,:,mm), PARAM.Capacities); 
        
        % Determine average score of students assigned to each school
        School_score_new = NaN(1,J);

        for jj = 1:J
            School_score_new(1,jj) = mean(MC0.Stu_score(Assigned == jj,mm),1);
        end

        % Test if average school scores have converged to equilibrium values
        if isequal(School_score_new, School_score_old) || isequal(School_score_new, School_score_old2)
            MC.School_scores(:,mm)=School_score_new';
            break
        else
            continue
        end    
    end
    
end

clear jj mm pp School_score_old School_score_old2 School_score_new U_ij V_ij Stu_rk_pref Submitted_ROL A1 A2 temp Assigned;

% Average school score across 100 preliminary sample
MC.School_mscores = nanmean(MC.School_scores,2)';

% -------------------- %
% STUDENT PREFERENCES  %
% -------------------- %

% Random utility model : U_ij = V_ij + E_ij
% Vij : Deterministic component
% E_ij: Random component (Type-I extreme value)

% V_ij: Deterministic component of utility %
% ---------------------------------------- %

% (i) MC samples #1
MC1.V_ij = repmat(PARAM.FE,[I,1,M]) + PARAM.coeff_dist .* MC1.distance_school + PARAM.coeff_score .* repmat(MC.School_mscores,[I,1,M]).* permute(repmat(MC1.Stu_score,[1 1 J]),[1 3 2]);

% (ii) MC samples #2
MC2.V_ij = repmat(PARAM.FE,[I,1,M]) + PARAM.coeff_dist .* MC2.distance_school + PARAM.coeff_score .* repmat(MC.School_mscores,[I,1,M]).* permute(repmat(MC2.Stu_score,[1 1 J]),[1 3 2]);

% E_ij: Idiosyncratic component of utility %
% ---------------------------------------- %
% N.B.: E_ij are random type-I errors

% (i) MC samples #1
rng(6)
MC1.E_ij = -log(-log(rand(I,J,M)));

% (ii) MC samples #2
rng(5)
MC2.E_ij = -log(-log(rand(I,J,M)));

% U_ij: Utility of each school %
% ---------------------------- %

% (i) MC samples #1

MC1.U_ij = MC1.V_ij + MC1.E_ij;

% Alternative structure (for compatibility with parfor)
MC1_U_ij = permute(MC1.U_ij, [3 2 1]);

% Check that all utilities are non-negative
MC1.negative_utilities = MC1.U_ij(MC1.U_ij<0);

% (ii) MC samples #2
MC2.U_ij = MC2.V_ij + MC2.E_ij;

% Alternative structure (for compatibility with parfor)
MC2_U_ij = permute(MC2.U_ij, [3 2 1]);

% Check that all utilities are non-negative
MC2.negative_utilities = MC2.U_ij(MC2.U_ij<0);

% If a utility is found to be negative, exit
if (isempty(MC1.negative_utilities)==0 || isempty(MC2.negative_utilities)==0)
    disp('Error: negative utilities')
    keyboard % To exit keyboard, type -dbquit-
end

disp(' ')
disp('MODEL PARAMETERS')
disp('----------------')
disp(' ')
disp(['# of students (I): ',num2str(I)])
disp(['# of schools (J): ',num2str(J)])
disp(['# of allowed choices (A): ',num2str(A)])
disp(['# of MC samples: ',num2str(M)])
disp(['Correlation between priorities across schools (rho): ',num2str(PARAM.rho)])
disp(['School capacities: ',mat2str(PARAM.Capacities')])
disp(['School scores: ',mat2str(MC.School_mscores,2)])
disp(['School fixed effects: ',mat2str(PARAM.FE)])
disp(['Coefficient on distance: ',num2str(PARAM.coeff_dist)])
disp(['Coefficient on own score x school score: ',num2str(PARAM.coeff_score)])
disp(['Unit cost of listing an extra school: ',num2str(PARAM.unit_cost)])
disp(['-- execution time: ',num2str(toc(TIMER.step_part1)),' seconds'])

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                               UNCONSTRAINED CHOICE SIMULATIONS                                  %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% Set timer Part 2
TIMER.step_part2 = tic;

disp(' ')
disp('UNCONSTRAINED CHOICE SIMULATIONS')
disp('--------------------------------')
disp(' ')

% ------------------------ %
% COMPLETE RANK-ORDER LIST %
% ------------------------- %

disp('1) Create uncontrained choice lists...')

% N.B.: list of school ranks(LSR) is different from rank-orderd list (ROL)
% ROL: (I,J) matrix where (i,j) = id of school ranked j-th in student i's preferences (from top = col 1 to bottom = col J)
% LSR: (I,J) matrix where (i,j) = rank order of school j in student i's preferences (highest rank = most preferred / 0 if unranked)
[~, DA_UNC1.ROL] = sort(MC1.U_ij,2,'descend');
[~, DA_UNC2.ROL] = sort(MC2.U_ij,2,'descend');

% ---------------------------------- %
% RUN DA TO DETERMINE SCHOOL CUTOFFS %
% ---------------------------------- %

disp(' ')
disp('2) Run DA to determine school cutoffs...')

DA_UNC1.Cutoffs = zeros(J,M);
DA_UNC1.Assigned = zeros(M,I);
DA_UNC1.Matched =  zeros(M,J,I);

DA_UNC2.Cutoffs = zeros(J,M);
DA_UNC2.Assigned = zeros(M,I);
DA_UNC2.Matched =  zeros(M,J,I);

% Replace priorities by ranks (to avoid ties)
[~, DA_UNC1.Ranks] = sort(MC1.Priorities,2);
[A1, A2, A3]=ndgrid( 1:size(DA_UNC1.Ranks,1), 1:size(DA_UNC1.Ranks,2), 1:size(DA_UNC1.Ranks,3) );
DA_UNC1.Ranks(sub2ind(size(DA_UNC1.Ranks), A1, DA_UNC1.Ranks, A3)) = A2;
clear A1 A2 A3;

[~, DA_UNC2.Ranks] = sort(MC2.Priorities,2);
[A1, A2, A3]=ndgrid( 1:size(DA_UNC2.Ranks,1), 1:size(DA_UNC2.Ranks,2), 1:size(DA_UNC2.Ranks,3) );
DA_UNC2.Ranks(sub2ind(size(DA_UNC2.Ranks), A1, DA_UNC2.Ranks, A3)) = A2;
clear A1 A2 A3;

% Compute students' and schools' ranked preferences 

% Stu_rk_pref: Students' preferences ranked over schools
% (J,I,M) matrix where, for a given MC sample, (j,i) is the preference of student i for school j, in ranks (from 0 if unranked to J for most preferred)
DA_UNC1.Stu_rk_pref = zeros(I,J,M);
temp = flipdim(DA_UNC1.ROL,2);
[A1, A2, A3]=ndgrid(1:I, 1:J, 1:M);
A1(temp == 0) = 0;
A2(temp == 0) = 0;
A3(temp == 0) = 0;
DA_UNC1.Stu_rk_pref(sub2ind(size(temp), nonzeros(A1), nonzeros(temp), nonzeros(A3))) = nonzeros(A2);
clear temp A1 A2 A3;

DA_UNC2.Stu_rk_pref = zeros(I,J,M);
temp = flipdim(DA_UNC2.ROL,2);
[A1, A2, A3]=ndgrid(1:I, 1:J, 1:M);
A1(temp == 0) = 0;
A2(temp == 0) = 0;
A3(temp == 0) = 0;
DA_UNC2.Stu_rk_pref(sub2ind(size(temp), nonzeros(A1), nonzeros(temp), nonzeros(A3))) = nonzeros(A2);
clear temp A1 A2 A3;

% Run student-proposing DA

DA_UNC1_Assigned=NaN(M,I);
DA_UNC1_Matched=NaN(M,J,I);
DA_UNC1_Cutoffs=NaN(J,M);

DA_UNC2_Assigned=NaN(M,I);
DA_UNC2_Matched=NaN(M,J,I);
DA_UNC2_Cutoffs=NaN(J,M);

parfor mm = 1:M % Loop over MC samples
   
    % Assigned: each student's assigned school / Matched: each school's assigned students
    [Assigned1, Matched1] = DA('stu', DA_UNC1.Stu_rk_pref(:,:,mm)', DA_UNC1.Ranks(:,:,mm), PARAM.Capacities); 
    [Assigned2, Matched2] = DA('stu', DA_UNC2.Stu_rk_pref(:,:,mm)', DA_UNC2.Ranks(:,:,mm), PARAM.Capacities); 
   
    % Find school cutoffs
    Cutoffs1=NaN(J,1);
    Cutoffs2=NaN(J,1);
    for jj = 1:J % loop over schools
        Cutoffs1(jj) = min( MC1.Priorities(jj, Matched1(jj,:) == 1, mm) );
        Cutoffs2(jj) = min( MC2.Priorities(jj, Matched2(jj,:) == 1, mm) );
        % If school capacity is not filled, then cutoff is 0
        if sum(Matched1(jj,:)) < PARAM.Capacities(jj)
            Cutoffs1(jj)=0;
        end
        % If school cutoff is the student with lowest priority, then cutoff is zero
        if Cutoffs1(jj) == min(MC1.Priorities(jj,:,mm))
            Cutoffs1(jj)=0;
        end
         % If school capacity is not filled, then cutoff is 0
        if sum(Matched2(jj,:)) < PARAM.Capacities(jj)
            Cutoffs2(jj)=0;
        end
        % If school cutoff is the student with lowest priority, then cutoff is zero
        if Cutoffs2(jj) == min(MC2.Priorities(jj,:,mm))
            Cutoffs2(jj)=0;
        end
    end
    DA_UNC1_Assigned(mm,:) = Assigned1;
    DA_UNC1_Matched(mm,:,:) = Matched1;
    DA_UNC1_Cutoffs(:,mm) = Cutoffs1;
    DA_UNC2_Assigned(mm,:) = Assigned2;
    DA_UNC2_Matched(mm,:,:) = Matched2;
    DA_UNC2_Cutoffs(:,mm) = Cutoffs2;
end

DA_UNC1.Assigned = DA_UNC1_Assigned;
DA_UNC1.Matched = DA_UNC1_Matched;
DA_UNC1.Cutoffs = DA_UNC1_Cutoffs;
DA_UNC2.Assigned = DA_UNC2_Assigned;
DA_UNC2.Matched = DA_UNC2_Matched;
DA_UNC2.Cutoffs = DA_UNC2_Cutoffs;

clear DA_UNC1_Assigned DA_UNC1_Matched DA_UNC1_Cutoffs DA_UNC2_Assigned DA_UNC2_Matched DA_UNC2_Cutoffs;

% ------------------------------------------------------------ %
% EMPIRICAL DISTRIBUTION OF CUTOFFS UNDER UNCONSTRAINED CHOICE %
% ------------------------------------------------------------ %

% disp(' ')
% disp('3) Empirical distribution of cutoffs under unconstrained choice...')
% 
% % Non parametric distribution of cutoffs
% figure('Name','Cutoffs under Unconstrained Choice','NumberTitle','off')
% xlabel('School Admission Cutoff'); 
% axis([0 1 0 30])
% xi = linspace(0,1,100);
% f=zeros(J,100);
% x=zeros(J,100);
% legend_names=cell(1,J);
% hold all;
% for jj=1:J
%   [f(jj,:),x(jj,:)] = ksdensity(DA_UNC1.Cutoffs(jj,:),xi); 
%    plot(x(jj,:),f(jj,:),'LineWidth',1.2);
%    legend_names{jj} = ['School' num2str(jj)];
% end
% legend(legend_names);
% hold off;
% close('Cutoffs under Unconstrained Choice (non-Parametric)');
% clear jj xi x f legend_names; 
% 
% disp(['-- execution time: ',num2str(toc(TIMER.step_part2)),' seconds'])

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                               CONSTRAINED CHOICE SIMULATIONS                                    %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% Set timer Part 3
TIMER.step_part3 = tic;

disp(' ')
disp('CONSTRAINED CHOICE SIMULATIONS')
disp('------------------------------')
disp(' ')

% ----------------------------------------------------------------------------------------------- %
%                              FIND ALL POTENTIAL CONSTRAINED ROLs                                %
% ----------------------------------------------------------------------------------------------- %

% Assumption : students can submit a maximum of A schools out of J
% K: number of rank-preserving ordering of A schools among J
% K = J!/(A!(J-A)!)  (e.g., J=5 & A=3 -> K=10)

% If application cost is zeo, only consider lists of size A
if PARAM.unit_cost == 0
    % included: Matrix of dummy variables that indicates the K possible ROLs of A schools
    included = (unique(perms( [ones(1,A), zeros(1,J-A) ] ),'rows'))';
    % Number of combinations for each student
    K = nchoosek(J,A);
end

% If application cost is strictly positive, consider all possible combination of schools of size 1 to A
if PARAM.unit_cost > 0
    included = [];
    for aa=1:A
        included_new = (unique(perms( [ones(1,aa), zeros(1,J-aa) ] ),'rows'))';
        included = horzcat(included, included_new);
    end
    % Number of combinations for each student
    K = size(included,2);
end

PARAM.K = K;

% Replace the INCLUDED matrix of 0 and 1s by the corresponding school positions of included schools (shifted to the left)
[row,col,~]=find(included);         
included_position = full(sparse(row,col,row,size(included,1),size(included,2)));
clear row col val included included_new;
[~, rows] = sort(included_position~=0, 'descend');
cols = repmat(1:size(included_position,2),size(included_position,1),1);
ind = sub2ind(size(included_position),rows(:),cols(:));
sol = NaN(size(included_position,1),size(included_position,2));
sol(:) = included_position(ind);
included_position = sol;
% Keep A schools at max
included_position = included_position(1:A,:);
clear rows cols ind sol;

% ------------------------------- %
% COMPLETE ROL (TRUE PREFERENCES) %
% ------------------------------- %

disp('1) True preferences of students over schools')

% True preferences over all J schools

% (i) MC samples #1
[~, DA_CONS1.True_ROL] = sort(MC1.U_ij,2,'descend');

% Alternative structure (for compatibility with parfor)
DA_CONS1_True_ROL = permute(DA_CONS1.True_ROL, [3 2 1]);

% (ii) MC samples #2
[~, DA_CONS2.True_ROL] = sort(MC2.U_ij,2,'descend');

% Alternative structure (for compatibility with parfor)
DA_CONS2_True_ROL = permute(DA_CONS2.True_ROL, [3 2 1]);


% ----------------------------------------------------------------------------------------------- %
%                    FIND BAYES-NASH EQUILIBRIUM OF CONSTRAINED SCHOOL CHOICE GAME                %
% ----------------------------------------------------------------------------------------------- %

disp('Simulate constrained choice lists')

% Outer loop: start with the cutoffs from unconstrained choice, then use cutoffs under constrained choice until a fixed point is found

% Convergence type
Conv_type=-1;

disp(' ')
disp('3) Find Bayes-Nash equilibrium of constrained choice game')

for zz = 1:PARAM.Z % loop over convergence iterations (20 max)
    % While convergence has not been achieved: use MC samples #1
    if Conv_type<=0
    
        disp(' ')
        disp(['Convergence iteration #', num2str(zz)]);

        U_ij = MC1_U_ij;
        Priorities = MC1.Priorities;
        True_ROL = DA_CONS1_True_ROL;
             
        % First iteration: beliefs are based on the distribution of cutoffs under unconstrained choice
        if zz==1
           disp('   First iteration: beliefs based on distribution of cutoffs under unconstrained choice...');
           Beliefs_Cutoffs = DA_UNC1.Cutoffs;
        % Subsequent iterations: use the distribution of cutoffs obtained under the previous constrained choice simulation
        else 
           disp('   Subsequent iteration: beliefs based on distribution of cutoffs under previous iteration...' );
           % Save previous cutoffs
           Beliefs_Cutoffs_old = Beliefs_Cutoffs;
           % Update beliefs about cutoffs
           Beliefs_Cutoffs=DA_CONS1.Cutoffs;
        end
    
    end
    
    % Once convergence has been achieved: use MC samples #2
    if Conv_type>0
        
        DA_CONS1.Convergence.distance_criterion = distance_criterion;
        DA_CONS1.Convergence.same_distribution = same_distribution;
        DA_CONS1.Convergence.pvalue = pvalue;
        DA_CONS1.Convergence.type = Conv_type;
        clear distance_criterion same_distribution pvalue same_mean same_cov same_distribution;
        
        disp(' ')
        disp('Find best response of students in MC samples #2');
        
        U_ij = MC2_U_ij;
        Priorities = MC2.Priorities;
        True_ROL = DA_CONS2_True_ROL;
               
        % Beliefs about cutoffs are those from final iteration
        Beliefs_Cutoffs=DA_CONS1.Cutoffs;        
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % i) EMPIRICAL ADMISSION PROBABILITIES %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Flatten the tail of the Cauchy cdf which is used to "smooth" probabilities of admission
    % Scaling parameter gamma for Cauchy cdf (left hand side)
    gamma_left=0.0025;
    % Scaling parameter for Cauchy cdf (right hand side)
    gamma_right=0.0025; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ii) SUBMITTED ROLs UNDER CONSTRAINED CHOICE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Artificially expand list of constrained choices by adding two columns of zeros (for the DA algorithm to run)
    Submitted_ROL = zeros(I,J,M);
    
    % Expected admission probabilities
    Submitted_probas = zeros(I,J,M);
    
    % Number of combinations yielding maximum utility
    Combinations_max = NaN(I,1,M);

    disp('   ii) Construct constrained choices...')
    TIMER.step_constrained = tic;   
    parfor mm = 1:M % loop over MC samples
        for ii = 1:I % loop over students

            % Student's true ranking of schools
            true_ranking = True_ROL(mm,:,ii)';
            % Student's priorities at each school
            stu_priorities = Priorities(:,ii,mm);
            % Student's utility for each school, starting with 0 utility for outside option 0 (i.e., index 2 corresponds to school 1)
            utilities = [0; U_ij(mm,:,ii)'];
       
            % All possible combinations of A schools among J, that preserve the ordering of schools (=K)
            [row,col,val] = find(included_position);
            true_ranking_expanded = repmat(true_ranking,1,K);
            combinations_possible = full(sparse(row,col,true_ranking_expanded(sub2ind(size(true_ranking_expanded), val, col)), size(true_ranking_expanded,1), K));
            % Restrict to ROLs of at most A schools
            combinations_possible=combinations_possible(1:A,:);

            % Key step : find the unconditional probabilities of admission to each school in the different ROLs
            
            % Priorities minus cutoffs for every school
            cutoffs_temp = Beliefs_Cutoffs;           
            cutoffs_temp(cutoffs_temp==0)=-0.1;
            prio_cut = repmat(stu_priorities,[1 M])- cutoffs_temp;          
            denominator = (prio_cut>=0).*gamma_right + (prio_cut<0).*gamma_left;
            admission_proba=(0.5+(1/pi)*atan(prio_cut./denominator));                       
            combinations_probabilities = NaN(size(combinations_possible));
            nb_combinations = size(combinations_possible,2);

            for cc = 1:nb_combinations
                combination_cc = combinations_possible(:,cc);
                admission_proba_cc = admission_proba(combination_cc(combination_cc>0),:);
                if A>1
                    for pp=2:size(admission_proba_cc,1)
                        admission_proba_cc(pp,:)=( 1 - sum(admission_proba_cc(1:pp-1,:),1) ) .* admission_proba_cc(pp,:);
                    end
                end
                admission_proba_cc = mean(admission_proba_cc,2)';
                if PARAM.unit_cost>0
                    admission_proba_cc = padarray(admission_proba_cc, [0 A-size(admission_proba_cc,2)], 'post');
                end
                combinations_probabilities(:,cc) = admission_proba_cc;
            end
             
            % if new simulation, compute admission probabilities from previous round and use average of old and new admission probabilities
            if zz > 1 && Conv_type <= 0
                 cutoffs_temp_old = Beliefs_Cutoffs_old;
                 cutoffs_temp_old(cutoffs_temp_old == 0) = -0.1;
                 prio_cut_old = repmat(stu_priorities,[1 M]) - cutoffs_temp_old;          
                 denominator_old = (prio_cut_old >= 0).*gamma_right + (prio_cut_old < 0).*gamma_left;
                 admission_proba_old=(0.5+(1/pi)*atan(prio_cut_old./denominator_old));                 
                 combinations_probabilities_old = NaN(size(combinations_possible));               
                
                for cc = 1:nb_combinations
                    combination_cc = combinations_possible(:,cc);
                    admission_proba_old_cc = admission_proba_old(combination_cc(combination_cc>0),:);
                    if A>1
                        for pp=2:size(admission_proba_old_cc,1)
                            admission_proba_old_cc(pp,:)=( 1 - sum(admission_proba_old_cc(1:pp-1,:),1) ) .* admission_proba_old_cc(pp,:);
                        end
                    end
                    admission_proba_old_cc = mean(admission_proba_old_cc,2)';
                    if PARAM.unit_cost>0
                        admission_proba_old_cc = padarray(admission_proba_old_cc, [0 A-size(admission_proba_old_cc,2)], 'post');
                    end
                    combinations_probabilities_old(:,cc) = admission_proba_old_cc;
                end
                 
                % Update admission probabilities from potential ROLs
                combinations_probabilities = (combinations_probabilities + combinations_probabilities_old)/2;
            end

            % Cumulative admission probabilities
            combinations_cumulative = cumsum(combinations_probabilities,1);

            % When the cumulative probability of admission reaches 1, truncate the list after
            to_remove=cumsum(combinations_probabilities >= 1,1)-1 >= 1;
            combinations_truncated = combinations_possible;
            combinations_truncated(to_remove == 1) = 0;

            % Utilities associated with each combination in the (possibly truncated) list
            % N.B. : use zero-indexing (1+prefix), i.e., if school = 2, search for 3rd row in utilities
            combinations_utilities =  utilities(1+(combinations_possible));       
            if A==1
                combinations_utilities=combinations_utilities';
            end
            % Expected utility from each combination
            combinations_expected_utilities = sum(combinations_probabilities .* combinations_utilities,1); 
            % Expected cost from each combination: start paying the cost if |L|>1
            cost = ((sum(combinations_possible~=0,1)-1))*PARAM.unit_cost;
            % Expected utility net of cost
            combinations_expected_utilities = combinations_expected_utilities - cost;
            % Find the combination that yields the highest expected utility
            index_best = find(combinations_expected_utilities == max(combinations_expected_utilities(:)));
            % Number of combinations yielding maximum utility
            Combinations_max(ii,:,mm) = numel(index_best);
            % There might be several combinations that yield the maximum expected utility: if this happens, pick one at random 
            index_single = index_best(randi(numel(index_best)));
            % Retrieve the corresponding combination of A schools as the student's submitted list under constrained choice
            list = zeros(1,J);
            final_prob = zeros(1,J);
            list(:,1:A) = combinations_truncated(:,index_single)';
            Submitted_ROL(ii,:,mm) = list;
            final_prob(:,1:A) = combinations_probabilities(:,index_single)';
            Submitted_probas(ii,:,mm) = final_prob;

        end % end of loop over students

     end % end of loop over samples
    
     disp(['   - execution time: ',num2str(toc(TIMER.step_constrained)),' seconds'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iii) CUTOFFS UNDER CONSTRAINED CHOICE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    disp('   iii) Run DA to compute cutoffs under constrained choices...');
    TIMER.step_constrained_cutoffs=tic;
         
    % Replace priorities by ranks (to avoid ties)
    [~, Ranks] = sort(Priorities,2);
    [B1, B2, B3] = ndgrid( 1:size(Ranks,1), 1:size(Ranks,2), 1:size(Ranks,3) );
    Ranks(sub2ind(size(Ranks), B1, Ranks,B3)) = B2;
    clear B1 B2 B3;
   
    % Compute students' and schools' rank-preferences 
    % Stu_rk_pref: Students' preferences ranked over schools
    Stu_rk_pref = zeros(I,J,M);
    temp = flip(Submitted_ROL,2);
    [A1, A2, A3]=ndgrid(1:I, 1:J, 1:M);
    A1(temp==0)=0;
    A2(temp==0)=0;
    A3(temp==0)=0;
    Stu_rk_pref(sub2ind(size(temp), nonzeros(A1), nonzeros(temp), nonzeros(A3))) = nonzeros(A2);
    clear temp A1 A2 A3;
    
    % Run student-proposing DA
    Constrained_Assigned=NaN(M,I);
    Constrained_Matched=NaN(M,J,I);
    Constrained_Cutoffs=NaN(J,M);

    parfor mm = 1:M % loop over MC samples
        % Assigned: each student's assigned school / Matched: each school's assigned students
        [Assigned, Matched] = DA('stu', Stu_rk_pref(:,:,mm)', Ranks(:,:,mm), PARAM.Capacities); 
        
        % Find school cutoffs
        Cutoffs=NaN(J,1);
        for jj = 1:J % loop over schools
            Cutoffs(jj) = min(Priorities(jj, Matched(jj,:) == 1, mm) );
            % If school capacity is not filled, then cutoff is 0
            if sum(Matched(jj,:)) < PARAM.Capacities(jj)
                Cutoffs(jj)=0;
            end
            % If school cutoff is the student with lowest priority, then cutoff is zero
            if Cutoffs(jj) == min(Priorities(jj,:,mm))
                Cutoffs(jj)=0;
            end
        end
        Constrained_Assigned(mm,:) = Assigned;
        Constrained_Matched(mm,:,:) = Matched;
        Constrained_Cutoffs(:,mm) = Cutoffs;
    end
    
    if Conv_type <= 0
        DA_CONS1.Beliefs.Cutoffs = Beliefs_Cutoffs;
        DA_CONS1.Stu_rk_pref = Stu_rk_pref;
        DA_CONS1.Submitted_ROL=Submitted_ROL;
        DA_CONS1.Combinations_max = Combinations_max;
        DA_CONS1.Assigned =Constrained_Assigned;
        DA_CONS1.Matched = Constrained_Matched;
        DA_CONS1.Cutoffs = Constrained_Cutoffs;       
    end
    if Conv_type > 0
        DA_CONS2.Beliefs.Cutoffs = Beliefs_Cutoffs;
        DA_CONS2.Stu_rk_pref = Stu_rk_pref;
        DA_CONS2.Submitted_ROL=Submitted_ROL;
        DA_CONS2.Combinations_max = Combinations_max;
        DA_CONS2.Assigned =Constrained_Assigned;
        DA_CONS2.Matched = Constrained_Matched;
        DA_CONS2.Cutoffs = Constrained_Cutoffs;       
    end
    
    clear Constrained_Assigned Constrained_Matched Constrained_Cutoffs Combinations_max
 
    disp(['   - execution time: ',num2str(toc(TIMER.step_constrained_cutoffs)),' seconds'])
   
    if Conv_type > 0
       break 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iv) COMPARE REALIZATION OF CUTOFFS WITH BELIEFS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % While convergence is not achieved
    if Conv_type <= 0
        disp('   iv) Compare current iteration cutoffs with previous iteration cutoffs')
       
        % 1/ Check whether convergence in distribution of cutoffs has been reached (same mean, same covariance)
        
        % Compare means of beliefs and of realized cutoffs
        same_mean = isequal( mean(DA_CONS1.Beliefs.Cutoffs,2), mean(DA_CONS1.Cutoffs,2) );
        % Compare covariance of beliefs and of realized cutoffs
        same_cov=isequal(cov(DA_CONS1.Beliefs.Cutoffs'),cov(DA_CONS1.Cutoffs'));

        % 2/ Compare moments of both distributions       

        moments_beliefs = [mean(DA_CONS1.Beliefs.Cutoffs,2); moment(DA_CONS1.Beliefs.Cutoffs,2,2); moment(DA_CONS1.Beliefs.Cutoffs,3,2)];
        moments_realization = [mean(DA_CONS1.Cutoffs,2); moment(DA_CONS1.Cutoffs,2,2); moment(DA_CONS1.Cutoffs,3,2)];
        distance_criterion = sum( abs(moments_realization - moments_beliefs) );
        disp(['   Moment Distance criterion: ',num2str( distance_criterion)]);

        % 3/ Run nearest-neighbor test to evaluate whether the beliefs/realization of distribution of cutoffs are equal
        [same_distribution, pvalue, ~] = Test_eq_distribution(DA_CONS1.Cutoffs',DA_CONS1.Beliefs.Cutoffs',10,1000,0.99);

        disp(['   Test for equality of distributions: ',num2str(same_distribution)]);
        disp(['   P-value for equality of distributions: ',num2str(pvalue)]);

        if (same_mean==1 && same_cov==1 && same_distribution==1)
            disp('   Exact convergence achieved for distribution of cutoffs')
            Conv_type=1;
            continue
        elseif (zz>1 && same_distribution==1)    
            disp('   Same distributions (nearest neighbor comparison)')
            Conv_type=2;
            continue        
        else
            disp('   Convergence not achieved for distribution of cutoffs: continue...')
            Conv_type=0;
            continue 
        end
    end
    
end % end of major loop on cutoff distributions (convergence)  

clear zz Submitted_probas Submitted_ROL included_position...
      moments_beliefs moments_realization DA_CONS_True_ROL MC_U_ij;

% ----------------------------------------------------------------------------------------------- %
%                    EMPIRICAL DISTRIBUTION OF CUTOFFS UNDER CONSTRAINED CHOICE                   %
% ----------------------------------------------------------------------------------------------- %

disp('Empirical distribution of cutoffs under constrained choice')

nb_points = 100;
x=zeros(J,nb_points);
g=zeros(J,nb_points);
xi = linspace(0,1,nb_points);
for jj=1:J
    [g(jj,:),x(jj,:)] = ksdensity(DA_CONS2.Cutoffs(jj,:),xi);
    legend_names{jj} = ['School' num2str(jj)];
end
clear jj;

if A<J
FIG.DA_CONS.Cutoffs = table( xi', g(1,:)', g(2,:)', g(3,:)', g(4,:)', g(5,:)', g(6,:)', 'VariableNames',{'Cutoff' 'School1' 'School2' 'School3' 'School4' 'School5' 'School6'});
end

if A==J
FIG.DA_CONS.Cutoffs = table( xi', g(1,:)', g(2,:)', g(3,:)', g(4,:)', g(5,:)', g(6,:)', 'VariableNames',{'Cutoff' 'School1' 'School2' 'School3' 'School4' 'School5' 'School6'});
end

clear nb_points xi x1 g1 x g xlabel legend_names h_legend;   

disp(['-- execution time: ',num2str(toc(TIMER.step_part3)),' seconds'])

% ----------------------------------------------------------------------------------------------- %
%                                                                                                 %
%                                   CONSTRAINED CHOICE DATASET                                    %
%                                                                                                 %
% ----------------------------------------------------------------------------------------------- %

% Set timer Part 4
TIMER.step_part4 = tic;

disp(' ')
disp('CONSTRAINED CHOICE DATASET')
disp('--------------------------')
disp(' ')

% Create constrained choice data set
% Number of entries = MxIxJ

% 1: sample_id = sample id
% 2: stu_id = student id
% 3: sch_id = school id
% 4: true_rk = true preference for school
% 5: choice_rk = rank of school in ROL
% 6: stu_priority = student priority
% 7: sch_cutoff = school cutoff
% 8: sch_feasible = whether school is feasible ex post
% 9: sch_assignment = assigned school
% 10: nb_comb_max = number of ROLs yielding maximum utility
% 11: distance = distance to school
% 12: nb_sch_ranked = size of ROL
% 13: true_rk_feasible = true preference order among feasible schools
% 14: choice_rk_feasible = rank order of school in ROL among feasible schools 
% 15: stu_assigned = whether student is assigned to a school
% 16: stu_assigned_pref_feas = whether student is assigned to preferred feasible school
% 17: stu_truthful_feas = whether student is truth-telling among feasible schools
% 18: stu_assigned_pref = whether student is assigned to preferred school
% 20: stu_topk = whether student is weakly truth-telling
% 21: stu_top1 = whether student top ranked most preferred school 
% 22: sch_assignment_unconstrained = assigned school under unconstrained DA
% 23: stu_score = student score

disp('1) Build constrained choice dataset')

% sample_id: sample_id
VAR.sample_id = repmat(1:M,I*J,1);
VAR.sample_id = VAR.sample_id(:);

% stu_id: student id
VAR.stu_id = repmat(1:I,J,M);
VAR.stu_id = VAR.stu_id(:);

% sch_id: school id
VAR.sch_id = repmat((1:J)',M*I,1);

% true_rk: student's true ranking of the school
% DA_CONS.True_ROL (I,J,M)->(I,M,J)
VAR.true_rk_long = permute(DA_CONS2.True_ROL, [1 3 2]);
VAR.true_rk_long = reshape(VAR.true_rk_long,M*I,J);
[row,col,val]=find(VAR.true_rk_long);
VAR.true_rk_long = full(sparse(row,val,col));
VAR.true_rk = reshape(VAR.true_rk_long',1,[])';
clear row col val;

% choice_rk: submitted rank of school (0 = not listed / 1 = first choice, 2 = second choice, etc.)
% DA_CONS.Submitted_ROL (I,J,M)->(I,M,J)
VAR.choice_rk_long = permute(DA_CONS2.Submitted_ROL,[1 3 2]);
VAR.choice_rk_long = reshape(VAR.choice_rk_long,M*I,J);
[row,col,val]=find(VAR.choice_rk_long);
VAR.choice_rk_long = full(sparse(row,val,col));
VAR.choice_rk = reshape(VAR.choice_rk_long',1,[])';
clear row col val;

% stu_priority: student's priority at each school (changes across samples)
% DA_CONS.Priorities (J,I,M)->(I,M,J)
VAR.stu_priority_long = permute(MC2.Priorities,[2 3 1]);
VAR.stu_priority_long = reshape(VAR.stu_priority_long,M*I,J);
VAR.stu_priority = reshape(VAR.stu_priority_long',1,[])';
 
% sch_cutoff: school cutoffs
% DA_CONS.Cutoffs (J,M)->(J,M,I)->(I,M,J)
VAR.sch_cutoff_long = repmat(DA_CONS2.Cutoffs,1,1,I);
VAR.sch_cutoff_long = permute(VAR.sch_cutoff_long,[3 2 1]);
VAR.sch_cutoff_long = reshape(VAR.sch_cutoff_long,M*I,J);
VAR.sch_cutoff = reshape(VAR.sch_cutoff_long',1,[])';

% sch_feas: whether school is feasible for the student
VAR.sch_feasible_long = (VAR.stu_priority_long >= VAR.sch_cutoff_long);
VAR.sch_feasible = reshape(VAR.sch_feasible_long',1,[])';

% sch_assignment: whether student is assigned to that school
VAR.sch_assignment_long = DA_CONS2.Assigned';
VAR.sch_assignment_long = VAR.sch_assignment_long(:);
[row,col,val]=find(VAR.sch_assignment_long);
VAR.sch_assignment_long = full(sparse(row,val,col,I*M,J));
VAR.sch_assignment = reshape(VAR.sch_assignment_long',1,[])';
clear row col val;
   
% nb_comb_max: number of combinations yielding maximum expected utility
VAR.nb_comb_max = repmat(DA_CONS2.Combinations_max,1,J,1);
VAR.nb_comb_max = permute(VAR.nb_comb_max, [2 1 3]);
VAR.nb_comb_max = VAR.nb_comb_max(:);

% distance: distance to school
VAR.distance_long = permute(MC2.distance_school, [1 3 2]);
VAR.distance_long = reshape(VAR.distance_long,M*I,J);
VAR.distance = reshape(VAR.distance_long',1,[])';

% nb_ranked: number of ranked schools
VAR.nb_sch_ranked_long = permute(DA_CONS2.Submitted_ROL,[1 3 2]);
VAR.nb_sch_ranked_long = reshape(VAR.nb_sch_ranked_long,M*I,J);
VAR.nb_sch_ranked_long = sum(VAR.nb_sch_ranked_long~=0,2);
VAR.nb_sch_ranked = reshape(repmat(VAR.nb_sch_ranked_long,1,J)',1,[])';

% true_rk_feas: true preference rank of the school among the ex post feasible
VAR.true_rk_feasible_long = VAR.true_rk_long .* VAR.sch_feasible_long;
VAR.true_rk_feasible_long(VAR.true_rk_feasible_long==0)=NaN;
[~, VAR.true_rk_feasible_long] = sort(VAR.true_rk_feasible_long,2);
[A1, A2]=ndgrid( 1:size(VAR.true_rk_feasible_long,1), 1:size(VAR.true_rk_feasible_long,2) );
VAR.true_rk_feasible_long(sub2ind(size(VAR.true_rk_feasible_long), A1, VAR.true_rk_feasible_long)) = A2;
VAR.true_rk_feasible_long = VAR.true_rk_feasible_long .* VAR.sch_feasible_long;
VAR.true_rk_feasible = reshape(VAR.true_rk_feasible_long',1,[])';
clear A1 A2;

% choice_rk_feas: submitted rank among the ex post feasible
VAR.choice_rk_feasible_long = VAR.choice_rk_long .* VAR.sch_feasible_long;
VAR.choice_rk_feasible_long(VAR.choice_rk_feasible_long==0)=NaN;
[~, VAR.choice_rk_feasible_long] = sort(VAR.choice_rk_feasible_long,2);
[A1, A2]=ndgrid( 1:size(VAR.choice_rk_feasible_long,1), 1:size(VAR.choice_rk_feasible_long,2) );
VAR.choice_rk_feasible_long(sub2ind(size(VAR.choice_rk_feasible_long), A1, VAR.choice_rk_feasible_long)) = A2;
VAR.choice_rk_feasible_long = VAR.choice_rk_feasible_long .* (VAR.choice_rk_feasible_long>0) .* (VAR.sch_feasible_long>0);
VAR.choice_rk_feasible = reshape(VAR.choice_rk_feasible_long',1,[])';
clear A1 A2;

% stu_assigned: whether student is assigned to a school
VAR.stu_assigned_long = sum(VAR.sch_assignment_long~=0,2);
VAR.stu_assigned = reshape(repmat(VAR.stu_assigned_long,1,J)',1,[])';

% stu_assigned_pref_feas: whether student is assigned to most preferred school among the ex post feasible
VAR.true_rk1_feasible_long = (VAR.true_rk_feasible_long==1);
VAR.stu_assigned_pref_feas_long = sum((VAR.sch_assignment_long == VAR.true_rk1_feasible_long),2)==J;
VAR.stu_assigned_pref_feas = reshape(repmat(VAR.stu_assigned_pref_feas_long,1,J)',1,[])';

% stu_truthful_feas: whether ex post feasible schools are ranked according to true preferences
VAR.stu_truthful_feas_long = sum( (VAR.true_rk_feasible_long .* (VAR.choice_rk_feasible_long>0)) == VAR.choice_rk_feasible_long, 2) == J;
VAR.stu_truthful_feas = reshape(repmat(VAR.stu_truthful_feas_long,1,J)',1,[])';

% stu_assigned_pref: whether student is assigned to most preferred school
VAR.true_rk1_long = (VAR.true_rk_long==1);
VAR.stu_assigned_pref_long = sum((VAR.sch_assignment_long == VAR.true_rk1_long),2)==J;
VAR.stu_assigned_pref = reshape(repmat(VAR.stu_assigned_pref_long,1,J)',1,[])';

% stu_topk: whether student has ranked top k most preferred schools (weak truth-telling)
VAR.stu_topk_long = sum( (VAR.true_rk_long .* (VAR.choice_rk_long>0)) == VAR.choice_rk_long, 2) == J;
VAR.stu_topk = reshape(repmat(VAR.stu_topk_long,1,J)',1,[])';

% stu_top1: Student has ranked most preferred school in first position
VAR.choice_rk1_long = (VAR.choice_rk_long==1);
VAR.stu_top1_long = sum(VAR.choice_rk1_long == VAR.true_rk1_feasible_long, 2) == J;
VAR.stu_top1 = reshape(repmat(VAR.stu_top1_long,1,J)',1,[])';

% sch_assignment_unconstrained: student's assignment under unconstrained choice
VAR.sch_assignment_unconstrained_long = DA_UNC2.Assigned';
VAR.sch_assignment_unconstrained_long = VAR.sch_assignment_unconstrained_long(:);
[row,col,val]=find(VAR.sch_assignment_unconstrained_long);
VAR.sch_assignment_unconstrained_long = full(sparse(row,val,col,I*M,J));
VAR.sch_assignment_unconstrained = reshape(VAR.sch_assignment_unconstrained_long',1,[])';
clear row col val;

% stu_score: student's score at the school
% (Ix1) -> (J,I) -> (J,I,M) -> (I,M,J)
VAR.stu_score_long = repmat(MC2.Stu_score,[1 1 J]);
VAR.stu_score_long = reshape(VAR.stu_score_long,M*I,J);
VAR.stu_score = reshape(VAR.stu_score_long',1,[])';

clear SIM_DATA;
SIM_DATA = NaN(M*I*J,22);

% Simulated choice data
SIM_DATA =  table( VAR.sample_id, VAR.stu_id, VAR.sch_id, VAR.true_rk, VAR.choice_rk, VAR.stu_priority, ...
                   VAR.sch_cutoff, VAR.sch_feasible, VAR.sch_assignment, VAR.nb_comb_max, VAR.distance, VAR.nb_sch_ranked, ...
                   VAR.true_rk_feasible, VAR.choice_rk_feasible, VAR.stu_assigned, VAR.stu_assigned_pref_feas, VAR.stu_truthful_feas, ...
                   VAR.stu_assigned_pref, VAR.stu_topk, VAR.stu_top1, VAR.sch_assignment_unconstrained, VAR.stu_score, ...
                   'VariableNames',{'Sample_id' 'stu_id' 'sch_id' 'true_rk' 'choice_rk' 'stu_priority', ...
                   'sch_cutoff' 'sch_feasible' 'sch_assignment' 'nb_comb_max' 'distance' 'nb_sch_ranked', ...
                   'true_rk_feasible' 'choice_rk_feasible' 'stu_assigned' 'stu_assigned_pref_feas' 'stu_truthful_feas' ...
                   'stu_assigned_pref' 'stu_topk' 'stu_top1' 'sch_assignment_unconstrained' 'stu_score'});
               
% ---------------------- %
% DESCRIPTIVE STATISTICS %
% ---------------------- %

disp(' ')
disp('2) Descriptive statistics')

% Tabulate number of ranked school
disp('  Number of ranked schools:')
tabulate(VAR.nb_sch_ranked_long)
STATS.nb_ranked_schools = tabulate(VAR.nb_sch_ranked_long);
STATS.mean_nb_ranked = sum(STATS.nb_ranked_schools(:,1).*STATS.nb_ranked_schools(:,2)./sum(STATS.nb_ranked_schools(:,2),1),1);

% Mean fraction of schools being ranked
STATS.size_ROL = mean(mean(reshape(VAR.nb_sch_ranked_long,I,M),1));
STATS.size_ROL_sd = std(mean(reshape(VAR.nb_sch_ranked_long,I,M),1));
disp(['  Size of ROL: ',sprintf('%.5f',STATS.size_ROL),' (',num2str(STATS.size_ROL_sd),')']);

% Fraction of students who ranked the maximum number of schools
STATS.pct_ranked_max = (sum(VAR.nb_sch_ranked_long==A) == numel(VAR.nb_sch_ranked_long));

% Average fraction of assigned students
STATS.pct_assigned = mean(mean(reshape(VAR.stu_assigned_long,I,M),1));
STATS.pct_assigned_sd = std(mean(reshape(VAR.stu_assigned_long,I,M),1));
STATS.ex_assigned= numel(VAR.stu_assigned_long)-sum(sum(VAR.stu_assigned_long==1));
disp(['  Fraction assigned: ',sprintf('%.5f',STATS.pct_assigned),' (',num2str(STATS.pct_assigned_sd),')','  (# exceptions: ',num2str(STATS.ex_assigned),')']);

% Average fraction of students assigned to favorite feasible school
STATS.pct_assigned_pref_feas = mean(mean(reshape(VAR.stu_assigned_pref_feas_long,I,M),1));
STATS.pct_assigned_pref_feas_sd = std(mean(reshape(VAR.stu_assigned_pref_feas_long,I,M),1));
STATS.ex_assigned_pref_feas= numel(VAR.stu_assigned_pref_feas_long)-sum(sum(VAR.stu_assigned_pref_feas_long==1));
disp(['  Fraction assigned to most preferred feasible school: ',sprintf('%.5f',STATS.pct_assigned_pref_feas),'(',sprintf('%.5f',STATS.pct_assigned_pref_feas_sd),')',' (# exceptions: ',num2str(STATS.ex_assigned_pref_feas),')']);

% Average fraction of students truthful among feasible
STATS.pct_truthful_feas = mean(mean(reshape(VAR.stu_truthful_feas_long,I,M),1));
disp(['  Fraction truthful among feasible: ',sprintf('%.5f',STATS.pct_truthful_feas)]);

% Average fraction of students assigned to preferred school
STATS.pct_assigned_pref = mean(mean(reshape(VAR.stu_assigned_pref_long,I,M),1));
disp(['  Fraction assigned to most preferred school: ',sprintf('%.5f',STATS.pct_assigned_pref)]);

% Average fraction of students who ranked most preferred school first
STATS.pct_top1 = mean(mean(reshape(VAR.stu_top1_long,I,M),1));
disp(['  Fraction ranked most preferred school first: ',sprintf('%.5f',STATS.pct_top1)]);

% Average fraction of students who ranked top k schools (weak truth-telling)
STATS.pct_topk = mean(mean(reshape(VAR.stu_topk_long,I,M),1));
STATS.pct_topk_sd = std(mean(reshape(VAR.stu_topk_long,I,M),1));
disp(['  Fraction weakly truth-telling: ',sprintf('%.5f',STATS.pct_topk),' (',sprintf('%.5f',STATS.pct_topk_sd),')']);

% Average fraction of students who are assigned to the same school as under unconstrained DA
SOSM = (VAR.sch_assignment_long == VAR.sch_assignment_unconstrained_long);
SOSM = (sum(SOSM,2) == size(SOSM,2));
STATS.pct_SOSM = mean(SOSM);
STATS.ex_SOSM = numel(SOSM)-sum(sum(SOSM==1));

disp(['  Fraction assigned to to the same school as under unconstrained DA: ', sprintf('%.5f',STATS.pct_SOSM), '  (# exceptions: ',num2str(STATS.ex_SOSM),')']);

% Compare utility of assigned school to utility of school under unconstrained DA
VAR.stu_utility_long = permute(MC2.U_ij,[1 3 2]);
VAR.stu_utility_long = reshape(VAR.stu_utility_long,M*I,J);
% Utility of assigned school
VAR.stu_utility_assigned_long = sum(VAR.stu_utility_long .* VAR.sch_assignment_long, 2);
% Utility of school under unconstrained DA
VAR.stu_utility_SOSM_long = sum(VAR.stu_utility_long .* VAR.sch_assignment_unconstrained_long, 2);
% Compare utilities: assigned school is less preferred than school under unconstrained DA = -1, same = 0, better = +1
VAR.diff_utility = VAR.stu_utility_assigned_long - VAR.stu_utility_SOSM_long;
VAR.diff_utility(VAR.diff_utility>0)=1;
VAR.diff_utility(VAR.diff_utility<0)=-1;

% Tabulate differences in utilities
disp('  Utility of assigned school is better, same or worse than utility of school under unconstrained DA:')
tabulate(VAR.diff_utility)

disp(' ')
disp(['-- execution time: ',num2str(toc(TIMER.step_part4)),' seconds'])
          
disp(' ')
disp(['TOTAL EXECUTION TIME: ',num2str(toc(TIMER.step_total_time)),' seconds'])

end