function [assigned, matched] = DA(type, pref_stu, pref_sch, capacities)
% DA   Implements Gale-Shapley Deferred Acceptance algorithm
%
% Inputs:
% - type: 'stu' for student-proposing DA / 'sch' for school-proposing DA  
% - pref_stu: students' ranking of schools ((J,I) matrix sorted by school_id and student_id). Highest rank = most preferred / 0 if unranked.
% - pref_sch: school's ranking of students ((J,I) matrix sorted by school_id and student_id). Highest rank = most preferred / 0 if unranked.
% - capacities: school capacities 
%
% Output:
% - assigned: (I,1) vector of each student's assigned school_id
% - matched: (J,I) matrix of school_id x student_id matches

% By Gabrielle Fack, Julien Grenet and Yinghua He
% 4 Oct. 2018

% Alternative 1 is faster if sample if small / alternative 2 is faster if sample is large
alternative = 2;

% ---------- %
% PARAMETERS %
% ---------- %

% # of students, # of schools
[nsch, nstu] = size(pref_stu);

% ------------------- %
% SCHOOL-PROPOSING DA %
% ------------------- %

if strcmp(type,'sch') == 1
    
    % Logical array to select last admitted student in each school
    offers_selector = zeros(nsch, nstu); 
    for ss=1:nsch
        offers_selector(ss,capacities(ss)) = 1;   
    end
    offers_selector = logical(offers_selector);
    
    stop = 0;
    while stop == 0
        % Sort students by priority in each school
        [priorities_stu, ~] = sort(pref_sch,2,'descend');
        priorities_stu = priorities_stu';
        % Priority of last admitted student
        prior_last_admitted_stu = priorities_stu(offers_selector');
        % Offers are made to all students with higher priority
        offers = (pref_sch >= repmat(prior_last_admitted_stu,[1 nstu]));
        % Add existing match to the matrix
        % Each student chooses the school with the highest utility out of offers
        matching.sch = (pref_stu.*offers == repmat(max(pref_stu.*offers),[nsch 1])).*offers;
        % Declined offers
        declined_offers = offers - matching.sch;
        % Stopping criterion: no new declined offers
        stop = 1 - any(declined_offers(:));
        % Update the preferences of schools that have been declined by students
        pref_sch(logical(declined_offers)) = 0;
    end

end

% -------------------- %
% STUDENT-PROPOSING DA %
% -------------------- %

if strcmp(type,'stu') == 1
    
    % create (J,I) matrix where (j,i) is the id of the student that school j ranks in i-th position (in descending order). col1 = most preferred students
    [~, ranking_sch] = sort(pref_sch,2,'descend');    
    
    % Logical array to select last admitted student in each school
    offers_selector =  zeros(nsch, nstu); 
    for ss=1:nsch
        offers_selector(ss,capacities(ss))=1;   
    end
    offers_selector = logical(offers_selector);

    stop = 0;
    while stop == 0     
        % School to which students propose: most preferred among non-rejected
        proposals = double(pref_stu ~= 0 & pref_stu == repmat(max(pref_stu),[nsch 1]));
        % Schools evaluate proposals
        if alternative == 1
            % Rank of all proposing students (col no = student id)
            stu_ranked = pref_sch.*proposals;
            [evaluate, ~] = sort(stu_ranked,2,'descend');
            evaluate = evaluate';
            % Priority of last admitted student
            cutoff = evaluate(offers_selector');
        elseif alternative == 2          
            for ss=1:nsch % Loop over schools
                % Id of proposing students
                id_proposing = ranking_sch(ss,:).*proposals(ss,ranking_sch(ss,:));
                % remove non-proposing students
                id_proposing(id_proposing == 0) = [];
                % Priority of last admitted student
                % If remaining seats
                if size(id_proposing,2) < capacities(ss)
                   cutoff(ss,1) = 0;                
                % If excess demand (or exactly filled)
                elseif size(id_proposing,2) >= capacities(ss)
                    % iD of last admitted student
                    id_last_admitted = id_proposing(capacities(ss));
                    % Cutoff = priority of last admitted student
                    cutoff(ss,1) = pref_sch(ss,id_last_admitted);
                end
            end
        end        
        % Any student below cutoff is rejected
        rejected_stu = proposals.*(pref_sch<repmat(cutoff,[1 nstu]));
        matching.sch = proposals - rejected_stu;
        % Update the preferences of rejected students
        pref_stu(logical(rejected_stu)) = 0;
        % Stopping criterion: no new declined offers
        stop = 1-any(rejected_stu(:));
    end
end

% ---------------- %
% MATCHING OUTCOME %
% ---------------- %

% matched: (J,I) matrix of school x student matches
matched = matching.sch;

% assigned: (1,I) vector of students' assigned school (0 if unassigned)
assigned = (max(matching.sch.*repmat((1:nsch)',[1 nstu]),[],1))';             