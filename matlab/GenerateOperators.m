
syms te td ta tb tc;
syms tdp;
syms ra rb rc rd;
syms rap rbp rcp rdp;
syms r_alpha r_beta r_gamma r_delta r_epsilon;

assume(ta,'positive');
assume(ta,'real');
assume(tb,'positive');
assume(tb,'real');
assume(tc,'positive');
assume(tc,'real');
assume(td,'positive');
assume(td,'real');
assume(te,'positive');
assume(te,'real');
assume(te > td);
assume(td > ta);
assume(td > tb);
assume(te > tc);
assume(ta >= 0);
assume(tc >= 0);
assume(tb >= 0);

assume(tdp,'positive');
assume(tdp,'real');

assume(ra,'positive');
assume(ra,'real');
assume(rb,'positive');
assume(rb,'real');
assume(rc,'positive');
assume(rc,'real');
assume(rd,'positive');
assume(rd,'real');

assume(rap,'positive');
assume(rap,'real');
assume(rbp,'positive');
assume(rbp,'real');
assume(rcp,'positive');
assume(rcp,'real');
assume(rdp,'positive');
assume(rdp,'real');


% Define the distance constraints
dAB = ra*(td - ta) + rb*(td - tb) ==  rap*(te - ta) + rdp*(te - tdp) + rbp*(tdp - tb);
dAC = ra*(td-ta) + rd*(te - td) + rc*(te - tc) == rap*(te-ta) + rdp*(te-tdp) + rcp*(tdp-tc);
dAE = ra*(td - ta) + rd*(te - td) ==  rap * (te - ta);

dBC = rb*(td-tb) + rd*(te - td) + rc*(te - tc) == rbp*(tdp-tb) + rcp*(tdp-tc);
dBE = rb*(td - tb) + rd*(te - td) == rbp*(tdp - tb) + rdp*(te - tdp);

dCE = rc*(te - tc) == rcp*(tdp - tc) + rdp*(te - tdp);


% Maps distances to distance names
dMap = containers.Map();
dMap('dAB') = dAB;
dMap('dAC') = dAC;
dMap('dAE') = dAE;
dMap('dBC') = dBC;
dMap('dBE') = dBE;
dMap('dCE') = dCE;


% Free parameters
free_params = [rap, rbp, rcp, rdp, tdp];
random_walks = [r_alpha, r_beta, r_gamma, r_delta, r_epsilon];


% All parameters
all_params = [ra, rb, rc, rd, ta, tb, tc, td, te, rap, rbp, rcp, rdp, tdp];
all_params_size = size(all_params);


% Indicates whether each combination of distances is solvable
cMap = containers.Map();

nsamples = 30;
threshold = 0;
sigma = 0.5;
rwalk_width = 0.1;


            
% Sample rates
s_ra = lognrnd(-0.5*sigma*sigma, sigma, [1,nsamples]);
s_rb = lognrnd(-0.5*sigma*sigma, sigma, [1,nsamples]);
s_rc = lognrnd(-0.5*sigma*sigma, sigma, [1,nsamples]);
s_rd = lognrnd(-0.5*sigma*sigma, sigma, [1,nsamples]);


% Sample times
s_ta = rand(1, nsamples);
s_tb = rand(1, nsamples);
s_tc = rand(1, nsamples);
s_td = max(max(s_ta, s_tb), s_tc) + rand(1, nsamples);





nSolvable = 0;
numDist = size(dMap.keys);
for n = 1: 1: numDist(2)
    
    diary searchUnsolvable.txt
    fprintf('-----------------------\n        n = %d\n-----------------------\n', n);
    diary off
    % Get all combinations of n distances
    combinations = nchoosek(dMap.keys, n);
    
    
    % Iterate through all combinations
    ncombo = size(combinations);
    for i = 1: 1: ncombo(1)
        
        
        dSet = sort(combinations(i,:));

        %disp(dSet)
        
        % Get key string for distance combinations
        setKey = join(dSet, '_');
        setKey = setKey{1};
        
        
        % Ensure that all subsets of this set are solvable
        subsets_are_solvable = true;
        if n > 1
            subset_combinations = nchoosek(dSet, n-1);
            n_subset_combo = size(subset_combinations);
           
            for k = 1: 1: n_subset_combo(1)

                dSubSet = sort(subset_combinations(k,:));
                subSetKey = join(dSubSet, '_');
                subSetKey = subSetKey{1};

                if (cMap(subSetKey) == false)
                    subsets_are_solvable = false;
                    break;
                end

            end
        end
        
        
        % Do not search for a solution if at least one of the immediate 
        % subsets do not have a solution
        if subsets_are_solvable == false
            cMap(setKey) = false;
            %fprintf('Can not solve %s because its subsets are unsolvable\n', setKey)
            continue;
        end
        
        
        
        % Get the distance constraint equations
        equations = [];
        for j = 1: 1: n
            key = dSet{j};
            eqn = dMap(key);
            equations = [equations eqn];
        end
        
        
        % Solve the equation
        S = solve(equations, free_params, 'ReturnConditions', 1, 'IgnoreAnalyticConstraints', 1, 'Real', 1);

        
        % Check for solution
        solSize = size(S.rap);
        
        % No solution
        if solSize(1) == 0
            cMap(setKey) = false;
            diary searchUnsolvable.txt
            fprintf('%s is unsolvable\n', setKey);
            diary off
        else
            cMap(setKey) = true;
            diary searchUnsolvable.txt
            fprintf('%s has a solution\n', setKey); 
            diary off
            nSolvable = nSolvable + 1;
        end
    end
   
end



fprintf("There are %d solvable distance subsets\n", nSolvable);

diary searchUnsolvable.txt
fprintf('\n\n###########################\n      SEARCH COMPLETE\n###########################\n\n\n');
diary off



% Find the maximal solvable subsets (ie. the sets of distances which have a
% solution but whose immediate supersets do not have a solution)
maximal_solvable = [];
for n = 1: 1: numDist(2)
    
    
    % Get all combinations of n distances
    combinations = nchoosek(dMap.keys, n);
    
    
    % Iterate through all combinations
    ncombo = size(combinations);
    for i = 1: 1: ncombo(1)
        dSet = sort(combinations(i,:));
        
        
        % Get key string for distance combinations
        setKey = join(dSet, '_');
        setKey = setKey{1};
        
        
        % Check if it is solvable
        if  cMap(setKey) == false
            continue;
        end
        
        
        % SKIP: ignore whether a set is maximally solvable or not
        maximal_solvable = [maximal_solvable, convertCharsToStrings(setKey)];
        continue
        
        
        % Get all the elements which are NOT in this set
        otherElements = setdiff(dMap.keys, dSet);
        other_size = size(otherElements);
        
        % Find all immediate supersets by adding one of these to the
        % current set
        is_maximal_solvable = true;
        for j = 1: 1: other_size(2)
            
            % Get key string for distance combinations
            otherEle = otherElements(1,j);
            otherEle = otherEle{1};
            superset = sort([dSet otherEle]);
            supersetKey = join(superset, '_');
            supersetKey = supersetKey{1};
            
            % Check if superset is solvable
            if cMap(supersetKey) == true
                is_maximal_solvable = false;
                break;
            end
            
            
        end
        
        
        % If maximally solvable add to list
        if is_maximal_solvable == true
            diary searchUnsolvable.txt
            fprintf('%s is maximally solvable\n', setKey);
            maximal_solvable = [maximal_solvable, convertCharsToStrings(setKey)];
            diary off
        end
        
    end
    
    
end


numMaxSolvable = size(maximal_solvable);
fprintf("There are %d maximally solvable sets\n", numMaxSolvable(2));







% Jacobians for the maximally solvable sets
numDets = 0;
for i = 1: 1: numMaxSolvable(2)
    
    
    % Get the equations from the set key and solve the equations
    setKey = maximal_solvable(1, i);
    equation_keys = split(setKey, "_");
    num_equations = size(equation_keys);
    equations = [];
    for e = 1: 1: num_equations(1)
        equations = [equations, dMap(equation_keys(e,1)) ];
    end
    S = solve(equations, free_params, 'ReturnConditions', 1, 'IgnoreAnalyticConstraints', 1, 'Real', 1);
    
    


    % Propose free parameters
    nFreeParams = size(S.parameters);
    if nFreeParams(2) > 0


        free_params_size = size(free_params);
        for p = 1: 1: free_params_size(2) 

            % Check if still a free parameter
            tree_param = char(free_params(p));
            eqn = eval(['S.' tree_param]);



            freeParamsSize = size(S.parameters);
            for fp = 1: 1: freeParamsSize(2)
                freeParam = S.parameters(fp);

                try
                    double(eval(subs(eqn, freeParam, -1000)));

                    % If no error was thrown then we have successfully
                    % substituted tree_param and it is a free
                    % parameter. Set its value to a proposal (based on
                    % the same variable without the 'p' on the end
                    random_walk_var = random_walks(p);
                    proposedValue = extractBetween(tree_param, 1, 2) + random_walk_var;

                    % Iterate through all proposal equations (rap, rbp, etc) and
                    % replace the unknown variable (eg. x) with the
                    % value of the variable before the proposal (eg.
                    % td)
                    for q = 1: 1: free_params_size(2) 
                        param = char(free_params(q));
                        eqn2 = eval(['S.' param]);
                        %disp(param);
                        S = setfield(S, param, subs(eqn2, freeParam, proposedValue));
                    end

                catch ME
                    %fprintf('%s is not a free parameter wrt %s\n', tree_param, freeParam);
                end

            end


        end


    end


    
    
    % Jacobian matrix
    J = [   [diff(S.rap, ra), diff(S.rap, rb), diff(S.rap, rc), diff(S.rap, rd), diff(S.rap, td)], 
            [diff(S.rbp, ra), diff(S.rbp, rb), diff(S.rbp, rc), diff(S.rbp, rd), diff(S.rbp, td)],
            [diff(S.rcp, ra), diff(S.rcp, rb), diff(S.rcp, rc), diff(S.rcp, rd), diff(S.rcp, td)],
            [diff(S.rdp, ra), diff(S.rdp, rb), diff(S.rdp, rc), diff(S.rdp, rd), diff(S.rdp, td)],
            [diff(S.tdp, ra), diff(S.tdp, rb), diff(S.tdp, rc), diff(S.tdp, rd), diff(S.tdp, td)]];
    
    
    Jsize = size(J);
    if Jsize(1) ~= Jsize(2)
        fprintf('%s does not allow unlinked free parameter proposals. Skipping. \n', setKey);
        continue
    end
        
    
    % Get log determinant (keep it expanded so to avoid logging a negative)
    JD = det(J);
    if JD == 0
        continue;
    end
    numDets = numDets + 1;
    
    
    

    outputFnStr = sprintf("package NERVariants;\n\n");
    outputFnStr = strcat(outputFnStr, sprintf("import operators.MetaNEROperator;\n"));
    outputFnStr = strcat(outputFnStr, sprintf("import beast.util.Randomizer;\n\n"));
    
    

    outputFnStr = strcat(outputFnStr, sprintf("public class NEROperator_%s extends MetaNEROperator {\n\n", setKey));

    outputFnStr = strcat(outputFnStr, sprintf("\t@Override\n"));
    outputFnStr = strcat(outputFnStr, sprintf("\tprotected double proposalRates(double rWindowSize, double ta, double tb, double tc, double td, double te, \n"));
    outputFnStr = strcat(outputFnStr, sprintf("\t\t\t\t\t\t\t\t\t\tdouble ra, double rb, double rc, double rd) {\n\n"));
    
    
    outputFnStr = strcat(outputFnStr, sprintf("\t\t// Random proposals\n"));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tdouble r_alpha = 0;\n"));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tdouble r_beta = 0;\n"));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tdouble r_gamma = 0;\n"));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tdouble r_delta = 0;\n"));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tdouble r_epsilon = this.getRandomWalkStepSize(rWindowSize);\n\n\n"));
    
    

    outputFnStr = strcat(outputFnStr, sprintf("\t\t// Propose new rates + times\n"));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tthis.rap = %s;\n", simplify(S.rap)));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tthis.rbp = %s;\n", simplify(S.rbp)));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tthis.rcp = %s;\n", simplify(S.rcp)));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tthis.rdp = %s;\n", simplify(S.rdp)));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tthis.tdp = %s;\n\n\n", simplify(S.tdp)));

    
    outputFnStr = strcat(outputFnStr, sprintf("\t\t// Jacobian determinant\n"));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tdouble JD = %s;\n", JD));
    outputFnStr = strcat(outputFnStr, sprintf("\t\tif (JD <= 0) return Double.NEGATIVE_INFINITY;\n"));
    outputFnStr = strcat(outputFnStr, sprintf("\t\treturn Math.log(JD);\n\n"));
    outputFnStr = strcat(outputFnStr, sprintf("\t}\n\n"));

    outputFnStr = strcat(outputFnStr, sprintf("}\n"));



    % Convert p^e into Math.pow(p, e)
    for p = 1: 1: all_params_size(2) 
        param_name = sprintf("%s", all_params(1,p));
        for e = 1: 1: 9
            outputFnStr = strrep(outputFnStr, sprintf("%s^%d", param_name, e), sprintf("Math.pow(%s, %d)", param_name, e));
        end
    end
        
      
    %outputFnStr
    
    % Print to java
    diary_filename = sprintf('NERVariants/NEROperator_%s.java', setKey);
    diary(diary_filename);
    fprintf(outputFnStr);
    diary off;
        
    fprintf('%s has a determinant of: %s \n', setKey, JD);
    
end

fprintf('%d out of %d have non-zero determinants\n\n', numDets, numMaxSolvable(2));


%save('Oct_15_2019_d15_nocons');



%load('Oct_14_2019');






%layers{1} = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J,' 'K', 'L', 'M', 'N'}

%x = size(S.rap)


