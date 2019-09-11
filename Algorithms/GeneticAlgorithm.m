function [solution, fval, fvals, runtime] = GeneticAlgorithm(fitness, n, ...
    constraints, pop_size, num_gen)
    stall_gen_limit = 100;
    num_elite = 2;
    fvals = [];

    options = gaoptimset(...
        'OutputFcn', @gaoutfun, ...
        'PopulationType', 'doubleVector', ...
        'Generations',num_gen-1,...
        'EliteCount', num_elite,...
        'TolFun', 0,...
        'TolCon', 0,...
        'PopulationSize', pop_size, ...
        'StallGenLimit', stall_gen_limit,...
        'Display', 'off');
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    LB = constraints(:,1);
    UB = constraints(:,2);
    tic;
    [solution, fval] = ga(fitness, n, A, b, Aeq, beq, LB, UB, nonlcon, options);
    runtime = toc;
    
    function [state,options,optchanged] = gaoutfun(options,state,flag)
        optchanged = false;
        % Find the best objective function
        if ~isempty(state.Best)
            ibest = state.Best(end);
            ibest = find(state.Score == ibest,1,'last');
            bestx = state.Population(ibest,:);
            bestf = fitness(bestx);
            fvals = [fvals; bestf];
        end
    end
end