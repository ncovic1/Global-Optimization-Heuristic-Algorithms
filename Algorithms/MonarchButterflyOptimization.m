
%% Monarch Butterfly Optimization (MBO)
% Author: Gai-Ge Wang
% Email: gaigewang@163.com
%             gaigewang@gmail.com
% Main paper:
% Gai-Ge Wang, Suash Deb, and Zhihua Cui, Monarch Butterfly Optimization.
% Neural Computing and Applications, in press.
% DOI: 10.1007/s00521-015-1923-y
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%% Notes:
% Different run may generate different solutions, this is determined by
% the the nature of metaheuristic algorithms.
%%
function [x_best, f_best, MinCost, runtime] = MonarchButterflyOptimization(costFun, numVar, bounds, popsize, Maxgen)
% Monarch Butterfly Optimization (MBO) software for minimizing a general function
% The fixed Function Evaluations (FEs) is considered as termination condition.
% INPUTS: ProblemFunction is the handle of the function that returns
%         the handles of the initialization, cost, and feasibility functions.
%         DisplayFlag = true or false, whether or not to display and plot results.
%         ProbFlag = true or false, whether or not to use probabilities to update emigration rates.
%         RandSeed = random number seed
% OUTPUTS: MinCost = array of best solution, one element for each generation
%          Hamming = final Hamming distance between solutions
% CAVEAT: The "ClearDups" function that is called below replaces duplicates with randomly-generated
%         individuals, but it does not then recalculate the cost of the replaced individuals.
tic
DisplayFlag = false;
RandSeed = round(sum(100*clock));
    
[OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
    MaxParValue, MinParValue, Population] = Init(DisplayFlag, @FUN, RandSeed, costFun, numVar, bounds, popsize, Maxgen);
nEvaluations = OPTIONS.popsize;
% % % % % % % % % % % %             Initial parameter setting          % % % % % % % % % % % %%%%
%% Initial parameter setting
Keep = 2; % elitism parameter: how many of the best habitats to keep from one generation to the next
maxStepSize = 1.0;        %Max Step size
partition = OPTIONS.partition;
numButterfly1 = ceil(partition*OPTIONS.popsize);  % NP1 in paper
numButterfly2 = OPTIONS.popsize - numButterfly1; % NP2 in paper
period = 1.2; % 12 months in a year
Land1 = zeros(numButterfly1, OPTIONS.numVar);
Land2 = zeros(numButterfly2, OPTIONS.numVar);
BAR = partition; % you can change the BAR value in order to get much better performance
% % % % % % % % % % % %       End of Initial parameter setting       % % % % % % % % % % % %%
%%
% % % % % % % % % % % %             Begin the optimization loop        % % % % % % % % % %%%%
% Begin the optimization loop
GenIndex = 1;
% for GenIndex = 1 : OPTIONS.Maxgen
while nEvaluations< OPTIONS.MaxFEs
    % % % % % % % % % % % %            Elitism Strategy           % % % % % % % % % % % %%%%%
    %% Save the best monarch butterflis in a temporary array.
    for j = 1 : Keep
        chromKeep(j,:) = Population(j).chrom;
        costKeep(j) = Population(j).cost;
    end
    % % % % % % % % % % % %       End of  Elitism Strategy      % % % % % % % % % % % %%%%
    %%
    
    % % % % % % % % % % % %    Divide the whole population into two subpopulations % % % %%%
    %% Divide the whole population into Population1 (Land1) and Population2 (Land2)
    % according to their fitness.
    % The monarch butterflis in Population1 are better than or equal to Population2.
    % Of course, we can randomly divide the whole population into Population1 and Population2.
    % We do not test the different performance between two ways.
    for popindex = 1 : OPTIONS.popsize
        if popindex <= numButterfly1
            Population1(popindex).chrom = Population(popindex).chrom;
        else
            Population2(popindex-numButterfly1).chrom = Population(popindex).chrom;
        end
    end
    % % % % % % % % % % %    End of Divide the whole population into two subpopulations  % % %%%
    %%
    
    % % % % % % % % % % % %%            Migration operator          % % % % % % % % % % % %%%%
    %% Migration operator
    for k1 = 1 : numButterfly1
        for parnum1 = 1 : OPTIONS.numVar
            r1 = rand*period;
            if r1 <= partition
                r2 = round(numButterfly1 * rand + 0.5);
                Land1(k1,parnum1) = Population1(r2).chrom(parnum1);
            else
                r3 = round(numButterfly2 * rand + 0.5);
                Land1(k1,parnum1) = Population2(r3).chrom(parnum1);
            end
        end %% for parnum1
        NewPopulation1(k1).chrom =  Land1(k1,:);
    end  %% for k1
    % % % % % % % % % % % %%%       End of Migration operator      % % % % % % % % % % % %%%
    %%
    
    % % % % % % % % % % % %             Evaluate NewPopulation1          % % % % % % % % % % % %%
    %% Evaluate NewPopulation1
    SavePopSize = OPTIONS.popsize;
    OPTIONS.popsize = numButterfly1;
    % Make sure each individual is legal.
    NewPopulation1 = FeasibleFunction(OPTIONS, NewPopulation1);
    % Calculate cost
    NewPopulation1 = CostFunction(OPTIONS, NewPopulation1, costFun);
    % the number of fitness evaluations
    nEvaluations = nEvaluations +  OPTIONS.popsize;
    OPTIONS.popsize = SavePopSize;
    % % % % % % % % % % % %       End of Evaluate NewPopulation1      % % % % % % % % % % % %%
    %%
    
    % % % % % % % % % % % %             Butterfly adjusting operator          % % % % % % % % % % % %%
    %% Butterfly adjusting operator
    for k2 = 1 : numButterfly2
        scale = maxStepSize/(GenIndex^2); %Smaller step for local walk
        StepSzie = ceil(exprnd(2*OPTIONS.Maxgen,1,1));
        delataX = LevyFlight(StepSzie,OPTIONS.numVar);
        for parnum2 = 1:OPTIONS.numVar,
            if (rand >= partition)
                Land2(k2,parnum2) = Population(1).chrom(parnum2);
            else
                r4 = round(numButterfly2*rand + 0.5);
                Land2(k2,parnum2) =  Population2(r4).chrom(1);
                if (rand > BAR) % Butterfly-Adjusting rate
                    Land2(k2,parnum2) =  Land2(k2,parnum2) +  scale*(delataX(parnum2)-0.5);
                end
            end
        end  %% for parnum2
        NewPopulation2(k2).chrom =  Land2(k2,:);
    end %% for k2
    % % % % % % % % % % % %       End of Butterfly adjusting operator      % % % % % % % % % % % %
    %%
    
    % % % % % % % % % % % %             Evaluate NewPopulation2          % % % % % % % % % % % %%
    %% Evaluate NewPopulation2
    SavePopSize = OPTIONS.popsize;
    OPTIONS.popsize = numButterfly2;
    % Make sure each individual is legal.
    NewPopulation2 = FeasibleFunction(OPTIONS, NewPopulation2);
    % Calculate cost
    NewPopulation2 = CostFunction(OPTIONS, NewPopulation2, costFun);
    % the number of fitness evaluations
    nEvaluations = nEvaluations +  OPTIONS.popsize;
    OPTIONS.popsize = SavePopSize;
    % % % % % % % % % % % %       End of Evaluate NewPopulation2      % % % % % % % % % % % %%
    %%
    
    % % % % % % %  Combine two subpopulations into one and rank monarch butterflis       % % % % % %
    %% Combine Population1 with Population2 to generate a new Population
    Population = CombinePopulation(OPTIONS, NewPopulation1, NewPopulation2);
    % Sort from best to worst
    Population = PopSort(Population);
    % % % % % %     End of Combine two subpopulations into one and rank monarch butterflis  % %% % %
    %%
    
    % % % % % % % % % % % %            Elitism Strategy          % % % % % % % % % % % %%% %% %
    %% Replace the worst with the previous generation's elites.
    n = length(Population);
    for k3 = 1 : Keep
        Population(n-k3+1).chrom = chromKeep(k3,:);
        Population(n-k3+1).cost = costKeep(k3);
    end % end for k3
    % % % % % % % % % % % %     End of  Elitism Strategy      % % % % % % % % % % % %%% %% %
    %%
    
    % % % % % % % % % %           Precess and output the results          % % % % % % % % % % % %%%
    % Sort from best to worst
    Population = PopSort(Population);
    % Compute the average cost
    [AverageCost, nLegal] = ComputeAveCost(Population);
    % Display info to screen
    MinCost = [MinCost Population(1).cost];
    AvgCost = [AvgCost AverageCost];
    if DisplayFlag
        disp(['The best and mean of Generation # ', num2str(GenIndex), ' are ',...
            num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
    end
    % % % % % % % % % % %    End of Precess and output the results     %%%%%%%%%% %% %
    %%
    
    %% Update generation number
    GenIndex = GenIndex+1;
    
end % end for GenIndex
f_best = min(MinCost);
x_best = Population(1).chrom;
runtime = toc;
% % % % % % % % % %     End of Monarch Butterfly Optimization implementation     %%%% %% %
%%
function [delataX] = LevyFlight(StepSize, Dim)
%Allocate matrix for solutions
delataX = zeros(1,Dim);
%Loop over each dimension
for i=1:Dim
    % Cauchy distribution
    fx = tan(pi * rand(1,StepSize));
    delataX(i) = sum(fx);
end



function [Population, indices] = PopSort(Population)
% Sort the population members from best to worst
popsize = length(Population);
Cost = zeros(1, popsize);
indices = zeros(1, popsize);
for i = 1 : popsize
    Cost(i) = Population(i).cost;
end
[Cost, indices] = sort(Cost, 2, 'ascend');
Chroms = zeros(popsize, length(Population(1).chrom));
for i = 1 : popsize
    Chroms(i, :) = Population(indices(i)).chrom;
end
for i = 1 : popsize
    Population(i).chrom = Chroms(i, :);
    Population(i).cost = Cost(i);
end



function [OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
    MaxParValue, MinParValue, Population] = Init(DisplayFlag, ProblemFunction, RandSeed, costFun, numVar, constraints, popsize, Maxgen)
% Initialize population-based optimization software.
% WARNING: some of the optimization routines will not work if population size is odd.
OPTIONS.popsize = popsize; % total population size
OPTIONS.Maxgen = Maxgen; % generation count limit
OPTIONS.numVar =  numVar; % number of variables in each population member
OPTIONS.MaxFEs = popsize*Maxgen; % number of Function Evaluations (FEs)
OPTIONS.partition = 5/12;  % the percentage of population for MBO
if ~exist('RandSeed', 'var')
    RandSeed = round(sum(100*clock));
end
rand('state', RandSeed); % initialize random number generator
if DisplayFlag
    disp(['random # seed = ', num2str(RandSeed)]);
end
% Get the addresses of the initialization, cost, and feasibility functions.
[InitFunction, CostFunction, FeasibleFunction] = ProblemFunction();
% Initialize the population.
[MaxParValue, MinParValue, Population, OPTIONS] = InitFunction(OPTIONS, constraints);
% Make sure the population does not have duplicates. 
Population = ClearDups(Population, MaxParValue, MinParValue);
% Compute cost of each individual  
Population = CostFunction(OPTIONS, Population, costFun);
% Sort the population from most fit to least fit
Population = PopSort(Population);
% Compute the average cost
AverageCost = ComputeAveCost(Population);
% Display info to screen
MinCost = [Population(1).cost];
AvgCost = [AverageCost];
if DisplayFlag
    disp(['The best and mean of Generation # 0 are ', num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
end
return;



function [AveCost, nLegal] = ComputeAveCost(Population)
% Compute the average cost of all legal individuals in the population.
% OUTPUTS: AveCost = average cost
%          nLegal = number of legal individuals in population
% Save valid population member fitnesses in temporary array
Cost = [];
nLegal = 0;
for i = 1 : length(Population)
    if Population(i).cost < inf
        Cost = [Cost Population(i).cost];
        nLegal = nLegal + 1;
    end
end
% Compute average cost.
AveCost = mean(Cost);
return;





function Population1 = CombinePopulation(OPTIONS, Population1, Population2)
numButterfly1 = ceil(OPTIONS.partition*OPTIONS.popsize);
for i = 1: OPTIONS.popsize - numButterfly1
    Population1(numButterfly1 + i).chrom = Population2(i).chrom;
    Population1(numButterfly1 + i).cost = Population2(i).cost;
end



function [Population] = ClearDups(Population, MaxParValue, MinParValue)
% Make sure there are no duplicate individuals in the population.
% This logic does not make 100% sure that no duplicates exist, but any duplicates that are found are
% randomly mutated, so there should be a good chance that there are no duplicates after this procedure.
for i = 1 : length(Population)
    Chrom1 = sort(Population(i).chrom);
    for j = i+1 : length(Population)
        Chrom2 = sort(Population(j).chrom);
        if isequal(Chrom1, Chrom2)
            parnum = ceil(length(Population(j).chrom) * rand);
            Population(j).chrom(parnum) = floor(MinParValue(parnum)...
                + (MaxParValue(parnum) - MinParValue(parnum) + 1) * rand); 
        end
    end
end
return;




function [InitFunction, CostFunction, FeasibleFunction] = FUN
InitFunction = @FUNInit;
CostFunction = @FUNCost;
FeasibleFunction = @FUNFeasible;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MaxParValue, MinParValue, Population, OPTIONS] = FUNInit(OPTIONS, constraints)
global MinParValue MaxParValue
Granularity = 0.1;
MinParValue = constraints(1,1)*ones(1,OPTIONS.numVar);
MaxParValue = constraints(1,2)*ones(1,OPTIONS.numVar);
% Initialize population
for popindex = 1 : OPTIONS.popsize
    chrom = MinParValue + (MaxParValue - MinParValue + 1) .* rand(1,OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = false;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = FUNCost(OPTIONS, Population, costFun)
% Compute the cost of each member in Population
global MinParValue MaxParValue
popsize = OPTIONS.popsize;
p = OPTIONS.numVar;
for popindex = 1 : popsize
    Population(popindex).cost = costFun(Population(popindex).chrom);    
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = FUNFeasible(OPTIONS, Population)
global MinParValue MaxParValue
for i = 1 : OPTIONS.popsize
    for k = 1 : OPTIONS.numVar
        Population(i).chrom(k) = max(Population(i).chrom(k), MinParValue(k));
        Population(i).chrom(k) = min(Population(i).chrom(k), MaxParValue(k));
    end
end
return;
