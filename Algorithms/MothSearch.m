
%% Moth Search (MS) Algorithm
% Author: Gai-Ge Wang
% Email: gaigewang@gmail.com
%             gaigewang@163.com
% Main paper:
% Gai-Ge Wang, Moth search algorithm: a bio-inspired metaheuristic
% algorithm for global optimization problems.
% Memetic Computing.
% DOI: 10.1007/s12293-016-0212-3
% http://rd.springer.com/article/10.1007%2Fs12293-016-0212-3
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%% Notes:
% Different run may generate different solutions, this is determined by
% the the nature of metaheuristic algorithms.
%%
function [x_best, f_best, MinCost, runtime] = MothSearch(costFun, numVar, bounds, popsize, Maxgen)
% Moth Search (MS) Algorithm software for minimizing a general function
% The fixed generation is considered as termination condition.
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
Keep = 2; % elitism parameter: how many of the best moths to keep from one generation to the next
maxStepSize = 1.0;        %Max Step size
partition = OPTIONS.partition;
numMoth1 = ceil(partition*OPTIONS.popsize);  % NP1 in paper
numMoth2 = OPTIONS.popsize - numMoth1; % NP2 in paper
TempChone = zeros(1, OPTIONS.numVar);
goldenRatio = (sqrt(5)-1)/2; % you can change this Ratio so as to get much better performance
% % % % % % % % % % % %       End of Initial parameter setting       % % % % % % % % % % % %%
%%
% % % % % % % % % % % %             Begin the optimization loop        % % % % % % % % % %%%%
% Begin the optimization loop
GenIndex = 1;
while nEvaluations< OPTIONS.MaxFEs
    
    % % % % % % % % % % % %            Elitism Strategy           % % % % % % % % % % % %%%%%
    %% Save the best monarch butterflis in a temporary array.
    for j = 1 : Keep
        chromKeep(j,:) = Population(j).chrom;
        costKeep(j) = Population(j).cost;
    end
    % % % % % % % % % % % %       End of  Elitism Strategy      % % % % % % % % % % % %%%%
    %%
    
    % % % % % % % % % % % %%            L¨¦vy flights          %% % % % % % % % % % % % % %%%%
    %% Migration operator
    for k1 = 1 : numMoth1
        scale = maxStepSize/(GenIndex^2); %Smaller step for local walk
        delataX = LevyWalk(OPTIONS.numVar);
        Population(k1).chrom = Population(k1).chrom + scale*delataX;
    end  %% for k1
    % % % % % % % % % % % %%%       End of L¨¦vy flights      % % % % % % % % % % % %%%%%
    %%
    
    % % % % % % % % % % % %             Flying in a straingt line          % % % % % % % % % % % %%
    %% Flying in a straingt line
    for k2 = 1 : numMoth2
        for parnum = 1:OPTIONS.numVar
            if (rand >= 0.5)
                TempChone(parnum) = Population(k1+k2).chrom(parnum)...
                    + goldenRatio*(Population(1).chrom(parnum) - Population(k1+k2).chrom(parnum));
            else
                TempChone(parnum) = Population(k1+k2).chrom(parnum)...
                    + (1/goldenRatio)*(Population(1).chrom(parnum) - Population(k1+k2).chrom(parnum));
            end
        end  %% for parnum
        Population(k1+k2).chrom =  rand*TempChone;
    end %% for k2
    % % % % % % % % % % % %       End of  Flying in a straingt line      % % % % % % % % % % % %
    %%
    
    % % % % % % % % % % % %          Evaluate new Population       % % % % % % %  % % %%  % % 
    % Make sure each individual is legal.
    Population = FeasibleFunction(OPTIONS, Population);
    % Calculate cost
    Population = CostFunction(OPTIONS, Population, costFun);
    % the number of fitness evaluations
    nEvaluations = nEvaluations +  OPTIONS.popsize;
    % Sort from best to worst
    Population = PopSort(Population);
    % % % % % %% % % % % %       End of Evaluate new Population   % %% % %% % % % %  % % 
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
% % % % % % % % % %     End of Moth Search (MS) Algorithm implementation     %%%% %% %
%%




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



function delataX=LevyWalk(d)
beta = 1.5;
%Eq. (2.23)
sigma=(gamma(1+beta)*sin(pi*(beta-1)/2)/(gamma((beta)/2)*(beta-1)*2^((beta-2)/2)))^(1/(beta-1));
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/(beta-1)); %Eq. (2.21)
delataX=0.01*step;



function [OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
    MaxParValue, MinParValue, Population] = Init(DisplayFlag, ProblemFunction, RandSeed, costFun, numVar, constraints, popsize, Maxgen)
% Initialize population-based optimization software.
% WARNING: some of the optimization routines will not work if population size is odd.

OPTIONS.popsize = popsize; % total population size
OPTIONS.Maxgen = Maxgen; % generation count limit
OPTIONS.numVar =  numVar; % number of variables in each population member
OPTIONS.MaxFEs = popsize*Maxgen; % number of Function Evaluations (FEs)

OPTIONS.partition = 0.5;  % the percentage of population for MBO
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



function [Population] = ClearDups(Population, MaxParValue, MinParValue)
% Make sure there are no duplicate individuals in the population.
% This logic does not make 100% sure that no duplicates exist, but any duplicates that are found are
% randomly mutated, so there should be a good chance that there are no duplicates after this procedure.
if length(MaxParValue) == 1
    for i = 1 : length(Population)
        Chrom1 = sort(Population(i).chrom);
        for j = i+1 : length(Population)
            Chrom2 = sort(Population(j).chrom);
            if isequal(Chrom1, Chrom2)
                parnum = ceil(length(Population(j).chrom) * rand);
                Population(j).chrom(parnum) = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand);
            end
        end
    end
else
    for i = 1 : length(Population)
        Chrom1 = sort(Population(i).chrom);
        for j = i+1 : length(Population)
            Chrom2 = sort(Population(j).chrom);
            if isequal(Chrom1, Chrom2)
                parnum = ceil(length(Population(j).chrom) * rand);
                Population(j).chrom(parnum) = floor(MinParValue(parnum) ...
                    + (MaxParValue(parnum) - MinParValue(parnum) + 1) * rand);
            end
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
