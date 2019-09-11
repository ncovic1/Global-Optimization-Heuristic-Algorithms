close all;
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%% Classical benchmark functions
% Unimodal and separable benchmark functions
% 1:  Step                  % n-dim
% 2:  Sphere                % n-dim
% 3:  SumSquares            % n-dim
% 4:  Quartic               % n-dim

% Unimodal and non-separable benchmark functions
% 5:  Beale                 % 2D
% 6:  Easom                 % 2D
% 7:  Matyas                % 2D
% 8:  Colville              % 4D
% 9:  Zakharov              % n-dim
% 10: Schwefel2_22          % n-dim
% 11: Schwefel1_2           % n-dim
% 12: DixonPrice            % n-dim

% Multimodal and separable benchmark functions
% 13: Bohachevsky1          % 2D
% 14: Booth                 % 2D
% 15: HolderTable           % 2D
% 16: Michalewicz2          % 2D 
% 17: Michalewicz5          % 5D
% 18: Michalewicz10         % 10D
% 19: Rastrigin             % n-dim

% Multimodal and non-separable benchmark functions
% 20: Schaffer2             % 2D
% 21: Schaffer4             % 2D
% 22: Schaffer6             % n-dim
% 23: SixHumpCamelBack      % 2D
% 24: Bohachevsky2          % 2D
% 25: Bohachevsky3          % 2D
% 26: Shubert               % 2D
% 27: DropWave              % 2D
% 28: Rosenbrock            % n-dim
% 29: Griewank              % n-dim
% 30: Ackley                % n-dim

%% CEC 2014 benchmark functions (see https://bee22.com/resources/Liang%20CEC2014.pdf)   
% Used functions: [1,2, 5,6,10,16, 17,22, 25,28]    

%% Settings 
n1 = [30,30,30,30, 2,2,2,4,10,30,30,30, 2,2,2,2,5,10,30, 2,2,30,2,2,2,2,2,30,30,30]; % dimensions for classical functions
n2 = 30*ones(1,10);             % dimensions for CEC 2014 functions
fun_type        = 1;            % 1: using classical functions; 2: using CEC 2014 functions
n               = 30;           % dimensionality
functions       = 19;         % numbers of used test functions
num_test        = 3;           % number of testing (repeating)
N               = 1000*n;       % number of points per iteration 
M               = 10;           % number of iterations

%% Algorithms
algorithms = [
    1,...   % Wingsuit Flying Search (WFS)
    0,...   % Genetic Algorithm (GA)
    0,...   % Particle Swarm Optimization (PSO)
    0,...   % Bat Algorithm (BA)
    0,...   % Grey Wolf Optimizer (GWO)
    0,...   % Butterfly Optimization Algorithm (BOA)
    0,...   % Whale Optimization Algorithm (WOA)
    0,...   % Moth Flame Optimization (MFO)
    1,...   % LSHADE
    1,...   % UMOEA
    ];      % 1 if used, 0 if not
algorithms = ones(10,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alg_names = {
    @WingsuitFlyingSearch
    @GeneticAlgorithm
    @ParticleSwarmOptimization
    @BatAlgorithm
    @GreyWolfOptimizer
    @ButterflyOptimizationAlgorithm
    @WhaleOptimizationAlgorithm
    @MothFlameOptimization
    @LSHADE
    @UMOEA
    };
global fun_num;

for ii = 1:length(functions)
    num_alg = sum(algorithms);    
    temp.algorithms = alg_names;
    temp.fun = functions(ii);
    if fun_type == 1
        temp.category = 'classical';
        [f_real, constraints, fun] = SetParametres(functions(ii), n(ii));
        disp(['Test function: ', func2str(fun)]);        
    elseif fun_type == 2
        temp.category = 'CEC2014';
        fun_num = functions(ii);
        f_real = 100*fun_num;
        constraints = 100*[-ones(n(ii),1), ones(n(ii),1)];
        fun = @CEC2014_functions;
        disp(['Test function: ', func2str(fun), ' no. ', num2str(fun_num)]);
    end
    temp.n = n(ii);
    temp.mean_values = zeros(1,num_alg);
    temp.std = zeros(1,num_alg);
    temp.best_values = zeros(1,num_alg);
    temp.worst_values = zeros(1,num_alg);
    temp.mean_runtimes = zeros(1,num_alg);
    temp.values = zeros(num_test, num_alg);
    temp.all_values = zeros(M, num_test, num_alg);
    
    runtimes = zeros(num_test, num_alg);
    num_eval_fun_values = N*M;
    disp(['Dimensionality: ', num2str(n(ii))]);
    disp(' ');

    for i = 1:num_test
        disp(['Testing: ', num2str(i),' of ', num2str(num_test)]);
        f_val = zeros(num_alg, 1);
        x_best = zeros(num_alg, n(ii));
        num = 1;
        for j = 1:length(algorithms)
            if algorithms(j)
                [x_best(num,:), temp.values(i,num), temp.all_values(:,i,num), runtimes(i,num)] = alg_names{j}(fun, n(ii), constraints, N(ii), M);
                temp.values(i,num) = temp.values(i,num)-f_real;
                if abs(temp.values(i,num)) < 1e-8
                    temp.values(i,num) = 0;
                end
                num = num + 1;
            end
        end
    end

    for i=1:num_alg
        temp.mean_values(1,i) = 1/num_test*sum(temp.values(:,i));
        temp.std(1,i) = sqrt(1/num_test*sum((temp.values(:,i)).^2));
        temp.best_values(1,i) = min(temp.values(:,i));
        temp.worst_values(1,i) = max(temp.values(:,i));
        temp.mean_runtimes(1,i) = 1/num_test*sum(runtimes(:,i));
    end

    format shortE;
    disp(' ');
    disp('Mean solution values: '); disp(temp.mean_values(1,:));
    disp('Standard deviation of solution values: '); disp(temp.std(1,:));
    disp('Best solution values: '); disp(temp.best_values(1,:));
    disp('Worst solution values: '); disp(temp.worst_values(1,:));
    disp('Mean runtimes: '); disp(temp.mean_runtimes(1,:));
    disp('---------------------------------------------------------------------------------------------');
    
    data{ii} = temp;
    save data2;
end