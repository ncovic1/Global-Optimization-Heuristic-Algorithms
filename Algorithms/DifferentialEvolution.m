function [best_x, best_val, BestCost, runtime] = DifferentialEvolution(CostFunction, nVar, constraints, nPop, MaxIt)

tic;
%% Problem Definition
%CostFunction=@(x) Sphere(x);    % Cost Function
%nVar=20;            % Number of Decision Variables
VarSize=[1 nVar];   % Decision Variables Matrix Size
VarMin= constraints(1,1);          % Lower Bound of Decision Variables
VarMax= constraints(1,2);          % Upper Bound of Decision Variables
%% DE Parameters
%MaxIt=1000;      % Maximum Number of Iterations
%nPop=50;        % Population Size
beta_min=0.2;   % Lower Bound of Scaling Factor
beta_max=0.8;   % Upper Bound of Scaling Factor
pCR=0.2;        % Crossover Probability
%% Initialization
empty_individual.Position=[];
empty_individual.Cost=[];
BestSol.Cost=inf;
pop=repmat(empty_individual,nPop,1);
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    pop(i).Cost=CostFunction(pop(i).Position);
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
    end
    
end
BestCost=zeros(MaxIt,1);
BestCost(1) = BestSol.Cost;
%% DE Main Loop
for it=2:MaxIt
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation
        %beta=unifrnd(beta_min,beta_max);
        beta=unifrnd(beta_min,beta_max,VarSize);
        y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
        y = max(y, VarMin);
		y = min(y, VarMax);
		
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        NewSol.Position=z;
        NewSol.Cost=CostFunction(NewSol.Position);
        
        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol=pop(i);
            end
        end
        
    end
    
    % Update Best Cost
    BestCost(it)=BestSol.Cost;
    best_x = BestSol.Position;
    best_val = BestSol.Cost;
    
    runtime = toc;
    
end