function [best_pos, fmin, fvals, runtime] = ButterflyOptimizationAlgorithm(fobj, dim, constraints, n, N_iter)

tic;
% n is the population size
% N_iter represnets total number of iterations
p=0.8;                       % probabibility switch
power_exponent=0.1;
sensory_modality=0.01;
Lb = constraints(1,1);
Ub = constraints(1,2);

%Initialize the positions of search agents
Sol=initialization(n,dim,Ub,Lb);

for i=1:n,
    Fitness(i)=fobj(Sol(i,:));
end

% Find the current best_pos
[fmin,I]=min(Fitness);
fvals = fmin;
best_pos=Sol(I,:);
S=Sol; 

% Start the iterations -- Butterfly Optimization Algorithm 
for t=1:N_iter-1,
  
        for i=1:round(n/2), % Loop over all butterflies/solutions
         
          %Calculate fragrance of each butterfly which is correlated with objective function
          Fnew=fobj(S(i,:));
          FP=(sensory_modality*(abs(Fnew)^power_exponent));   
    
          %Global or local search
          if rand<p,    
            dis = rand * rand * best_pos - Sol(i,:);        %Eq. (2) in paper
            S(i,:)=Sol(i,:)+dis*FP;
           else
              % Find random butterflies in the neighbourhood
              epsilon=rand;
              JK=randperm(n);
              dis=epsilon*epsilon*Sol(JK(1),:)-Sol(JK(2),:);
              S(i,:)=Sol(i,:)+dis*FP;                         %Eq. (3) in paper
          end
           
            % Check if the simple limits/bounds are OK
            S(i,:)=simplebounds(S(i,:),Lb,Ub);
          
            % Evaluate new solutions
            if t<N_iter
                Fnew=fobj(S(i,:));  %Fnew represents new fitness values
            end
            
            % If fitness improves (better solutions found), update then
            if (Fnew<=Fitness(i)),
                Sol(i,:)=S(i,:);
                Fitness(i)=Fnew;
            end
           
            % Update the current global best_pos
            if Fnew<=fmin,
                best_pos=S(i,:);
                fmin=Fnew;
            end           
         end
                     
         %Update sensory_modality
         sensory_modality=sensory_modality_NEW(sensory_modality, N_iter);
          
        fvals = [fvals; fmin];
end
    runtime = toc;
end

% Boundary constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb;
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub;
  % Update this new move 
  s=ns_tmp;
end
  
function y=sensory_modality_NEW(x,Ngen)
y=x+(0.025/(x*Ngen));
end

function [X]=initialization(N,dim,up,down)

if size(up,1)==1
    X=rand(N,dim).*(up-down)+down;
end
if size(up,1)>1
    for i=1:dim
        high=up(i);low=down(i);
        X(:,i)=rand(1,N).*(high-low)+low;
    end
end
end

