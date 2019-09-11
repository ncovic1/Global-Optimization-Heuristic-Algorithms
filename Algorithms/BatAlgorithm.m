function  [BestPositions, fmin, Convergence_curve, runtime] = BatAlgorithm(fobj, dim, constraints, N, Max_iter)
%%
%Max_iter=20;            % maximum generations
%N=10;                   % BAT numbers
%dim=10;
lb=constraints(:,1)';
ub=constraints(:,2)';
Fmax=2;                 %maximum frequency
Fmin=0;                 %minimum frequency
A=rand(N,1);            %loudness for each BAT
r=rand(N,1);            %pulse emission rate for each BAT
alpha=0.5;              %constant for loudness update
gamma=0.5;              %constant for emission rate update
ro=0.001;                 %initial pulse emission rate

tic;
% Initializing arrays
F=zeros(N,1);           % Frequency
v=zeros(N,dim);           % Velocities
% Initialize the population
x=initializationb(N,Max_iter,dim,ub,lb);
Convergence_curve=zeros(1,Max_iter);
%calculate the initial solution for initial positions

for ii=1:N
    fitness(ii)=fobj(x(ii,:));
end

[fmin,index]=min(fitness);          %find the initial best fitness value
Convergence_curve(1) = fmin;
bestsol=x(index,:);                 %find the initial best solution for best fitness value
%%
iter=1;             % start the loop counter
while iter<Max_iter                              %start the loop for iterations
    for ii=1:size(x)
        F(ii)=Fmin+(Fmax-Fmin)*rand;              %randomly chose the frequency
        v(ii,:)=v(ii,:)+(x(ii,:)-bestsol)*F(ii);  %update the velocity
        x(ii,:)=x(ii,:)+v(ii,:);                  %update the BAT position
        %         x(ii,:)=round(x(ii,:));
        % Apply simple bounds/limits
        Flag4up=x(ii,:)>ub;
        Flag4low=x(ii,:)<lb;
        x(ii,:)=(x(ii,:).*(~(Flag4up+Flag4low)))+ub.*Flag4up+lb.*Flag4low;
        %check the condition with r
        if rand>r(ii)
            % The factor 0.001 limits the step sizes of random walks
            %               x(ii,:)=bestsol+0.001*randn(1,dim);
            eps=-1+(1-(-1))*rand;
            x(ii,:)=bestsol+eps*mean(A);
        end
        fitnessnew=fobj(x(ii,:));  % calculate the objective function
        % Update if the solution improves, or not too loud
        if (fitnessnew<=fitness(ii)) && (rand<A(ii)) ,
            
            fitness(ii)=fitnessnew;
            A(ii)=alpha*A(ii);
            r(ii)=ro*(1-exp(-gamma*iter));
        end
        if fitnessnew<=fmin,
            bestsol=x(ii,:);
            fmin=fitnessnew;
        end
        
    end
    
    iter=iter+1;                                  % update the while loop counter    
    Convergence_curve(iter)=  fmin;
end
%
[bestfit]=(fmin);
BestPositions=bestsol;

runtime = toc;
end


% This function initialize the first population of search agents
function x=initializationb(N,Max_iter,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    x=rand(N,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        x(:,i)=rand(N,1).*(ub_i-lb_i)+lb_i;
    end
end
end
