function [bestnest, fmin, fvals, runtime] = CuckooSearch(fun, dim, constraints, ps, genn)
%% Simple bounds of the search domain  dim
pa=0.25;
Lb = constraints(1,1);
Ub = constraints(1,2);
opti=ones(1,genn);
nest=zeros(ps,dim);
new_nest=zeros(ps,dim);
fitness=10^10*ones(ps,1);
av=zeros(1,genn);
Niter=1;
s=zeros(1,dim);
%% Change this if you want to get better results
% Tolerance
%Tol=1.0e-5;
tic;
% Random initial solutions
nest=rand(ps,dim)*(Ub-Lb)+Lb; 
% Get the current best
[fmin,bestnest,nest,fitness]=get_best_nest(fun, nest,nest,fitness,ps);
fvals = fmin;
av(Niter)=sum(fitness)/ps;
opti(1)=fmin; 
%% Starting iterations
while(Niter < genn)
    % Generate new solutions (but keep the current best)
    new_nest=get_cuckoos(nest,bestnest,Lb,Ub);
    [fnew,best,nest,fitness]=get_best_nest(fun, nest,new_nest,fitness,ps);     
    % Update the counter
    Niter=Niter+1; 
    av(Niter)=sum(fitness)/ps; 
    opti(Niter)=opti(Niter-1);    
    if opti(Niter)>=fnew
        opti(Niter)=fnew;
        bestnest=best;
    end
    
    fvals = [fvals; opti(Niter)];
    % Discovery and randomization
  
    % A fraction of worse nests are discovered with a probability pa
    % Discovered or not -- a status vector
    K=rand(size(nest))>pa;
    % In the real world, if a cuckoo's egg is very similar to a host's eggs, then
    % this cuckoo's egg is less likely to be discovered, thus the fitness should 
    % be related to the difference in solutions.  Therefore, it is a good idea 
    % to do a random walk in a biased way with some random step sizes. 
    % New solution by biased/selective random walks
    stepsize=rand*(nest(randperm(ps),:)-nest(randperm(ps),:));
    new_nest=nest+stepsize.*K;
    for j=1:ps
        s=new_nest(j,:);
        new_nest(j,:)=simplebounds(s,Lb,Ub); 
    end
    if (Niter < genn)
        [fnew,best,nest,fitness]=get_best_nest(fun, nest,new_nest,fitness,ps);

        % Update the counter
        Niter=Niter+1; 
        av(Niter)=sum(fitness)/ps; 
        opti(Niter)=opti(Niter-1);    
        if opti(Niter)>=fnew  % Find the best objective so far  
            opti(Niter)=fnew;
            bestnest=best;
        end
        fvals = [fvals; opti(Niter)];
    end
end %% End of iterations

%% Post-optimization processing
%% Display all the nests
%disp(strcat('Total number of iterations=',num2str(Niter)));
fmin=opti(Niter);
runtime = toc;
end


%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(fun, nest,newnest,fitness,ps)
% Evaluating all new solutions
for j=1:ps
    fnew=fun(newnest(j,:));
    if fnew<=fitness(j),
       fitness(j)=fnew;
       nest(j,:)=newnest(j,:);
    end
end
% Find the current best
[fmin,K]=min(fitness);
best=nest(K,:);
end


%% --------------- All subfunctions are list below ------------------
%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best,Lb,Ub)
% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n,
    s=nest(j,:);
    % This is a simple way of implementing Levy flights
    % For standard random walks, use step=1;
    %% Levy flights by Mantegna's algorithm
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
  
    % In the next equation, the difference factor (s-best) means that 
    % when the solution is the best solution, it remains unchanged.     
    stepsize=0.01*step.*(s-best);
    % Here the factor 0.01 comes from the fact that L/100 should the typical
    % step size of walks/flights where L is the typical lenghtscale; 
    % otherwise, Levy flights may become too aggresive/efficient, 
    % which makes new solutions (even) jump out side of the design domain 
    % (and thus wasting evaluations).
    % Now the actual random walks or flights
    s=s+stepsize.*randn(size(s));
   % Apply simple bounds/limits
   nest(j,:)=simplebounds(s,Lb,Ub);
end
end


% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  
  II=ns_tmp<Lb;
  ns_tmp(II)=Lb;
  
  % Apply the upper bounds 
  JJ=ns_tmp>Ub;
  ns_tmp(JJ)=Ub;
  % Update this new move 
  s=ns_tmp;
end
