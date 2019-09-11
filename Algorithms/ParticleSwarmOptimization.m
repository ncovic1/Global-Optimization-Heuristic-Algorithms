function [gbest, fmin0, fvals, runtime] = ParticleSwarmOptimization(cost_fun, dimensionality, constraints,...
    num_particles, maxite)
    t = cputime;
    
    % objective function (minimization) 
    ofun = cost_fun; 
    
    LB = constraints(:,1);      %lower bounds of variables 
    UB = constraints(:,2);      %upper bounds of variables 

    % pso parameters values 
    m = dimensionality;           % number of variables 
    n = num_particles;       % population size 
    wmax = 0.9;       % inertia weight 
    wmin = 0.4;       % inertia weight 
    c1 = 2;           % acceleration factor 
    c2 = 2;           % acceleration factor 

    % pso main program----------------------------------------------------start 
    maxrun = 1;      % set maximum number of runs need to be 
    for run=1:maxrun 
        % pso initialization----------------------------------------------start 
        for i=1:n 
            for j=1:m 
                x0(i,j)=LB(j)+rand()*(UB(j)-LB(j)); 
            end 
        end 
        x=x0;       % initial population 
        v=0.1*x0;   % initial velocity 
        for i=1:n 
            f0(i,1)=ofun(x0(i,:)); 
        end 
        [fmin0,index0]=min(f0); 
        fvals = fmin0;
        pbest=x0;               % initial pbest 
        gbest=x0(index0,:);     % initial gbest 
        % pso initialization------------------------------------------------end 

        % pso algorithm---------------------------------------------------start 
        ite=1;     
        tolerance=1; 
        while ite<=maxite-1

            w=wmax-(wmax-wmin)*ite/maxite; % update inertial weight 

            % pso velocity updates 
            for i=1:n 
                for j=1:m 
                    v(i,j)=w*v(i,j)+c1*rand()*(pbest(i,j)-x(i,j))... 
                            +c2*rand()*(gbest(1,j)-x(i,j)); 
                end 
            end 

            % pso position update 
            for i=1:n 
                for j=1:m 
                    x(i,j)=x(i,j)+v(i,j); 
                end 
            end 

            % handling boundary violations 
            for i=1:n 
                for j=1:m 
                    if x(i,j)<LB(j) 
                        x(i,j)=LB(j); 
                    elseif x(i,j)>UB(j) 
                        x(i,j)=UB(j); 
                    end 
                end 
            end 

            % evaluating fitness 
            for i=1:n 
                f(i,1)=ofun(x(i,:)); 
            end 

            % updating pbest and fitness 
            for i=1:n 
                if f(i,1)<f0(i,1) 
                    pbest(i,:)=x(i,:); 
                    f0(i,1)=f(i,1); 
                end 
            end 

            [fmin,index]=min(f0);   % finding out the best particle 
            ffmin(ite,run)=fmin;    % storing best fitness 
            ffite(run)=ite;         % storing iteration count 

            % updating gbest and best fitness 
            if fmin<fmin0 
                gbest=pbest(index,:); 
                fmin0=fmin; 
            end     
            fvals = [fvals; fmin0];

            % calculating tolerance 
            if ite>100; 
                tolerance=abs(ffmin(ite-100,run)-fmin0); 
            end 

            ite=ite+1; 
        end 
        % pso algorithm-----------------------------------------------------end 
    end
    
    
    runtime = cputime-t;
    %%
end