function [x_min_m, f_min_m, f_vals, runtime] = WingsuitFlyingSearch(fun, n, bounds, N, M)
    % --------------------------------------------------------------------------
    % Wingsuit Flying Search algorithm
    % Inputs:
        % 'fun' - function needs to be minimized
        % 'n' - dimensionality of the search space
        % 'bounds' - bounds of the search space
        % 'N' - total number of points at each iteration
        % 'M' - total number of the algorithm iterations
    % Outputs:
        % 'x_min_m' - best solution
        % 'f_min_m' - best solution value
        % 'f_vals' - best solutions from each iteration (convergence curve)
        % 'runtime' - total runtime of the algorithm
    % --------------------------------------------------------------------------
    
    % --------------------------------------------------------------------------
    % Wingsuit Flying Search Settings
    v           = 90*rand+10;   % flier's velocity (random between 10 and 100)
    delta_x_min = zeros(1,n);   % minimal discretization step
    plot_res    = 0;            % plot results
    % --------------------------------------------------------------------------
    
    tic;			
    num_x = 1;		% aux. variable for plotting
    N = N - 2;		% number of considered points (two more are added in the end of each iteration)
    
    N0 = ceil((N)^(1/n));	% number of points in regular grid per one dimension
    step = zeros(1, n);		% initial discretization step
    for i = 1:n
        step(i) = (bounds(i,2)-bounds(i,1))/N0;
    end

    % Generating initial Halton points
    halton_set = haltonset(n, 'Skip', 0, 'Leap', 0);
    halton_set = scramble(halton_set, 'RR2');
    x_temp = ceil(rand*1e6);
    x = halton_set(x_temp:x_temp+N-1, :);
    for i = 1:n
        x(:,i) = x(:,i).*(bounds(i,2)-bounds(i,1))+bounds(i,1);
    end        
    
    f = GetValues(x, n, bounds, fun);	% values of the initial points
    m = 1;                      % number of current iteration number
    f_max = max(f);				% maximal found value		
    a_m = f_max;				% flier's altitude
    
    [f_min_m, index] = min(f);  % current solution value
    x_min_m = x(index,:);       % current solution
    
    % Generating centroid and random point and updating current solution if necessary
    [x_centr, x_rand] = GenerateCRP(x, n, f, f_min_m, a_m); 
    f_temp = GetValues([x_centr; x_rand], n, bounds, fun); 
    x = [x; x_centr; x_rand];
    f = [f; f_temp];
    [x_min_m, f_min_m] = UpdateMin(x(N+1:N+2,:), f(N+1:N+2), x_min_m, f_min_m);
    f_vals(m) = f_min_m;    
    
    % Plotting results
    if plot_res     
        num_x = PlotRes(x, f, n, num_x, bounds, f_min_m);
    end
    
    % Another format of discretization step
    delta_x1 = zeros(3, n);     
    delta_x1(1,:) = step;
    delta_x1(2,:) = zeros(1,n);
    delta_x1(3,:) = -step;
    
    SC = CheckSC(delta_x1, delta_x_min, m, M);  % checking stopping condition
    %%    
    while(SC == false)
        m = m+1;    
        alpha_m = 1-v^(-(m-1)/(M-1));       % search sharpness
        P_max_m = ceil(alpha_m*N);          % maximal number of neighborhood points
        N_m = ceil(2*N/P_max_m);            % number of considered points
        if N_m > length(f)
            N_m = length(f);
        end
        delta_x_m = (1-alpha_m)*delta_x1;   % discretization step
        
        % Sorting points in ascending order w.r.t. their solution values
        for i = 1:N_m
            for j = i+1:length(f)
                if f(j)<f(i)
                    f_temp = f(j);
                    x_temp = x(j,:);
                    f(j) = f(i);
                    x(j,:) = x(i,:);
                    f(i) = f_temp;
                    x(i,:) = x_temp;
                end
            end
        end
        f = f(1:N_m);       % taking only first N_m values
        x = x(1:N_m,:);     % taking only first N_m points
        a_m = f(N_m);       % flier's current altitude

        % Generating neighborhood size for each point
        P_m = zeros(N_m,1);     
        S = 0; 
        for i = 1:N_m
            P_m_temp = ceil(P_max_m*(1-(i-1)/(N_m-1)));
            if S+P_m_temp <= N
                P_m(i) = P_m_temp;
            else
                P_m(i) = N-S;
            end
            S = S+P_m(i);
        end
        if N_m-S > 0
            P_m(i) = N_m-S;
        end
        
        % Generating new points
        x_size = size(x,1);
        for i = 1:length(f)
            if P_m(i)>0
                direction = x_min_m-x(i,:);
                row = zeros(n, 2);
                for j = 1:n
                    if direction(j)>0
                        direction(j) = 1;
                        row(j,:) = [1,2];
                    elseif direction(j)<0
                        direction(j) = -1;
                        row(j,:) = [2,3];
                    else
                        row(j,:) = [1,3];
                    end
                end 
                x_selected = zeros(P_m(i), n);
                x_temp = ceil(rand*1e6);                        
                neighborhood = halton_set(x_temp:x_temp+P_m(i)-1, :); 
                for j = 1:n
                    x_selected(:,j) = (2*neighborhood(:,j)-1+direction(j))*delta_x_m(1,j)+...
                        x(i,j)*ones(P_m(i),1);
                end
                x = [x; x_selected];
            end
        end
        
        f_temp = GetValues(x(x_size+1:size(x,1),:), n, bounds, fun);   % new points values
        f = [f; f_temp];    % all values: old and new
        
        [f_min_m, index] = min(f);  % new solution value
        x_min_m = x(index, :);      % new solution

        % Generating centroid and random point and updating current solution if necessary
        [x_centr, x_rand] = GenerateCRP(x, n, f, f_min_m, a_m); 
        f_temp = GetValues([x_centr; x_rand], n, bounds, fun);
        x = [x; x_centr; x_rand];
        f = [f; f_temp];
        [x_min_m, f_min_m] = UpdateMin(x(N+1:N+2,:), f(N+1:N+2), x_min_m, f_min_m);
        f_vals(m) = f_min_m;
        
        % Plotting results
        if plot_res     
            num_x = PlotRes(x, f, n, num_x, bounds, f_min_m);
        end
        
        SC = CheckSC(delta_x_m, delta_x_min, m, M);     % checking stopping condition
    end
    runtime = toc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Get values 'f' for 'n'-dim. points 'x' of function 'fun' with constraints 'constraints' of serach space
function f = GetValues(x, n, constraints, fun)
    f = zeros(size(x,1),1);
    for i = 1:size(x,1)
        for j = 1:n
            if x(i,j) <= constraints(j,1)
                x(i,j) = 2*constraints(j,1) - x(i,j);
            elseif x(i,j) >= constraints(j,2)
                x(i,j) = 2*constraints(j,2) - x(i,j); 
            end
        end        
        f(i) = fun(x(i,:));
    end
end

% Update minimum 'x_min_m' and its value 'f_min_m' considering points 'x' and its values 'f'
function [x_min_m, f_min_m] = UpdateMin(x, f, x_min_m, f_min_m)
    for i = 1:length(f)
        if f(i) < f_min_m
            f_min_m = f(i);
            x_min_m = x(i,:);
            break;
        end
    end
end

% Generate centroid 'x_centr' and random point 'x_rand' from 'n'-dim.
% points 'x', its values 'f', current solution value 'f_min_m', and current
% flier's attitude 'a_m'
function [x_centr, x_rand] = GenerateCRP(x, n, f, f_min_m, a_m)
    gama = 1-(f-f_min_m)/(a_m-f_min_m);
    gama = (gama > 0).*gama;
    x_centr = zeros(1,n);
    current_constr = zeros(n,2);
    for i = 1:n
        x_centr(i) = (x(:,i))'*gama/sum(gama);
        current_constr(i,:) = [min(x(:,i)); max(x(:,i))];
    end
    x_rand = ((current_constr(:,2)-current_constr(:,1)).*rand(n,1)+current_constr(:,1))';
end

% Check stopping condition 'SC' comparing current discretization step
% 'delta_x_m' with its minimal 'delta_x_min', and current iteration number
% 'm' with total ones 'M'
function SC = CheckSC(delta_x_m, delta_x_min, m, M)
    SC = false; 
    for i = 1:length(delta_x_min)
        if abs(delta_x_m(1,i)) < delta_x_min(i)
            SC = true;
            disp('Discretization step is less than minimal. Algorithm terminated.');             
            break;
        end
    end
    if m == M
        SC = true;
        %disp('*********** algorithm terminated ***********'); disp(' ');
    end
end

% Plot results using 'n'-dim. points 'x', its values 'f', current
% considered point number 'num_x', constraints of search space
% 'constraints', and current solution value 'f_min_m'
function num_x = PlotRes(x, f, n, num_x, constraints, f_min_m)
    global x_real;
    for i = 1:length(f)
        if i < length(f)-3
            color = 'b';
            marker_size = 5;
        elseif i == length(f)-3
            color = 'r';
            marker_size = 20;
        elseif i == length(f)-2
            color = 'g';
            marker_size = 20;
        elseif i == length(f)-1
            color = 'm';
            marker_size = 20;
        else
            color = 'k';
            marker_size = 20;
        end
        if f(i) < 1e10
            if mod(n,2)==0
                n_temp = n;
            else
                n_temp = n-3;
            end
            for k = 1:2:n_temp
                figure((k+1)/2);
                if mod(num_x, 500) == 1
                    for ii = 1:size(x_real,1)
                        plot(x_real(ii,k), x_real(ii,k+1), 'green x', 'MarkerSize', 20); hold on;                        
                    end
                    xlabel(strjoin({'$x_{', num2str(k),'}$'})); 
                    ylabel(strjoin({'$x_{', num2str(k+1),'}$'}));
                    axis([constraints(k,:), constraints(k+1,:)]);
                    grid on; 
                end
                plot(x(i,k), x(i,k+1), [color, '.'], 'MarkerSize', marker_size);
                hold on;  
                drawnow;
            end
            if n > n_temp
                figure(n_temp+1);
                if mod(num_x, 50) == 1
                    for ii = 1:size(x_real,1)
                        plot3(x_real(ii,n_temp+1), x_real(ii,n_temp+2), x_real(ii,n_temp+3), ...
                            'green x', 'MarkerSize', 20); hold on;
                    end
                    xlabel(strjoin({'$x_{', num2str(n_temp+1),'}$'}));
                    ylabel(strjoin({'$x_{', num2str(n_temp+2),'}$'}));
                    zlabel(strjoin({'$x_{', num2str(n_temp+3),'}$'}));
                    axis([constraints(n_temp+1,:), constraints(n_temp+2,:), constraints(n_temp+3,:)]);
                    grid on; 
                end
                plot3(x(i,n_temp+1), x(i,n_temp+2), x(i,n_temp+3), [color, '.'], 'MarkerSize', marker_size);
                hold on;  
                drawnow;
            end
            figure(n_temp+2)
            plot(num_x, f(i),'r.');
            xlabel('Point ordinal number');
            ylabel('Function value');
            axis([0, num_x, f_min_m, f_min_m+10*abs(f_min_m)]);
            grid on; hold on;  
            drawnow;
        end
        num_x = num_x+1;
    end
end