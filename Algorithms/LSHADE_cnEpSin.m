%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB/Octave source code of LSHADE_cnEpSin which is a new version of LSHADE-EpSin.
%% Please see the following papers:
%% 1. LSHADE_cnEpSin:
%%     Noor H. Awad, Mostafa Z. Ali, Ponnuthurai N. Suganthan, Ensemble Sinusoidal Differential Covariance Matrix Adaptation with Euclidean Neighborhood  for Solving CEC2017 Benchmark Problems, in Proc. IEEE Congr. Evol. Comput. CEC 2017, June, Donostia - San Sebastián, Spain

%% 2. LSHADE-EpSin:
%%    Noor H. Awad, Mostafa Z. Ali, Ponnuthurai N. Suganthan and Robert G. Reynolds: An Ensemble Sinusoidal Parameter Adaptation incorporated with L-SHADE for Solving CEC2014 Benchmark Problems, in Proc. IEEE Congr. Evol. Comput. CEC 2016, Canada, July, 2016

%% About L-SHADE, please see following papers:
%% Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.
%%  J. Zhang, A.C. Sanderson: JADE: Adaptive differential evolution with optional external archive,” IEEE Trans Evol Comput, vol. 13, no. 5, pp. 945–958, 2009

function [bsf_solution, bsf_fit_var, run_funcvals, runtime] = LSHADE_cnEpSin(fhd, problem_size, bounds, pop_size, M)

tic;
%problem_size = 50;

%%% change freq
freq_inti = 0.5;

max_nfes = M * pop_size;

rand('seed', sum(100 * clock));

val_2_reach = 10^(-8);
max_region = bounds(1,2);
min_region = bounds(1,1);
lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];
%fhd=@cec17_func;
pb = 0.4;
ps = .5;

S.Ndim = problem_size;
S.Lband = ones(1, S.Ndim)*bounds(1,1);
S.Uband = ones(1, S.Ndim)*bounds(1,2);

%%%% Count the number of maximum generations before as NP is dynamically
%%%% decreased
% G_Max = 0;
% if problem_size == 10
%     G_Max = 2163;
% end
% if problem_size == 30
%     G_Max = 2745;
% end
% if problem_size == 50
%     G_Max = 3022;
% end
% if problem_size == 100
%     G_Max = 3401;
% end
G_Max = M;

% num_prbs = 30;
% runs = 51;
run_funcvals = [];
% RecordFEsFactor = ...
%     [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, ...
%     0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
% progress = numel(RecordFEsFactor);

% allerrorvals = zeros(progress, runs, num_prbs);
% result=zeros(num_prbs,5);
% 
% fprintf('Running LSHADE_cnEpSin on D= %d\n', problem_size)
% for func = 1 : num_prbs
%     optimum = func * 100.0;
%     S.FuncNo = func;
%     
%     %% Record the best results
%     outcome = [];
%     
%     fprintf('\n-------------------------------------------------------\n')
%     fprintf('Function = %d, Dimension size = %d\n', func, problem_size)
%     
%     parfor run_id = 1 : runs
        
        run_funcvals = [];
        col=1;              %% to print in the first column in all_results.mat
        
        %%  parameter settings for L-SHADE
        p_best_rate = 0.11;    %0.11
        arc_rate = 1.4;
        memory_size = 5;
        %pop_size = 18 * problem_size;   %18*D
        SEL = round(ps*pop_size);
        
        max_pop_size = pop_size;
        min_pop_size = 4.0;
        
        nfes = 0;
        %% Initialize the main population
        popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
        pop = popold; % the old population becomes the current population
        
        fitness = [];
        for brojac = 1:size(pop,1)
            fitness(brojac) = feval(fhd, pop(brojac,:));
        end
%         fitness = feval(fhd,pop',func);
         fitness = fitness';
        
        bsf_fit_var = 1e+30;
        bsf_index = 0;
        bsf_solution = zeros(1, problem_size);
        
        %%%%%%%%%%%%%%%%%%%%%%%% for out
        for i = 1 : pop_size
            nfes = nfes + 1;
            
            if fitness(i) < bsf_fit_var
                bsf_fit_var = fitness(i);
                bsf_solution = pop(i, :);
                bsf_index = i;
            end
            
            if nfes > max_nfes; break; end
        end
        %%%%%%%%%%%%%%%%%%%%%%%% for out
        
        memory_sf = 0.5 .* ones(memory_size, 1);
        memory_cr = 0.5 .* ones(memory_size, 1);
        
        memory_freq = freq_inti*ones(memory_size, 1);
        memory_pos = 1;
        
        archive.NP = arc_rate * pop_size; % the maximum size of the archive
        archive.pop = zeros(0, problem_size); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions
        
        %% main loop
        gg=0;  %%% generation counter used For Sin
        igen =1;  %%% generation counter used For LS
        
        flag1 = false;
        flag2 = false;
        
        goodF1all = [];
        goodF2all =[];
        badF1all = [];
        badF2all = [];
        goodF1 = [];
        goodF2 = [];
        badF1 = [];
        badF2 = [];
        
        while nfes < max_nfes
            gg=gg+1;
            
            pop = popold; % the old population becomes the current population
            [temp_fit, sorted_index] = sort(fitness, 'ascend');
            
            mem_rand_index = ceil(memory_size * rand(pop_size, 1));
            mu_sf = memory_sf(mem_rand_index);
            mu_cr = memory_cr(mem_rand_index);
            mu_freq = memory_freq(mem_rand_index);
            
            %% for generating crossover rate
            cr = normrnd(mu_cr, 0.1);
            term_pos = find(mu_cr == -1);
            cr(term_pos) = 0;
            cr = min(cr, 1);
            cr = max(cr, 0);
            
            %% for generating scaling factor
            sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
            pos = find(sf <= 0);
            
            while ~ isempty(pos)
                sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
                pos = find(sf <= 0);
            end
            
            
            freq = mu_freq + 0.1 * tan(pi*(rand(pop_size, 1) - 0.5));
            pos_f = find(freq <=0);
            while ~ isempty(pos_f)
                freq(pos_f) = mu_freq(pos_f) + 0.1 * tan(pi * (rand(length(pos_f), 1) - 0.5));
                pos_f = find(freq <= 0);
            end
            
            sf = min(sf, 1);
            freq = min(freq, 1);
            
            LP = 20;
            flag1 = false;
            flag2 = false;
            if(nfes <= max_nfes/2)
                flag1 = false;
                flag2 = false;
                if (gg <= LP)
                    %% Both have the same probability
                    %% Those generations are the learning period
                    %% Choose one of them randomly
                    p1 = 0.5;
                    p2 = 0.5;
                    c=rand;
                    if(c < p1)
                        sf = 0.5.*( sin(2.*pi.*freq_inti.*gg+pi) .* ((G_Max-gg)/G_Max) + 1 ) .* ones(pop_size,problem_size);
                        flag1 = true;
                    else
                        sf = 0.5 *( sin(2*pi .* freq(:, ones(1, problem_size)) .* gg) .* (gg/G_Max) + 1 ) .* ones(pop_size,problem_size);
                        flag2 = true;
                    end
                    
                else
                    %% compute the probability as used in SaDE 
                    ns1 = size(goodF1,1);
                    ns1_sum = 0;
                    nf1_sum = 0;
                    %               for hh = 1 : size(goodF1all,2)
                    for hh = gg-LP : gg-1
                        ns1_sum = ns1_sum + goodF1all(1,hh);
                        nf1_sum = nf1_sum + badF1all(1,hh);
                    end
                    sumS1 = (ns1_sum/(ns1_sum + nf1_sum)) + 0.01;
                    
                    
                    ns2 = size(goodF2,1);
                    ns2_sum = 0;
                    nf2_sum = 0;
                    %             for hh = gg-LP : gg-1
                    %               for hh = 1 : size(goodF2all,2)
                    for hh = gg-LP : gg-1
                        ns2_sum = ns2_sum + goodF2all(1,hh);
                        nf2_sum = nf2_sum + badF2all(1,hh);
                    end
                    sumS2 = (ns2_sum/(ns2_sum + nf2_sum)) + 0.01;
                    
                    p1 = sumS1/(sumS1 + sumS2);
                    p2 = sumS2/(sumS2 + sumS1);
                    
                    if(p1 > p2)
                        sf = 0.5.*( sin(2.*pi.*freq_inti.*gg+pi) .* ((G_Max-gg)/G_Max) + 1 ) .* ones(pop_size,problem_size);
                        flag1 = true;
                        %                   size(goodF1,1)
                    else
                        sf = 0.5 *( sin(2*pi .* freq(:, ones(1, problem_size)) .* gg) .* (gg/G_Max) + 1 ) .* ones(pop_size,problem_size);
                        flag2 = true;
                        %                   size(goodF2,1)
                    end
                end
            end
            
            r0 = [1 : pop_size];
            popAll = [pop; archive.pop];
            [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
            
            pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
            randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions
            
            vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
            vi = boundConstraint(vi, pop, lu);
            
            %       %% Bin Crx
            %       mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
            %       rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
            %       jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
            %       ui = vi; ui(mask) = pop(mask);
            %       %%
            
            %% Bin crossover according to the Eigen coordinate system
            J_= mod(floor(rand(pop_size, 1)*problem_size), problem_size) + 1;
            J = (J_-1)*pop_size + (1:pop_size)';
            crs = rand(pop_size, problem_size) < cr(:, ones(1, problem_size));
            if rand<pb
                %% coordinate ratation
                
                %%%%% Choose neighbourhood region to the best individual
                best = pop(sorted_index(1), :);
                Dis = pdist2(pop,best,'euclidean'); % euclidean distance
                %D2 = sqrt(sum((pop(1,:) - best).^2, 2));
                
                %%%% Sort
                [Dis_ordered idx_ordered] = sort(Dis, 'ascend');
                Neighbour_best_pool = pop(idx_ordered(1:SEL), :); %%% including best also so start from 1
                Xsel = Neighbour_best_pool;
                %            sizz = size(Xsel)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %Xsel = pop(sorted_index(1:SEL), :);
                xmean = mean(Xsel);
                % covariance matrix calculation
                C =  1/(SEL-1)*(Xsel - xmean(ones(SEL,1), :))'*(Xsel - xmean(ones(SEL,1), :));
                C = triu(C) + transpose(triu(C,1)); % enforce symmetry
                [R,D] = eig(C);
                % limit condition of C to 1e20 + 1
                if max(diag(D)) > 1e20*min(diag(D))
                    tmp = max(diag(D))/1e20 - min(diag(D));
                    C = C + tmp*eye(problem_size);
                    [R, D] = eig(C);
                end
                TM = R;
                TM_=R';
                Xr = pop*TM;
                vi = vi*TM;
                %% crossover according to the Eigen coordinate system
                Ur = Xr;
                Ur(J) = vi(J);
                Ur(crs) = vi(crs);
                %%
                ui = Ur*TM_;
                
            else
                
                ui = pop;
                ui(J) = vi(J);
                ui(crs) = vi(crs);
                
            end
            %%%%%%%%
            
            children_fitness = [];
            for brojac = 1:size(ui,1)
                children_fitness(brojac) = feval(fhd, ui(brojac,:));
            end
            
%             children_fitness = feval(fhd, ui', func);
             children_fitness = children_fitness';
            
            
            %%%% To check stagnation
            flag = false;
            bsf_fit_var_old = bsf_fit_var;
            %%%%%%%%%%%%%%%%%%%%%%%% for out
            for i = 1 : pop_size
                nfes = nfes + 1;
                
                if children_fitness(i) < bsf_fit_var
                    bsf_fit_var = children_fitness(i);
                    bsf_solution = ui(i, :);
                    bsf_index = i;
                end
                
                if nfes > max_nfes; break; end
            end
            %%%%%%%%%%%%%%%%%%%%%%%% for out
            
            dif = abs(fitness - children_fitness);
            
            
            %% I == 1: the parent is better; I == 2: the offspring is better
            I = (fitness > children_fitness);
            goodCR = cr(I == 1);
            goodF = sf(I == 1);
            goodFreq = freq(I == 1);
            dif_val = dif(I == 1);
            
            %% chnage here also
            %% recored bad too
            badF = sf(I == 0);
            
            %% Change Noor
            if flag1 == true
                goodF1 = goodF;
                goodF1all = [goodF1all size(goodF1,1)];
                
                badF1 = badF;
                badF1all = [badF1all size(badF1,1)];
                
                %% Add zero for other one  or add 1 to prevent the case of having NaN
                goodF2all = [goodF2all 1];
                badF2all = [badF2all 1];
                
            end
            if flag2 == true
                goodF2 = goodF;
                goodF2all = [goodF2all size(goodF2,1)];
                
                badF2 = badF;
                badF2all = [badF2all size(badF2,1)];
                
                %% Add zero for other one
                goodF1all = [goodF1all 1];
                badF1all = [badF1all 1];
            end
            %%%%%%
            
            
            %      isempty(popold(I == 1, :))
            archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));
            
            [fitness, I] = min([fitness, children_fitness], [], 2);
            
            run_funcvals = [run_funcvals; fitness];
            
            popold = pop;
            popold(I == 2, :) = ui(I == 2, :);
            
            num_success_params = numel(goodCR);
            
            if num_success_params > 0
                sum_dif = sum(dif_val);
                dif_val = dif_val / sum_dif;
                
                %% for updating the memory of scaling factor
                memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
                
                %% for updating the memory of crossover rate
                if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
                    memory_cr(memory_pos)  = -1;
                else
                    memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
                end
                
                %% for updating the memory of freq
                if max(goodFreq) == 0 || memory_freq(memory_pos)  == -1
                    memory_freq(memory_pos)  = -1;
                else
                    memory_freq(memory_pos) = (dif_val' * (goodFreq .^ 2)) / (dif_val' * goodFreq);
                end
                
                memory_pos = memory_pos + 1;
                if memory_pos > memory_size;  memory_pos = 1; end
            end
            
            %% for resizing the population size
            plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);
            
            if pop_size > plan_pop_size
                reduction_ind_num = pop_size - plan_pop_size;
                if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end
                
                pop_size = pop_size - reduction_ind_num;
                SEL = round(ps*pop_size);
                for r = 1 : reduction_ind_num
                    [valBest indBest] = sort(fitness, 'ascend');
                    worst_ind = indBest(end);
                    popold(worst_ind,:) = [];
                    pop(worst_ind,:) = [];
                    fitness(worst_ind,:) = [];
                end
                
                archive.NP = round(arc_rate * pop_size);
                
                if size(archive.pop, 1) > archive.NP
                    rndpos = randperm(size(archive.pop, 1));
                    rndpos = rndpos(1 : archive.NP);
                    archive.pop = archive.pop(rndpos, :);
                end
            end
            
        end %%%%%%%%nfes
run_funcvals = run_funcvals(round(linspace(1,length(run_funcvals),M)));
runtime = toc;
        
%         bsf_error_val = bsf_fit_var - optimum;
%         if bsf_error_val < val_2_reach
%             bsf_error_val = 0;
%         end
%         
%         fprintf('%d th run, best-so-far error value = %1.8e\n', run_id , bsf_error_val)
%         outcome = [outcome bsf_error_val];
        
        %%% Noor: To print files
        %     errorVals= [];
        %     for w = 1 : progress
        %         bestold = run_funcvals(RecordFEsFactor(w) * max_nfes) - optimum;
        %         if abs(bestold)>1e-8
        %             errorVals(w)= abs(bestold);
        % %             col=col+1;
        %         else
        %             bestold=0;
        %             col=col+1;
        %             errorVals(w)= bestold;
        % %             current_eval=max_eval;
        %         end
        %     end
        %     allerrorvals(:, run_id, func) = errorVals;
        
%     end %% end 1 run
%     
%     fprintf('\n')
%     fprintf('min error value = %1.8e, max = %1.8e, median = %1.8e, mean = %1.8e, std = %1.8e\n', min(outcome), max(outcome), median(outcome), mean(outcome), std(outcome))
%     
%     result(func,1)=  min(outcome);
%     result(func,2)=  max(outcome);
%     result(func,3)=  median(outcome);
%     result(func,4)=  mean(outcome);
%     result(func,5)=  std(outcome);
%     
% end %% end 1 function run
% 
% disp(result);
% 
% name1 = 'results_stat_';
% name2 = num2str(problem_size);
% name3 = '.txt';
% f_name=strcat(name1,name2,name3);
% 
% save(f_name, 'result', '-ascii');

% %%% To print files
% for i =1 : num_prbs
%     name1 = 'LSHADE_cnEpSin';
%     name2 = num2str(i);
%     name3 = '_';
%     name4 = num2str(problem_size);
%     name5 = '.txt';
%     f_name=strcat(name1,name2,name3,name4,name5);
%     res = allerrorvals(:,:,i);
%     save(f_name, 'res', '-ascii');
% end




function [r1, r2] = gnR1R2(NP1, NP2, r0)

% gnA1A2 generate two column vectors r1 and r2 of size NP1 & NP2, respectively
%    r1's elements are choosen from {1, 2, ..., NP1} & r1(i) ~= r0(i)
%    r2's elements are choosen from {1, 2, ..., NP2} & r2(i) ~= r1(i) & r2(i) ~= r0(i)
%
% Call:
%    [r1 r2 ...] = gnA1A2(NP1)   % r0 is set to be (1:NP1)'
%    [r1 r2 ...] = gnA1A2(NP1, r0) % r0 should be of length NP1
%
% Version: 2.1  Date: 2008/07/01
% Written by Jingqiao Zhang (jingqiao@gmail.com)

NP0 = length(r0);

r1 = floor(rand(1, NP0) * NP1) + 1;
%for i = 1 : inf
for i = 1 : 99999999
    pos = (r1 == r0);
    if sum(pos) == 0
        break;
    else % regenerate r1 if it is equal to r0
        r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
    end
    if i > 1000, % this has never happened so far
        error('Can not genrate r1 in 1000 iterations');
    end
end

r2 = floor(rand(1, NP0) * NP2) + 1;
%for i = 1 : inf
for i = 1 : 99999999
    pos = ((r2 == r1) | (r2 == r0));
    if sum(pos)==0
        break;
    else % regenerate r2 if it is equal to r0 or r1
        r2(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
    end
    if i > 1000, % this has never happened so far
        error('Can not genrate r2 in 1000 iterations');
    end
end





%This function is used for L-SHADE bound checking 
function vi = boundConstraint (vi, pop, lu)

% if the boundary constraint is violated, set the value to be the middle
% of the previous value and the bound
%

[NP, D] = size(pop);  % the population size and the problem's dimension

%% check the lower bound
xl = repmat(lu(1, :), NP, 1);

pos = vi < xl;
vi(pos) = (pop(pos) + xl(pos)) / 2;

%% check the upper bound
xu = repmat(lu(2, :), NP, 1);
pos = vi > xu;
vi(pos) = (pop(pos) + xu(pos)) / 2;





function archive = updateArchive(archive, pop, funvalue)
% Update the archive with input solutions
%   Step 1: Add new solution to the archive
%   Step 2: Remove duplicate elements
%   Step 3: If necessary, randomly remove some solutions to maintain the archive size
%
% Version: 1.1   Date: 2008/04/02
% Written by Jingqiao Zhang (jingqiao@gmail.com)

if archive.NP == 0, return; end

if size(pop, 1) ~= size(funvalue,1), error('check it'); end

% Method 2: Remove duplicate elements
popAll = [archive.pop; pop ];
funvalues = [archive.funvalues; funvalue ];
[dummy IX]= unique(popAll, 'rows');
if length(IX) < size(popAll, 1) % There exist some duplicate solutions
  popAll = popAll(IX, :);
  funvalues = funvalues(IX, :);
end

if size(popAll, 1) <= archive.NP   % add all new individuals
  archive.pop = popAll;
  archive.funvalues = funvalues;
else                % randomly remove some solutions
  rndpos = randperm(size(popAll, 1)); % equivelent to "randperm";
  rndpos = rndpos(1 : archive.NP);
  
  archive.pop = popAll  (rndpos, :);
  archive.funvalues = funvalues(rndpos, :);
end
