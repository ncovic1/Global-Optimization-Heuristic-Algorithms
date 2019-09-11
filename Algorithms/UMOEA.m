function [bestx, bestold, Convergence_curve, com_time]= UMOEA(fun, dim, constraints, N, Max_iter)

Par.n_opr=3;  %% number of operators in MODE
Par.n=dim;     %% number of decision vriables

if Par.n==10
    Par.CS=50; %% cycle
elseif Par.n==30
    Par.CS=100; %% cycle
elseif Par.n==50
    Par.CS=150; %% cycle
else
    Par.CS=150; %% cycle
end
% Par.CS=1; %% cycle

Par.xmin= constraints(:,1)';
Par.xmax= constraints(:,2)';
Par.Max_FES=N*Max_iter;
Par.PopSize=N; %% population size
Par.MinPopSize=4;
Par.prob_ls=0.1;

%% printing the detailed results- this will increase the computational time
Par.Printing=0; %% 1 to print; 0 otherwise


iter=0;             %% current generation


stream = RandStream('mt19937ar','Seed',sum(100*clock));
seed_run=stream.Seed; %% to record seeds for further validation, if needed
RandStream.setGlobalStream(stream);

%% define variables
current_eval=0;             %% current fitness evaluations
PS1=Par.PopSize;            %% define PS1
PS2=4+floor(3*log(Par.n));  %% define PS2
PS1 = PS1-PS2;
Par.PopSize=PS1+PS2;        %% PS = PS1+PS2

%% ====================== Initalize x ==================================
x=repmat(Par.xmin,Par.PopSize,1)+repmat((Par.xmax-Par.xmin),Par.PopSize,1).*rand(Par.PopSize,Par.n);

%% calc. fit. and update FES
tic;
fitx = zeros(Par.PopSize,1);
for ii = 1:Par.PopSize
    fitx(ii) = fun(x(ii,:));    
end
current_eval =current_eval+Par.PopSize;
res_det= min(repmat(min(fitx),1,Par.PopSize), fitx); %% used to record the convergence

%% ====================== store the best ==================
[bestold, bes_l]=min(fitx);     bestx= x(bes_l,:);
Convergence_curve = bestold;
%% ================== fill in each MOEA ===================================
%% DE
EA_1= x(1:PS1,:);    EA_obj1= fitx(1:PS1);
%% ES
EA_2= x(PS1+1:size(x,1),:);    EA_obj2= fitx(PS1+1:size(x,1));
%% ================ define CMA-ES parameters ==============================
setting=[];fitness=[];bnd=[];
[setting,fitness,bnd]= init_cmaes_par(setting,fitness,bnd,EA_obj2,EA_2, Par.n, PS2, Par.xmin, Par.xmax);

%% ===== prob. of each DE variant
probDE1=1./Par.n_opr .* ones(1,Par.n_opr);
%% ===================== archive data ====================================
arch_rate=2.6;
archive.NP = arch_rate * PS1; % the maximum size of the archive
archive.pop = zeros(0, Par.n); % the solutions stored in te archive
archive.funvalues = zeros(0, 1); % the function value of the archived solutions
%% ==================== to adapt CR and F =================================
hist_pos=1;
memory_size=6;
archive_f= ones(1,memory_size).*0.5;
archive_Cr= ones(1,memory_size).*0.5;
%%
stop_con=0; avgFE=Par.Max_FES; InitPop=PS1;

cy=0;indx = 0; Probs=ones(1,2);

%% main loop
while stop_con==0
    
    iter=iter+1;
    cy=cy+1; % to control CS
    %  ================ determine the best MOEA ===========================
    if(cy==ceil(Par.CS+1))
        
        %%calc normalized qualit -- NQual
        qual(1) = EA_obj1(1);qual(2) = EA_obj2(1);
        norm_qual = qual./sum(qual);
        norm_qual=1-norm_qual; %% to satisfy the bigger is the better
        
        %%Normalized diversity
        D(1) = mean(pdist2(EA_1(2:PS1,:),EA_1(1,:)));
        D(2) = mean(pdist2(EA_2(2:PS2,:),EA_2(1,:)));
        norm_div= D./sum(D);
        
        %%Total Imp
        Probs=norm_qual+norm_div;
        %%Update Prob_MODE and Prob_CMAES
        Probs = max(0.1, min(0.9,Probs./sum(Probs)));
        
        [~,indx]=max(Probs);
        if Probs(1)==Probs(2)
            indx=0;%% no sharing of information
        end
        
    elseif cy==2*ceil(Par.CS)
        
        %% share information
        if indx==1
            list_ind = randperm(PS1);
            list_ind= list_ind(1:(min(PS2,PS1)));
            EA_2(1:size(list_ind,2),:)= EA_1(list_ind,:);
            EA_obj2(1:size(list_ind,2))= EA_obj1(list_ind);
            [setting,fitness,bnd]= init_cmaes_par(setting,fitness,bnd,EA_obj2,EA_2, Par.n, PS2, Par.xmin, Par.xmax);
            setting.sigma= setting.sigma*(1- (current_eval/Par.Max_FES));
        else
            if (min (EA_2(1,:)))> -100 && (max(EA_2(1,:)))<100 %% share best sol. in EA_2 if it is feasible
                EA_1(PS1,:)= EA_2(1,:);
                EA_obj1(PS1)= EA_obj2(1);
                [EA_obj1, ind]=sort(EA_obj1);
                EA_1=EA_1(ind,:);
            end
            
        end
        %% reset cy and Probs
        cy=1;   Probs=ones(1,2);
    end
    %% ====================== MOEA1 <-- SAMODE ============================
    if (current_eval<Par.Max_FES)
        if rand<Probs(1)
            
            %% =============================== LR of PS ===================================================
            UpdPopSize = round((((Par.MinPopSize - InitPop) / Par.Max_FES) * current_eval) + InitPop);
            if PS1 > UpdPopSize
                reduction_ind_num = PS1 - UpdPopSize;
                if PS1 - reduction_ind_num <  Par.MinPopSize;
                    reduction_ind_num = PS1 - Par.MinPopSize;
                end
                %% remove the worst ind.
                for r = 1 : reduction_ind_num
                    vv=PS1;
                    EA_1(vv,:)=[];
                    EA_obj1(vv)=[];
                    PS1 = PS1 - 1;
                end
                archive.NP = round(arch_rate * PS1);
                if size(archive.pop, 1) > archive.NP
                    rndpos = randperm(size(archive.pop, 1));
                    rndpos = rndpos(1 : archive.NP);
                    archive.pop = archive.pop(rndpos, :);
                end
            end
            
            %% apply MODE
            [EA_1, EA_obj1,probDE1,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,current_eval,res_det] = ...
                SAMO_DE(fun, EA_1, EA_obj1,probDE1,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,...
                Par.xmin, Par.xmax,  Par.n,  PS1,  current_eval,res_det,Par.Printing,Par.Max_FES);
        end
    end
    %% ====================== MOEA2 <-- SAMO-ES ======================
    if (current_eval<Par.Max_FES)
        if   rand<Probs(2)
            [ EA_2, EA_obj2, setting,bestold,bestx,bnd,fitness,current_eval,res_det] = ...
                SAMO_ES(fun, EA_2, EA_obj2, setting, iter,bestold,bestx,fitness,bnd,...
                Par.xmin,Par.xmax,Par.n,PS2,current_eval,res_det,Par.Printing,Par.Max_FES);
            
        end
    end
    %% ============================ LS ====================================
    if current_eval>0.75*Par.Max_FES
        if rand<Par.prob_ls
            old_fit_eva=current_eval;
            [bestx,bestold,current_eval,succ] = LS (fun, bestx,bestold,Par,current_eval,Par.Max_FES,Par.xmin,Par.xmax);
            if succ==1 %% if LS was successful
                EA_1(PS1,:)=bestx';
                EA_obj1(PS1)=bestold;
                [EA_obj1, sort_indx]=sort(EA_obj1);
                EA_1= EA_1(sort_indx,:);
                
                EA_2=repmat(EA_1(1,:), PS2, 1);
                [setting,fitness,bnd]= init_cmaes_par(setting,fitness,bnd,EA_obj2,EA_2, Par.n, PS2, Par.xmin, Par.xmax);
                setting.sigma=1e-05;
                EA_obj2(1:PS2)= EA_obj1(1);
                Par.prob_ls=0.1;
            else
                Par.prob_ls=0.01; %% set p_LS to a small value it  LS was not successful
            end
            %% record best fitness -- set Par.Printing==0 if not
            if Par.Printing==1
                res_det= [res_det repmat(bestold,1,(current_eval-old_fit_eva))];
            end
            
        end
    end
    
    Convergence_curve = [Convergence_curve, bestold];
    
    %% ====================== stopping criterion check ====================
    if (current_eval>=Par.Max_FES)
        stop_con=1;
    end
    
    %% =============================== Print ==============================
    %          fprintf('current_eval\t %d fitness\t %d \n', current_eval, abs(Par.f_optimal-bestold));
    if stop_con
        L = length(Convergence_curve);
        ind = round((1:Max_iter)*L/Max_iter);
        Convergence_curve = Convergence_curve(ind);
        com_time= toc;%cputime-start_time;
    end
end
end



function [setting,fitness,bnd]= init_cmaes_par(setting,fitness,bnd,EA_obj2, EA_2, n, n2, xmin, xmax)
% Initialize dynamic internal state parameters

setting.flgActiveCMA=0;
setting.flgresume=0;
setting.flgDiagonalOnly=0;
%% remember that EA_2 was uniformly generated between xmin and xmax. 
%% So, mean(EA_2) does not violate the initialization condition in the competition
setting.xmean = mean(EA_2);
setting.xmean=setting.xmean';
setting.insigma=0.3;
setting.sigma = setting.insigma;
setting.maxdx = Inf;
setting.mindx = 0;

setting.sigma = max(setting.insigma);              % overall standard deviation
setting.pc = zeros(n,1); setting.ps = zeros(n,1);  % evolution paths for setting.C and setting.sigma


setting.idx = (xmin > -Inf) | (xmax < Inf);
if length(setting.idx) == 1
    setting.idx = setting.idx * ones(n,1);
end
if length(setting.insigma) == 1
    setting.insigma = setting.insigma * ones(n,1) ;
end
setting.diagD = setting.insigma/max(setting.insigma);      % diagonal matrix D defines the scaling
setting.diagC = setting.diagD.^2;
if setting.flgDiagonalOnly ~= 1            % use at some point full covariance matrix
    setting.B = eye(n,n);                      % setting.B defines the coordinate system
    setting.BD = setting.B.*repmat(setting.diagD',n,1);        % setting.B*D for speed up only
    setting.C = diag(setting.diagC);                   % covariance matrix == setting.BD*(setting.BD)'
end

if setting.flgDiagonalOnly
    setting.B = 1;
end
  setting.D = ones(n,1); 
  setting.invsqrtC =   setting.B * diag(  setting.D.^-1) *   setting.B';

fitness.hist=NaN*ones(1,10+ceil(3*10*n/n2)); % history of fitness values
fitness.histsel=NaN*ones(1,10+ceil(3*10*n/n2)); % history of fitness values
fitness.histbest=[]; % history of fitness values
fitness.histmedian=[]; % history of fitness values

% Initialize boundary handling
bnd.isactive = any(xmin > -Inf) || any(xmax < Inf);
if bnd.isactive
    bnd.weights = zeros(n,1);         % setting.weights for bound penalty
    % scaling is better in axis-parallel case, worse in rotated
    bnd.flgscale = 0; % scaling will be omitted if zero
    if bnd.flgscale ~= 0
        bnd.scale = setting.diagC/mean(setting.diagC);
    else
        bnd.scale = ones(n,1);
    end
    setting.idx = (xmin > -Inf) | (xmax < Inf);
    if length(setting.idx) == 1
        setting.idx = setting.idx * ones(n,1);
    end
    setting.idx = (xmin > -Inf) | (xmax < Inf);
    if length(setting.idx) == 1
        setting.idx = setting.idx * ones(n,1);
    end
    bnd.isbounded = zeros(n,1);
    bnd.isbounded(setting.idx) = 1;
    setting.maxdx = min(setting.maxdx, (xmax - xmin)/2);
    if any(setting.sigma*sqrt(setting.diagC) > setting.maxdx')
         setting.fac = min(setting.maxdx ./ sqrt(setting.diagC))/setting.sigma;
        setting.sigma = min(setting.maxdx ./ sqrt(setting.diagC));
    end
    
    setting.dd = setting.diagC;
    bnd.dfithist = 1;              % delta fit for setting setting.weights
    bnd.aridxpoints = [];          % remember complete outside points
    bnd.arfitness = [];            % and their fitness
    bnd.validfitval = 0;
    bnd.iniphase = 1;
    % scaling is better in axis-parallel case, worse in rotated
    bnd.flgscale = 0; % scaling will be omitted if zero
    if bnd.flgscale ~= 0
        bnd.scale = setting.diagC/mean(setting.diagC);
    else
        bnd.scale = ones(n,1);
    end
end

% Initialize dynamic (internal) strategy parameters and constants
% usually close to 1
% Initialize dynamic (internal) strategy parameters and constants


setting.endinvsqrtC = 0;%setting.B * diag(D.^-1) * setting.B';    % setting.C^-1/2
setting.eigeneval = 0;                      % track update of setting.B and D
setting.chiN=n^0.5*(1-1/(4*n)+1/(21*n^2));  % expectation of
%   ||n(0,I)|| == norm(randn(n,1))

setting.mu = ceil(n2/2);               % number of parents/points for recombination
setting.weights = log(max(setting.mu, n/2) + 1/2)-log(1:setting.mu)'; % muXone array for weighted recombination setting.mu = floor(setting.mu);
setting.mueff=sum(setting.weights)^2/sum(setting.weights.^2); % variance-effective size of setting.mu
setting.weights = setting.weights/sum(setting.weights);     % normalize recombination setting.weights array

% Strategy parameter setting: Adaptation
setting.cc = (4 + setting.mueff/n) / (n+4 + 2*setting.mueff/n); % time constant for cumulation for setting.C
setting.cs = (setting.mueff+2) / (n+setting.mueff+3);  % t-const for cumulation for setting.sigma control
setting.ccov1 = 2 / ((n+1.3)^2+setting.mueff);    % learning rate for rank-one update of setting.C
setting.ccovmu = 2 * (setting.mueff-2+1/setting.mueff) / ((n+2)^2+setting.mueff);  % and for rank-setting.mu update
setting.damps = 0.5 + 0.5*min(1, (0.27*n2/setting.mueff-1)^2) + 2*max(0,sqrt((setting.mueff-1)/(n+1))-1) + setting.cs; % damping for setting.sigma

if setting.flgDiagonalOnly
    setting.ccov1_sep = min(1, setting.ccov1 * (n+1.5) / 3);
    setting.ccovmu_sep = min(1-setting.ccov1_sep, setting.ccovmu * (n+1.5) / 3);
else
    setting.ccov1_sep=0;
    setting.ccovmu_sep=0;
end

setting.stopOnEqualFunctionValues= 2 + n/3;
setting.arrEqualFunvals = zeros(1, 10+n);
if isempty(setting.insigma)
    if all(size(( setting.xmean)) > 1)
        setting.insigma = std( setting.xmean, 0, 2);
    else
    end
end


fitness.hist(1)= EA_obj2(1);%cec14_func( setting.xmean,I_fno);
fitness.histsel(1)=fitness.hist(1);

setting.xold =  setting.xmean;
end




function [x, fitx,prob,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,current_eval,res_det ] = ...
    SAMO_DE(fun, x, fitx,prob,bestold,bestx,archive,hist_pos,memory_size, archive_f,archive_Cr,xmin, xmax,  n,...
    PopSize,  current_eval,res_det,Printing,Max_FES)


vi=zeros(PopSize,n);

%% calc CR and F
mem_rand_index = ceil(memory_size * rand(PopSize, 1));
mu_sf = archive_f(mem_rand_index);
mu_cr = archive_Cr(mem_rand_index);

%% ========================= generate CR ==================================
cr = normrnd(mu_cr, 0.1)';
cr(mu_cr == -1) = 0;
cr = min(cr, 1);

%% ========================= generate F ===================================
F = mu_sf + 0.1 * tan(pi * (rand( 1,PopSize) - 0.5));
pos = find(F <= 0);
while ~ isempty(pos)
    F(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand( 1,length(pos)) - 0.5));
    pos = find(F <= 0);
end
F = min(F, 1)';

%% ======================== generate new x =================================
popAll = [x;archive.pop]; %% set archive
r0 = 1 : PopSize;
%% generate random integer numbers
[r1, r2,r3] = gnR1R2(PopSize, size(popAll, 1), r0);

%% mutation
bb= rand(PopSize, 1);
probiter = prob(1,:);
l2= sum(prob(1:2));
op_1 = bb <=  probiter(1)*ones(PopSize, 1);
op_2 = bb > probiter(1)*ones(PopSize, 1) &  bb <= (l2*ones(PopSize, 1)) ;
op_3 = bb > l2*ones(PopSize, 1) &  bb <= (ones(PopSize, 1)) ;

pNP = max(round(0.1 * PopSize), 2); %% choose at least two best solutions
randindex = ceil(rand(1, PopSize) .* pNP); %% select from [1, 2, 3, ..., pNP]
randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
phix = x(randindex, :); %% randomly choose one of the top 10% solutions
%% DE2 -- current-to-pbest/archive
vi(op_1==1,:) = x(op_1==1,:)+ F(op_1==1, ones(1, n)) .*(phix(op_1==1,:) - x(op_1==1,:) + x(r1(op_1==1), :) - popAll(r2(op_1==1), :));
%% DE2 --  current-to-pbest/without archive
vi(op_2==1,:) =  x(op_2==1,:)+ F(op_2==1, ones(1, n)) .*(phix(op_2==1,:) - x(op_2==1,:) + x(r1(op_2==1), :) - x(r3(op_2==1), :));
%% DE3 -- weighted-phi-rand
pNP = max(round(0.5 * PopSize), 2); %% choose at least two best solutions
randindex = ceil(rand(1, PopSize) .* pNP); %% select from [1, 2, 3, ..., pNP]
randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
phix = x(randindex, :); %% randomly choose one of the top 50% solutions
%% DE3 -- weighted DE
vi(op_3==1,:) = F(op_3==1, ones(1, n)).* x(r1(op_3==1), :) + phix(op_3==1,:) - x(r3(op_3==1), :);

%% handle boundaries
vi = han_boun(vi, xmax, xmin, x,PopSize,1);
%% crossover
mask = rand(PopSize, n) > cr(:, ones(1, n)); % mask is used to indicate which elements of ui comes from the parent
rows = (1 : PopSize)'; cols = floor(rand(PopSize, 1) * n)+1; % choose one position where the element of ui doesn't come from the parent
jrand = sub2ind([PopSize n], rows, cols); mask(jrand) = false;
ui = vi; ui(mask) = x(mask);

%% evaluate
fitx_new = zeros(PopSize,1);
for ii = 1:PopSize
    fitx_new(ii) = fun(ui(ii,:));
end
%% update FITNESS EVALUATIONS
current_eval =current_eval+PopSize;

%% calc. imprv. for Cr and F
diff = abs(fitx - fitx_new);
I =(fitx_new < fitx);
goodCR = cr(I == 1);
goodF = F(I == 1);

%% ========================= update archive ===============================
archive = updateArchive(archive, x(I == 1, :), fitx(I == 1)');
%% ==================== update Prob. of each DE ===========================
diff2 = max(0,(fitx - fitx_new))./abs(fitx);
count_S(1)=max(0,mean(diff2(op_1==1)));
count_S(2)=max(0,mean(diff2(op_2==1)));
count_S(3)=max(0,mean(diff2(op_3==1)));

%% update probs.
%%% Althouth count_S~=0 may slow down the convergence, it gives more
%%% diversity. In case you look for a better convergence you can set it to
%%% sum(count_S)~=0 
if count_S~=0 
prob= max(0.1,min(0.9,count_S./(sum(count_S))));
else
    prob=1/3 * ones(1,3);
end
%% ==================== update x and fitx =================================
fitx(I==1)= fitx_new(I==1);
x(I == 1, :) = ui(I == 1, :);
%% =================== update memory cr and F =============================
num_success_params = numel(goodCR);
if num_success_params > 0
    weightsDE = diff(I == 1)./ sum(diff(I == 1));
    %% for updating the memory of scaling factor
    archive_f(hist_pos) = (weightsDE' * (goodF .^ 2))./ (weightsDE' * goodF);
    
    %% for updating the memory of crossover rate
    if max(goodCR) == 0 || archive_Cr(hist_pos)  == -1
        archive_Cr(hist_pos)  = -1;
    else
        archive_Cr(hist_pos) = (weightsDE' * (goodCR .^ 2)) / (weightsDE' * goodCR);
    end
    
    hist_pos= hist_pos+1;
    if hist_pos > memory_size;  hist_pos = 1; end

end

%% sort new x, fitness
[fitx, ind]=sort(fitx);
x=x(ind,:);

%% record the best value after checking its feasiblity status
if fitx(1)<bestold  && min(x(ind(1),:))>=xmin(1) && max(x(ind(1),:))<=xmax(1)
    bestold=fitx(1);
     bestx= x(1,:);
end
%% check to print
if Printing==1
    res_det= [res_det repmat(bestold,1,PopSize)];
end

end



function [r1, r2,r3] = gnR1R2(NP1, NP2, r0)

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
    if i > 1000 % this has never happened so far
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
    if i > 1000 % this has never happened so far
        error('Can not genrate r2 in 1000 iterations');
    end
end

r3= floor(rand(1, NP0) * NP1) + 1;
%for i = 1 : inf
for i = 1 : 99999999
    pos = ((r3 == r0) | (r3 == r1) | (r3==r2));
    if sum(pos)==0
        break;
    else % regenerate r2 if it is equal to r0 or r1
         r3(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
    end
    if i > 1000 % this has never happened so far
        error('Can not genrate r2 in 1000 iterations');
    end
end
end


%% ============United Multi-Operator Evolutionary AlgorithmsII ============
% Should you have any queries, please contact
% Dr. Saber Elsayed. University of New South Wales at Canberra
% s.elsayed@adfa.edu.au
% www.saberelsayd.net or
% https://sites.google.com/site/saberelsayed3/home
% =========================================================================

function x = han_boun (x, xmax, xmin, x2, PopSize,hb)

switch hb
    case 1 % for DE
        x_L = repmat(xmin, PopSize, 1);
        pos = x < x_L;
        x(pos) = (x2(pos) + x_L(pos)) / 2;
        
        x_U = repmat(xmax, PopSize, 1);
        pos = x > x_U;
        x(pos) = (x2(pos) + x_U(pos)) / 2;
        
    case 2 % for CMA-ES
        x_L = repmat(xmin, PopSize, 1);
        pos = x < x_L;
        x_U = repmat(xmax, PopSize, 1);
        x(pos) = min(x_U(pos),max(x_L(pos),2*x_L(pos)-x2(pos)))  ;
        pos = x > x_U;
        x(pos) = max(x_L(pos),min(x_U(pos),2*x_L(pos)-x2(pos)));
        
end  
end


function archive = updateArchive(archive, pop, funvalue)
% Update the archive with input solutions
%   Step 1: Add new solution to the archive
%   Step 2: Remove duplicate elements
%   Step 3: If necessary, randomly remove some solutions to maintain the archive size
%
% Version: 1.1   Date: 2008/04/02
% Written by Jingqiao Zhang (jingqiao@gmail.com)

if archive.NP == 0, return; end

% Method 2: Remove duplicate elements
popAll = [archive.pop; pop ];
funvalues = [archive.funvalues; funvalue' ];
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
end



function[ x, fitx, setting,bestold,bestx,bnd,fitness,current_eval,res_det] = ...
    SAMO_ES(fun, x, ~, setting, iter,bestold,bestx,fitness,bnd,xmin,xmax,n,PopSize,current_eval,res_det,Printing,Max_FES)

stopOnWarnings=0;
noiseReevals = 0;
fitness.raw = NaN(1, PopSize + noiseReevals);
fitness.raw(PopSize + find(isnan(fitness.raw(1:noiseReevals)))) = NaN;

arz = randn(n,PopSize);
arx = repmat(setting.xmean, 1, PopSize) + setting.sigma * (setting.BD * arz);

%% ignore handling the boundaries constraints during the first 50% evolutionary process
%%-this is based on our earlier analysis carried out on UMOEAs in 2014.
handle_limit=0.5;
if current_eval >=handle_limit*Max_FES
    arxvalid =han_boun(arx', xmax, xmin, x,PopSize,2);
    arxvalid=arxvalid';
else
    arxvalid=arx;
end
%% evaluate and update cfe
fitness.raw = zeros(size(arxvalid,2),1);
for ii = 1:size(arxvalid,2)
    fitness.raw(ii) = fun(arxvalid(:,ii)');
end
current_eval=current_eval+PopSize; %% increase the fitness evaluations

fitness.sel= fitness.raw ;
[fitness.sel, fitness.idxsel] = sort(fitness.sel);

fitness.raw= fitness.raw(fitness.idxsel);
arxvalid= arxvalid(:, fitness.idxsel);
arx= arx(:, fitness.idxsel);
arz=arz(:, fitness.idxsel);
[~,pos_ro]=min(fitness.raw);

%% record the best value after checking its feasiblity status
if fitness.raw(pos_ro) < bestold && (min(arxvalid(:,pos_ro)))>=xmin(1) && (max(arxvalid(:,pos_ro)))<=xmax(1)
    bestold=fitness.raw(pos_ro);
    bestx= arxvalid(:,pos_ro)';
end
if Printing==1
    res_det= [res_det repmat(bestold,1,PopSize)];
end


% Calculate new setting.xmean, this is selection and recombination
setting.xold = setting.xmean; % for speed up of Eq. (2) and (3)
cmean =1;% 1/min(max((PopSize-1*n)/2, 1), n);  % == 1/kappa
setting.xmean = (1-cmean) * setting.xold + cmean * arx(:,(1:setting.mu))*setting.weights;
if  current_eval >=handle_limit*Max_FES
    % setting.xmean = xintobounds(setting.xmean, xmin', xmax');
    setting.xmean =han_boun(setting.xmean', xmax, xmin, x(1,:),1,2);
    setting.xmean=setting.xmean';
end
zmean = arz(:,(1:setting.mu))*setting.weights;%==D^-1*setting.B'*(setting.xmean-setting.xold)/setting.sigma
% Cumulation: update evolution paths
setting.ps = (1-setting.cs)*setting.ps + sqrt(setting.cs*(2-setting.cs)*setting.mueff) * (setting.B*zmean);          % Eq. (4)
hsig = norm(setting.ps)/sqrt(1-(1-setting.cs)^(2*iter))/setting.chiN < 1.4 + 2/(n+1);

setting.pc = (1-setting.cc)*setting.pc ...
    + hsig*(sqrt(setting.cc*(2-setting.cc)*setting.mueff)/setting.sigma/cmean) * (setting.xmean-setting.xold);     % Eq. (2)
if hsig == 0
    % disp([num2str(iter) ' ' num2str(counteval) ' setting.pc update stalled']);
end
% Adapt covariance matrix
neg.ccov = 0;  % TODO: move parameter setting upwards at some point
if setting.ccov1 + setting.ccovmu > 0                                                    % Eq. (3)
    if setting.flgDiagonalOnly % internal linear(?) complexity
        setting.diagC = (1-setting.ccov1_sep-setting.ccovmu_sep+(1-hsig)*setting.ccov1_sep*setting.cc*(2-setting.cc)) * setting.diagC ... % regard old matrix
            + setting.ccov1_sep * setting.pc.^2 ...               % plus rank one update
            + setting.ccovmu_sep ...                      % plus rank setting.mu update
            * (setting.diagC .* (arz(:,(1:setting.mu)).^2 * setting.weights));
        %             * (repmat(setting.diagC,1,setting.mu) .* arz(:,(1:setting.mu)).^2 * setting.weights);
        setting.diagD = sqrt(setting.diagC); % replaces eig(setting.C)
    else
        arpos = (arx(:,(1:setting.mu))-repmat(setting.xold,1,setting.mu)) / setting.sigma;
        setting.C = (1-setting.ccov1-setting.ccovmu+(1-hsig)*setting.ccov1*setting.cc*(2-setting.cc)) * setting.C ... % regard old matrix
            + setting.ccov1 * setting.pc*setting.pc' ...     % plus rank one update
            + setting.ccovmu ...             % plus rank setting.mu update
            * arpos * (repmat(setting.weights,1,n) .* arpos');
        % is now O(setting.mu*n^2 + setting.mu*n), was O(setting.mu*n^2 + setting.mu^2*n) when using diag(setting.weights)
        %   for setting.mu=30*n it is now 10 times faster, overall 3 times faster
        
        setting.diagC = diag(setting.C);
    end
end


% Adapt setting.sigma
setting.sigma = setting.sigma * exp(min(1, (sqrt(sum(setting.ps.^2))/setting.chiN - 1) * setting.cs/setting.damps));             % Eq. (5)
% disp([iter norm(setting.ps)/setting.chiN]);

if 11 < 3   % testing with optimal step-size
    setting.sigma = 0.04 * setting.mueff * sqrt(sum(setting.xmean.^2)) / n; % 20D,lam=1000:25e3
    setting.sigma = 0.3 * setting.mueff * sqrt(sum(setting.xmean.^2)) / n; % 20D,lam=(40,1000):17e3
    %      75e3 with def (1.5)
    %      35e3 with setting.damps=0.25
end
if 11 < 3
    
    setting.xmean = ones(n,1);
end

% Update setting.B and D from setting.C

if ~setting.flgDiagonalOnly && (setting.ccov1+setting.ccovmu+neg.ccov) > 0 && mod(iter, 1/(setting.ccov1+setting.ccovmu+neg.ccov)/n/10) < 1
    setting.C=triu(setting.C)+triu(setting.C,1)'; % enforce symmetry to prevent complex numbers
    [setting.B,tmp] = eig(setting.C);     % eigen decomposition, setting.B==normalized eigenvectors
    % effort: approx. 15*n matrix-vector multiplications
    setting.diagD = diag(tmp);
    
    % limit condition of setting.C to 1e14 + 1
    if min(setting.diagD) <= 0
        
        setting.diagD(setting.diagD<0) = 0;
        tmp = max(setting.diagD)/1e14;
        setting.C = setting.C + tmp*eye(n,n); setting.diagD = setting.diagD + tmp*ones(n,1);
        
    end
    if max(setting.diagD) > 1e14*min(setting.diagD)
        
        tmp = max(setting.diagD)/1e14 - min(setting.diagD);
        setting.C = setting.C + tmp*eye(n,n); setting.diagD = setting.diagD + tmp*ones(n,1);
        
    end
    
    setting.diagC = diag(setting.C);
    setting.diagD = sqrt(setting.diagD); % D contains standard deviations now
    % setting.diagD = setting.diagD / prod(setting.diagD)^(1/n);  setting.C = setting.C / prod(setting.diagD)^(2/n);
    setting.BD = setting.B.*repmat(setting.diagD',n,1); % O(n^2)
end % if mod

% Align/rescale order of magnitude of scales of setting.sigma and setting.C for nicer output
% TODO: interference with sigmafacup: replace 1e10 with 2*sigmafacup
% not a very usual case
if 1 < 2 && setting.sigma > 1e10*max(setting.diagD) && setting.sigma > 8e14 * max(setting.insigma)
    fac = setting.sigma; % / max(setting.diagD);
    setting.sigma = setting.sigma/fac;
    setting.pc = fac * setting.pc;
    setting.diagD = fac * setting.diagD;
    if ~setting.flgDiagonalOnly
        setting.C = fac^2 * setting.C; % disp(fac);
        setting.BD = setting.B .* repmat(setting.diagD',n,1); % O(n^2), but repmat might be inefficient todo?
    end
    setting.diagC = fac^2 * setting.diagC;
end

if setting.flgDiagonalOnly > 1 && iter > setting.flgDiagonalOnly
    % full covariance matrix from now on
    setting.flgDiagonalOnly = 0;
    setting.B = eye(n,n);
    setting.BD = diag(setting.diagD);
    setting.C = diag(setting.diagC); % is better, because correlations are spurious anyway
end

% ----- numerical error management -----
% Adjust maximal coordinate axis deviations
if any(setting.sigma*sqrt(setting.diagC) > setting.maxdx')
    setting.sigma = min(setting.maxdx ./ sqrt(setting.diagC'));
    %warning(['Iteration ' num2str(iter) ': coordinate axis std ' ...
    %         'deviation at upper limit of ' num2str(setting.maxdx)]);
    % stopflag(end+1) = {'maxcoorddev'};
end
% Adjust minimal coordinate axis deviations
if any(setting.sigma*sqrt(setting.diagC) < setting.mindx)
    setting.sigma = max(setting.mindx ./ sqrt(setting.diagC)) * exp(0.05+setting.cs/setting.damps);
    %warning(['Iteration ' num2str(iter) ': coordinate axis std ' ...
    %         'deviation at lower limit of ' num2str(setting.mindx)]);
    % stopflag(end+1) = {'mincoorddev'};;
end
% Adjust too low coordinate axis deviations
if any(setting.xmean == setting.xmean + 0.2*setting.sigma*sqrt(setting.diagC))
    if stopOnWarnings
        %         stopflag(end+1) = {'warnnoeffectcoord'};
    else
        %         warning(['Iteration ' num2str(iter) ': coordinate axis std ' ...
        %             'deviation too low' ]);
        if setting.flgDiagonalOnly
            setting.diagC = setting.diagC + (setting.ccov1_sep+setting.ccovmu_sep) * (setting.diagC .* ...
                (setting.xmean == setting.xmean + 0.2*setting.sigma*sqrt(setting.diagC)));
        else
            setting.C = setting.C + (setting.ccov1+setting.ccovmu) * diag(setting.diagC .* ...
                (setting.xmean == setting.xmean + 0.2*setting.sigma*sqrt(setting.diagC)));
        end
        setting.sigma = setting.sigma * exp(0.05+setting.cs/setting.damps);
    end
end
% Adjust step size in case of (numerical) precision problem
if setting.flgDiagonalOnly
    tmp = 0.1*setting.sigma*setting.diagD;
else
    tmp = 0.1*setting.sigma*setting.BD(:,1+floor(mod(iter,n)));
end
if all(setting.xmean == setting.xmean + tmp)
    %     ii = 1+floor(mod(iter,n));
    if stopOnWarnings
    else
        setting.sigma = setting.sigma * exp(0.2+setting.cs/setting.damps);
    end
end
% Adjust step size in case of equal function values (flat fitness)
% isequalfuncvalues = 0;
if fitness.sel(1) == fitness.sel(1+ceil(0.1+PopSize/4))
    % isequalfuncvalues = 1;
    if setting.stopOnEqualFunctionValues
        setting.arrEqualFunvals = [iter setting.arrEqualFunvals(1:end-1)];
        % stop if this happens in more than 33%
        if setting.arrEqualFunvals(end) > iter - 3 * length(setting.arrEqualFunvals)
            %             stopflag(end+1) = {'equalfunvals'};
        end
    else
        if flgWarnOnEqualFunctionValues
            %             warning(['Iteration ' num2str(iter) ...
            %                 ': equal function values f=' num2str(fitness.sel(1)) ...
            %                 ' at maximal main axis setting.sigma ' ...
            %                 num2str(setting.sigma*max(setting.diagD))]);
        end
        setting.sigma = setting.sigma * exp(0.2+setting.cs/setting.damps);
    end
end
% Adjust step size in case of equal function values
if iter > 2 && myrange([fitness.hist fitness.sel(1)]) == 0
    if stopOnWarnings
        % 	stopflag(end+1) = {'warnequalfunvalhist'};
    else
        %         warning(['Iteration ' num2str(iter) ...
        %             ': equal function values in history at maximal main ' ...
        %             'axis setting.sigma ' num2str(setting.sigma*max(setting.diagD))]);
        setting.sigma = setting.sigma * exp(0.2+setting.cs/setting.damps);
    end
end

%% print out final results
x= arxvalid';
fitx= fitness.raw;
end

function res=myrange(x)
res = max(x) - min(x);
end


function [x,f,current_eval,succ] = LS (fun, bestx,f,Par,current_eval,Max_FES,xmin,xmax)


Par.LS_FE=ceil(20.0000e-003*Max_FES ); %% Max FFEs_LS

options=optimset('Display','off','algorithm','interior-point',...
    'UseParallel','never','MaxFunEvals',Par.LS_FE) ;

[Xsqp, FUN , ~ , details]=fmincon(fun, bestx(1,:),[],[],[],[],xmin,xmax, [],options);

%% check if there is an improvement in the fitness value and update P_{ls}
if (f-FUN)>0
    succ=1;
    f = FUN;
    x(1,:)=Xsqp;
else
    succ=0;
    x=bestx;
   
end
%% update FFEs
current_eval=current_eval+details.funcCount;
end