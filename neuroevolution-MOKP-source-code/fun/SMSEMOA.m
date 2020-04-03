function SMSEMOA(params, varargin)
%% Reference
% N. Beumea, B. Naujoksa and M. Emmerich.
% SMS-EMOA: Multiobjective selection based on dominated hypervolume.
% European Journal of Operational Research, 181(3): 1653-1669, 2007.
%% Input
% Problem.m - number of objective functions
% Problem.n - number of variables
% Problem.lb - lower bounds of variables (column vector)
% Problem.ub - upper bounds of variables (column vector)
% Problem.objFunction - objective function (column vector is returned)
%% Output
% A.x - approximation in the decision space
% A.f - approximation in the objective space
% Statistics - structure with running data
%% Implementation
% Roman Denysiuk

%% Initialization Procedure

global Problem Algorithm Population Offspring

% structure of algorithm parameters
Algorithm = struct('seed', 1000, ... seed for random number generator
                   'mu', 100, ... Population size
                   'maxEval', 100*1000, ... stopping criterion
                   'variation', 'UniformCrossover', ... type of variation operator (GA ES DE CMA) OnePointCrossover TwoPointCrossover UniformCrossover
                   'plotPop', false, ... population vizualization using plot (true, false)
                   'savePop', false ... save population during evolution (true, false)
                   );
               
% individual in the Population
Individual = struct('x', [], ... real chromosome                    
                    'f', [], ... objective vector
                    'rank', 0, ... individual's rank
                    'cv', 0, ... constraint violation
                    's', [] ... binary string 
                    );
                
if isfield(Problem, 'topology')
  Individual.topology = Problem.topology; % network topology
end

% initialize Population structure
Population = repmat(Individual, Algorithm.mu, 1);

% input params
initial = [];
if nargin > 0 && ~isempty(params)
    for str = {'seed' 'mu' 'maxEval' 'variation' 'plotPop' 'savePop'}
        if isfield(params, str{1})
            Algorithm.(str{1}) = params.(str{1});
        end
    end
    if isfield(params, 'popFile')
        initial = load(params.popFile);
    end
end

RandStream.setGlobalStream( RandStream('mt19937ar', 'Seed', Algorithm.seed) );
% RandStream.setDefaultStream( RandStream('mt19937ar', 'Seed', params.seed) ); % older version

% initialize strategy parameters
global cmaparams
if strcmp(Algorithm.variation, 'CMA')
    n = Problem.n;
    d = 1+n/2;
    p_target = 1/(5+sqrt(1/2));
    c_p = p_target/(2+p_target);
    c_c = 2/(n+2);
    c_cov = 2/(n^2+6);
    p_thresh = 0.44;
    cmaparams = struct('d', d, 'p_target', p_target, 'c_p', c_p, 'c_c', c_c, 'c_cov', c_cov, 'p_thresh', p_thresh);
    for i = 1:Algorithm.mu
        Population(i).p_succ = p_target;
        Population(i).sigma = 1;
        Population(i).p_c = 0;
        Population(i).C = eye(n);
    end
elseif strcmp(Algorithm.variation, 'ES')
    % init strategy params
    for i = 1:Algorithm.mu        
        Population(i).sigma = sqrt((Problem.ub-Problem.lb).^2/(12*Problem.n));
    end
end

% randomly generate initial population
for i = 1:Algorithm.mu
    
    if isempty(initial)
        
        if ~isfield(Problem, 'lb') && ~isfield(Problem, 'ub')
            Population(i).x = round(rand(Problem.n,1));
        else
            Population(i).x = Problem.lb + rand(Problem.n,1).*(Problem.ub - Problem.lb);
        end
        
        [Population(i)] = Evaluation(Population(i), varargin{:});
    else
        %Population(i).topology = Problem.topology;
        %Population(i).x = initial.Population(i).x;
        %Population(i).f = initial.Population(i).f;
        Population(i) = initial.Population(i);
        %[Population(i)] = Evaluation(Population(i), varargin{:});
    end
    
end
clear initial Individual

% vizualize current population
if (Algorithm.savePop)
    SavePopulation(0);
end
if (Algorithm.plotPop)
    ptr_plot = PlotPopulation([Population.f], [], 0, varargin{:});
end

% start evolutionary process
for i = 1:Algorithm.maxEval
    
    % generate offspring
    eval( sprintf('Variation%s(varargin{:});', Algorithm.variation) );
    
    % offspring is evaluated and added to the population   
    Population(end+1) = Evaluation(Offspring, varargin{:});
    
    % perform environmental selection
    EnvironmentalSelection(varargin{:});
    
    % vizualize current population
    if mod(i, 100*Algorithm.mu) == 0
        g = i/Algorithm.mu;
        if (Algorithm.savePop)
            SavePopulation(g);
        end
        if Algorithm.plotPop
            ptr_plot = PlotPopulation([Population.f], ptr_plot, g, varargin{:});
        end
    end
    
end

% SavePopulation('final');

return


%% Visualization Procedure
function ptr_plot = PlotPopulation(data, ptr_plot, g, varargin)

if isrow(data);
    data = data(:);
end

if isempty(ptr_plot)
    title(sprintf('Population at generation: %d', g));
    xlabel('f_1(x)');
    ylabel('f_2(x)');
    grid on;
    drawnow;
    hold on;
    ptr_plot = plot(data(1,:), data(2,:),'.');
else
    pause(0.01)
    title(sprintf('Population at generation: %d', g));
    set(ptr_plot, 'XData', data(1,:), 'YData', data(2,:), 'color', 'r');
    drawnow;
end

return

function SavePopulation(g)

global Population

% sort population
[~, idx] = sort([Population.f],2);
Population = Population(idx(1,:));

% save population
if isnumeric(g)
    save(sprintf('pop_%d', g), 'Population');
else
    save(sprintf('pop_%s', g), 'Population');
end

return
        
        
%% Evaluation Procedure
function [individual] = Evaluation(individual, varargin)

global Problem

% calc objectives
individual = feval(Problem.ObjFunction, Reparation(individual), varargin{:});

% if strcmp(Problem.ObjFunction, 'MOKS')
%     [individual.f, individual.x] = feval(Problem.ObjFunction, individual.x, varargin{:});
% elseif strcmp(Problem.ObjFunction, 'MOKS2')
%     individual = Reparation(individual);
%     individual = MOKS2(individual);
% else
%     individual = feval(Problem.ObjFunction, Reparation(individual), varargin{:});
% end

persistent neval
if isempty(neval); neval = 0; end
neval = neval + 1;

% print to screen
% fprintf('feval: %d\t cv: %.3f\t f1: %.4f\t f2: %.4f\t time: %.3f\n', neval, individual.cv, individual.f(1), individual.f(2), t);

return

function [individual] = Reparation(individual)

global Problem

if isfield(individual, 'topology') % repair topology
    individual.x(end-individual.topology(2)+1:end) = Repair( individual.x(end-individual.topology(2)+1:end) );
elseif isfield(Problem, 'ub') && isfield(Problem, 'lb')
    individual.x = min(Problem.ub, max(Problem.lb, individual.x));
end

return


%% Environmental Selection Procedure
function EnvironmentalSelection(varargin)

global Problem Population

% if any(all(cell2mat({Population(1:end-1).f})'==repmat(Population(end).f', numel(Population)-1, 1), 2)) 
%     Population(end) = []; % remove offspring if pop contain individual with identical objectives
%     return
% end

[value, index] = max([Population.cv]); % check constraint violation
if value > 0
    % worst is one having largest cv
    toRemove = index(1);
else
    % select based on dominance and hv
    NondominatedSorting(varargin{:}); % perform nondominated sorting
    index = find([Population.rank]==max([Population.rank])); % seelct last front
    frontSize = numel(index);
    if frontSize > Problem.m
        [~, s] = hypervolumeContribution(Population(index), varargin{:});
        [~, idx] = min(s);
        toRemove = index(idx);
    else
        toRemove = index(randi(frontSize));
    end
end
Population(toRemove) = []; % remove worst

% parameters update for CMA
global parent
if isfield(Population, 'C') 
    % update step size    
    if toRemove ~= parent
        UpdateStepSize( double(toRemove <= numel(Population)) ); % success if true, offspring is added to pop
    end
end

return



%% Nondominated Sorting Procedure
function NondominatedSorting(varargin)

global Population

% get objectives
score = [Population.f]';

% Population size
popSize = size(score,1);

% Boolean to track if individuals are ranked
rankedIndiv = false(popSize,1);
% Initialization: rank of individuals (denotes which front)
nondominatedRank = inf(popSize,1);

numObj = size(score,2);
dominationMatrix = false(popSize);
% First test for domination is to check which points have better function
% values in at least one of the objectives
for count = 1:numObj
    dominationMatrix = dominationMatrix | bsxfun(@lt,score(:,count),score(:,count)');
end
% Now, check to see if those points that pass the first test, if they are
% at least as good as others in all the objectives
for count = 1:numObj
    dominationMatrix = dominationMatrix & bsxfun(@le,score(:,count),score(:,count)');
end
% We will do the test along the column
dominationMatrix = ~dominationMatrix;

rankCounter = 1;
while ~all(rankedIndiv) && nnz(isfinite(nondominatedRank)) <= popSize
    dominates = all(dominationMatrix);
    nondominatedRank(dominates) = rankCounter;
    rankCounter = rankCounter + 1;
    
    dominationMatrix(dominates,:) = true;
    
    dominationMatrix(dominates,dominates) = false;
    
    rankedIndiv(dominates) = true;
end

% rank of individuals in the Population
cnondominatedRank = num2cell(nondominatedRank,2);
[Population.rank] = deal(cnondominatedRank{:});

return


%% Distance Metric Procedure
function [front, s] = hypervolumeContribution(front, varargin)

% if numel(front) == 1
%     front.I = inf;
%     return
% end

% get objectives
f = [front.f]';

% get size
[K, M] = size(f);

% scale each objective to the interval [0,1]
f = (f - repmat(min(f),K,1))./repmat(max(f)-min(f),K,1);

% init s metric values
s = zeros(K,1);

if M == 2
    
    % two objective case
    [~, index] = sort(f(:,1));
    for i = 2:K-1
        s(index(i))=(f(index(i+1),1)-f(index(i),1))*(f(index(i-1),2)-f(index(i),2));
    end
    s([index(1) index(end)]) = 1;
    
% elseif M == 3
%     
%     % three objective case
%     a=sort(f(:,1));
%     b=sort(f(:,2));
%     r = [1.1 1.1 1.1]; %max(f)+1;
%     a=[a;r(1)];
%     b=[b;r(2)];
%     
%     best1_f3=r(3)*ones(K);
%     best2_f3=r(3)*ones(K);
%     
%     for i=1:K
%         for j=1:K
%             for k=1:K
%                 if f(k,1)<=a(i) && f(k,2)<=b(j)
%                     % s_k dominates cell (i,j) conc. f1, f2
%                     if f(k,3)<best1_f3(i,j)
%                         best2_f3(i,j)=best1_f3(i,j);
%                         best1_f3(i,j)=f(k,3);
%                     elseif f(k,3)<best2_f3(i,j) && f(k,3)>best1_f3(i,j)
%                         best2_f3(i,j)=f(k,3);
%                     end
%                 end
%             end
%         end
%     end
%     
%     for i=1:K
%         for j=1:K
%             ownerNumber=0;
%             owner=0;
%             for k=1:K
%                 if f(k,1)<=a(i) && f(k,2)<=b(j) && best1_f3(i,j)==f(k,3)
%                     % s_k dominates cell (i,j) conc. f1, f2
%                     ownerNumber=ownerNumber+1;
%                     owner=k;
%                 end
%             end
%             if ownerNumber==1
%                 % cell (i,j) is dominated disjoint
%                 s(owner)=s(owner)+(a(i+1)-a(i))*(b(j+1)-b(j))*(best2_f3(i,j)-best1_f3(i,j));
%             end
%         end
%     end
    
elseif M > 2 && M < 7
    % find extrem
    [~, extrem] = max(f);
    s(extrem) = inf;
    % set ref point
    r = 1.1*ones(1, M);
    % calc all hypervolume
    hv = Indicator_HV_exact_mex(f, r);
    % calc contributions
    idx = true(K, 1);
    for i = 1:K
        if s(i) > 0
            continue
        end
        idx(i) = false;
        s(i) = hv - Indicator_HV_exact_mex(f(idx,:), r);
        idx(i) = true;
    end
    
else
    s = hv_contributions_mex(f', [5; 0.0; 1; 1e6]); % approximate hv, using hype
    [~, extrem] = max(f);
    s(extrem) = 1;
end

% assign s metric values to Population members
I = num2cell(s);
[front.I] = deal(I{:});

return


%% uniform parent selection
function [parents] = UniformSelection(n, varargin)

global Algorithm

parents = zeros(1, n);
for i = 1:n
   p = randi(Algorithm.mu);
   while any(ismember(p, parents))
       p = randi(Algorithm.mu);
   end
   parents(i) = p;
end

return


%% genetic algorithm variation
function VariationGA(varargin)
                                
global Problem Population Offspring
                                
% randomly select parents
[parents] = UniformSelection(2);
p1 = parents(1);
p2 = parents(2);

% init offspring
Offspring = Population(p1);

% real croms
Offspring.x = SBXcrossover(Population(p1).x, Population(p2).x, Problem.lb, Problem.ub);

return

            
%% evolution strategy variation
function VariationES(varargin)
                                
global Problem Population Offspring
                                
% select parent
[parents] = UniformSelection(1);
p = parents(1);

% init offspring
Offspring = Population(p);

% strategy parameters mutation
Offspring.sigma = ESMutationParameters(Offspring.sigma);

% real chromosome mutation
Offspring.x = ESMutationChromosome(Offspring.x, Offspring.sigma, Problem.lb, Problem.ub);

return


%% differential evolution variation  
function VariationDE(varargin)

global Problem Population Offspring

% select parent
[parents] = UniformSelection(4);
p = parents(1);
p1 = parents(2);
p2 = parents(3);
p3 = parents(4);

% init offspring
Offspring = Population(p);

% real chromosome
Offspring.x = DEoperator(Population(p).x, Population(p1).x, Population(p2).x, Population(p3).x, Problem.lb, Problem.ub);

return


%% CMA based operator (possible bug with cma), do not remember
function VariationCMA(varargin)
                                
global Problem Population Offspring

global parent
                                
% randomly select parents
parent = randi(numel(Population));

% init offspring
Offspring = Population(parent);

% real croms
[B, ~] = eig(Offspring.C);
Offspring.x = Offspring.x + Offspring.sigma*( B*randn(Problem.n,1) );
% Offspring.x = min(Problem.ub, max(Problem.lb, Offspring.x));

return


function UpdateStepSize(succ)

global Population Offspring parent
global cmaparams 

% params
for str = {'d' 'p_target' 'c_p' 'c_c' 'c_cov' 'p_thresh'}
    eval( sprintf('%s = cmaparams.%s;', str{1}, str{1}) );
end

% step size
p_succ = Population(parent).p_succ;
sigma = Population(parent).sigma;

Population(parent).p_succ = (1-c_p)*p_succ + c_p*succ;
Population(parent).sigma = sigma * exp(1/d*(p_succ-p_target)/(1-p_target));

% covariance matrix
p_c = Population(parent).p_c;
C = Population(parent).C;

if p_succ < p_thresh
    x_step = Offspring.x - Population(parent).x;    
    Population(parent).p_c = (1-c_c)*p_c + sqrt(c_c*(2-c_c))*x_step/sigma;
    Population(parent).C = (1-c_cov)*C + c_cov*(p_c*p_c');
else
    Population(parent).p_c = (1-c_c)*p_c;
    Population(parent).C = (1-c_cov)*C + c_cov*( (p_c*p_c') + c_c*(2-c_c)*C );
end

return


%% One Point Crossover 
function VariationOnePointCrossover(varargin)

global Problem Population Offspring

% select parent
[parents] = UniformSelection(2);
p1 = parents(1);
p2 = parents(2);

% init offspring
Offspring = Population(p1);

% real chromosome
Offspring.x = OnePointCrossover(Population(p1).x, Population(p2).x);

return

%% Two Point Crossover 
function VariationTwoPointCrossover(varargin)

global Problem Population Offspring

% select parent
[parents] = UniformSelection(2);
p1 = parents(1);
p2 = parents(2);

% init offspring
Offspring = Population(p1);

% real chromosome
Offspring.x = TwoPointCrossover(Population(p1).x, Population(p2).x);

return

%% Uniform Crossover 
function VariationUniformCrossover(varargin)

global Problem Population Offspring

% select parent
[parents] = UniformSelection(2);
p1 = parents(1);
p2 = parents(2);

% init offspring
Offspring = Population(p1);

% real chromosome
Offspring.x = UniformCrossover(Population(p1).x, Population(p2).x);

return