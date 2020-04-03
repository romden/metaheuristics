clear all
clear
clc

addpath(genpath('fun'));


%% multiple runs

global Problem knapsack 

typeOfAlgorithm = 'neuroevolutionary'; % 'evolutionary'; % 
params.variation = 'ES'; % OnePointCrossover TwoPointCrossover UniformCrossover

objs = 2;%3:5;
probs = 250;%[500 1000 2000 5000 10000]; % [10000 5000 2000 1000 500]; %
seeds =  1000:1000:21000; %  21000:-1000:1000; %

for numKnapsacks = objs % dimensions
    
    for numItems = probs % problems
        
        % define knapsack
        inpFolder = sprintf('input/knapsack_%d_%d', numKnapsacks, numItems);        
        knapsack = struct('profit', [], 'weight', [], 'capacity', [], 'indexes', []);        
        knapsack.profit = importdata(sprintf('%s/%s', inpFolder, 'profit'));
        knapsack.weight = importdata(sprintf('%s/%s', inpFolder, 'weight'));
        knapsack.capacity = importdata(sprintf('%s/%s', inpFolder, 'capacity'));
        [~, knapsack.indexes] = sort( max(knapsack.profit./knapsack.weight, [], 2) );
        
        % define problem
        Problem.ObjFunction = 'MOKSmatlab';
        Problem.m = numKnapsacks;
        
        switch typeOfAlgorithm
            case 'evolutionary'
                Problem.n = numItems; % number of variables			
                
            case 'neuroevolutionary'
                knapsack.input = [knapsack.profit./repmat(max(knapsack.profit),numItems,1) knapsack.weight./repmat(max(knapsack.weight),numItems,1)];
                
                Problem.topology = [size(knapsack.input,2) 20 1];
                Problem.n = Problem.topology(2);
                for l = 2:numel(Problem.topology)
                    Problem.n = Problem.n + (Problem.topology(l-1) + 1) * Problem.topology(l);
                end
                Problem.lb = -1*ones(Problem.n, 1);
                Problem.ub = +1*ones(Problem.n, 1);                
        end

        
        for seed = seeds % runs
            
            % params for algorithm
            params.seed = seed;
            params.savePop = true;
            
            % save results
            folder = sprintf('results/%s_%s/knapsack_%d_%d/seed_%d', typeOfAlgorithm, params.variation, numKnapsacks, numItems, seed);
            if exist(folder,'dir') == 7;
                continue
            else
                mkdir(folder);
            end
            
            % run moea
            tic
            SMSEMOA(params);
            time = toc;
            
            for g = 0:100:1000
                SOURCE = sprintf('pop_%d.mat', g);
                if exist(SOURCE, 'file') > 0                    
                    DESTINATION = sprintf('%s/%s', folder, SOURCE);
                    movefile(SOURCE, DESTINATION, 'f');
                end
            end
            
            % output information
            fprintf('numKnapsacks: %d\t numItems: %d\t seed: %d\t Time: %.2f \n', numKnapsacks, numItems, seed, time);
            
        end
    end
    
end
