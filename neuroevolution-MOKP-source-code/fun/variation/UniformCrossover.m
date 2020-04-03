function [child] = UniformCrossover(parent1, parent2, pm, varargin)

% params
n = numel(parent1); % chromosome length
pc = 0.8; % crossover probability
if nargin < 3
    pm = 2*1/n; % mutation probability
end

% decide whether to apply crossover
flag = false; 
if rand <= pc
    flag = true;    
end

% init child
child = parent1;

for j = 1:n
    % apply uniform crossover
    if flag
        if rand <= 0.5
            child(j) = parent2(j);
        end
    end
    
    % bit flip mutation
    if rand <= pm
        if islogical(child(j))
            child(j) = ~child(j);
        else
            child(j) = double(~child(j));
        end
    end
end

return