function [child] = TwoPointCrossover(parent1, parent2, pm, varargin)

% params
n = numel(parent1); % chromosome length
pc = 1; % crossover probability
if nargin < 3
    pm = 1/n; % mutation probability
end

% decide whether to apply crossover
flag = false; 
if rand <= pc
    flag = true;    
end

% generate crossover points
point1 = randi(n);
point2 = randi(n);
while point2 == point1
    point2 = randi(n);
end

if point2 < point1
    temp = point1;
    point1 = point2;
    point2 = temp;
end

% init child
child = parent1;

for j = 1:n    
    % apply crossover
    if flag
        if point1 <= j && j <= point2
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