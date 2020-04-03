function [individual] = MOKSmatlab(individual, varargin)

global knapsack

% define binary string
if isfield(individual, 'topology')
    x = SimulateIndividual(individual, knapsack.input) > 0.5;
else
    x = logical(individual.x);
end

if sum(x) == 0
    x(randi(numel(x))) = true;
end

% repair, items are considered in increasing order of profit/weight ratio,
% those achieving the lowest profit per weight unit are removed first
while any( sum( knapsack.weight(x,:) ) > knapsack.capacity )    
    x( knapsack.indexes( find(x(knapsack.indexes),1) ) ) = false;
end

if isfield(individual, 'topology')
    individual.s = x;
else
    individual.x = double(x);
end
individual.f = -sum( knapsack.profit(x,:), 1 )';


return