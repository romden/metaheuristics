function [child] = SBXcrossover(parent1, parent2, lb, ub, varargin)

% params
n = numel(parent1); % chromosome length
pc = 1; % crossover probability
etac = 20; % crossover distribution
pm = 1/n; % mutation probability
etam = 20; % mutation distribution index

% apply SBX (Simulated Binary Crossover) with probability pc
flag = false;
if rand <= pc
    flag = true;
end

child = parent1;

for j = 1:n
    
    % perform crossover
    if flag        
        % calculate beta
        u = rand;
        if u <= 0.5
            beta = power(2.0*u, 1.0/(etac+1.0));
        else
            beta = power(1.0/(2.0 - 2.0*u), 1.0/(etac+1.0));
        end
        
        % calculate Children gens, swaping them with probability 0.5
        if rand <= 0.5
            child(j) = 0.5*((1 + beta)*parent1(j) + (1 - beta)*parent2(j));
        else
            child(j) = 0.5*((1 - beta)*parent1(j) + (1 + beta)*parent2(j));
        end        
    end
    
    % polynomial mutation applied with probability pm
    if rand <= pm
        u = rand;
        if u <= 0.5
            delta = power(2.0*u, 1.0/(etam+1.0)) - 1.0;
        else
            delta = 1.0 - power(2.0 - 2.0*u, 1.0/(etam+1.0));
        end
        child(j) = child(j) + delta*(ub(j) - lb(j));
    end
    
    % bounds
%     if child(j) < lb(j),
%         child(j) = lb(j);
%     elseif child(j) > ub(j)
%         child(j) = ub(j);
%     end    
end

return