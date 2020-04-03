function [x] = PolynomialMutation(x, lb, ub, varargin)

% params
n = numel(x); % chromosome length
pm = 1/n; % mutation probability
etam = 20; % mutation distribution index

% apply polynomial mutation
for j = 1:n    
    % polynomial mutation applied with probability pm
    if rand <= pm
        u = rand;
        if u <= 0.5
            delta = power(2.0*u, 1.0/(etam+1.0)) - 1.0;
        else
            delta = 1.0 - power(2.0 - 2.0*u, 1.0/(etam+1.0));
        end
        x(j) = x(j) + delta*(ub(j) - lb(j));
    end
    
    % bounds
%     if x(j) < lb(j),
%         x(j) = lb(j);
%     elseif x(j) > ub(j)
%         x(j) = ub(j);
%     end    
end

return