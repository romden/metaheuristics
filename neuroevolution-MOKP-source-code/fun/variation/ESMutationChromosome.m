function [x] = ESMutationChromosome(x, sigma, lb, ub)

n = numel(x);

for j = 1:n
    % mutation
    x(j) = x(j) + sigma(j)*randn(1);
    
    % bounds
%     if x(j) < lb(j),
%         x(j) = lb(j);
%     elseif x(j) > ub(j)
%         x(j) = ub(j);
%     end
end

return