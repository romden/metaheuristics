% function child = DEoperator(parent1, parent2, parent3, lb, ub)
% 
% % params
% n = numel(x); % chromosome length
% CR = 1.0; % crossover probability
% F = 0.5; % scaling factor
% pm = 1/n; % mutation probability
% etam = 20; % mutation distribution index
% 
% % init
% child = parent1;
% 
% n = numel(parent1);
% jr = randi(n);
% 
% for j = 1:n
%     
%     % DE operator
%     if rand <= CR || j == jr
%         child(j) = parent1(j) + F*(parent2(j)-parent3(j));
%     end
%     
%     % bounds
%     if child(j) < lb(j)
%         child(j) = lb(j) + rand*(parent1(j) - lb(j));
%     elseif child(j) > ub(j)
%         child(j) = ub(j) - rand*(ub(j) - parent1(j));
%     end
%     
%     % apply polynomial mutation, mutate offspring with probability pm
%     if rand <= pm
%         u = rand;
%         if u <= 0.5
%             delta = (2*u)^(1/(etam+1)) - 1;
%         else
%             delta = 1 - (2*(1 - u))^(1/(etam+1));
%         end
%         child(j) = child(j) + (ub(j)-lb(j))*delta;
%     end
%     
%     % bounds
%     if child(j) < lb(j)
%         child(j) = lb(j);
%     end
%     if child(j) > ub(j)
%         child(j) = ub(j);
%     end
% end
% 
% return


function [x] = DEoperator(x, x1, x2, x3, lb, ub)

% params
n = numel(x); % chromosome length
CR = 1.0; % crossover probability
F = 0.5; % scaling factor
pm = 1/n; % mutation probability
etam = 20; % mutation distribution index

% DE/rand/1/bin operator
jr = randi(n);
for j = 1:n
    if (rand <= CR || j == jr)
        x(j) = x1(j) + F*(x2(j) - x3(j));
    end
    
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
%     if x(j) < lb(j)
%         x(j) = lb(j);
%     elseif x(j) > ub(j)
%         x(j) = ub(j);
%     end
end

return