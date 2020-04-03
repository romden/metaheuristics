function [child] = BitFlip(child)

% params
n = numel(child); % chromosome length
pm = 1/n; % mutation probability

for j = 1:n
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
