function [W, b, neurons] = DecodeIndividual(individual, varargin)

numLayers = numel(individual.topology);

W = cell(numLayers-1, 1);
b = cell(numLayers-1, 1);
neurons  = (individual.x(end-individual.topology(2)+1:end) >= 0);

head = 0;
tail = 0;
for l = 2:numLayers
    % decode weights
    head = tail+1;
    tail = tail + individual.topology(l-1)*individual.topology(l);
    
    W{l-1} = reshape(individual.x(head:tail), individual.topology(l-1), individual.topology(l));
    
    % decode biases
    head = tail+1;
    tail = tail + individual.topology(l);
    b{l-1} = reshape(individual.x(head:tail), 1, individual.topology(l));
end

return
