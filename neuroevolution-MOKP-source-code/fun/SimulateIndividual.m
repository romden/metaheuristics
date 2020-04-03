function [output] = SimulateIndividual(individual, input, varargin)

[cW, cb, neurons] = DecodeIndividual(individual);

% output data
output = SimulateANN(input, cW{1}(:,neurons), cW{2}(neurons,:), cb{1}(neurons), cb{2} );

e = 1e-4;
output(output < e) = 0;
output(output > 1-e) = 1;

return
