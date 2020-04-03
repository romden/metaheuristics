function [output] = SimulateANN(data, W1, W2, b1, b2, varargin)

% activation function
sigm = @(z) 1./(1+exp(-z));

% transformation function
mapping1 = @(input) sigm( input*W1 + repmat(b1, size(input,1), 1) ); % to hidden layer
mapping2 = @(input) sigm( input*W2 + repmat(b2, size(input,1), 1) ); % to logistic output unit
% mapping2 = @(input) input*W2 + repmat(b2, size(input,1), 1) ; % to linear output unit

% output data
output = mapping2( mapping1( data ) );

return
