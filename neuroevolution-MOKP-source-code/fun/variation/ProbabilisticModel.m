function child = ProbabilisticModel(parent, data, varargin)

[k, n] = size(data);

% build a model
data = data - repmat(mean(data), k, 1);
C = data'*data/(k-1);
[B, D] = eig(C);

% generate offspring
child = parent + B*(sqrt(diag(abs(D))).*randn(n, 1));

return