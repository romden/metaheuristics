function [sigma] = ESMutationParameters(sigma)

n = numel(sigma);

SAdaptType = 'N'; % adaptation type
Tol = 1e-10; % inferior bound for step sizes

% initialize learning coefficents for self-adaptation rules
Eta1 = 1/sqrt(2*n);
Eta2 = 1/sqrt(2*sqrt(n));

lognormal = Eta1*randn(1);
for j = 1:n
    old_value = sigma(j);
    switch SAdaptType
        case 'I'
            % Isotropic rule
            sigma(j) = sigma(j)*exp(lognormal);
        case 'N'
            % Nonisotropic rule
            sigma(j) = sigma(j)*exp(lognormal + Eta2*randn(1));
    end
    % Test if inferior bound for step sizes is reached...
    if sigma(j) < Tol
        sigma(j) = old_value;
    end
end

return