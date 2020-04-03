function [x] = Repair(x, varargin)

n = numel(x);
Cmin = 3;
Cmax = n;

nrOnes = sum(x >= 0);
if nrOnes < Cmin
    idx = find(x < 0);
    idx = idx(randperm(numel(idx), Cmin-nrOnes));
    x(idx) = abs(x(idx));
elseif nrOnes > Cmax
    idx = find(x >= 0);
    idx = idx(randperm(numel(idx), nrOnes-Cmax));
    x(idx) = -x(idx);
    x(idx(x(idx)==0)) = -1e-4;
end

return



% function [x] = Repair(x, varargin)
% 
% n = numel(x);
% Cmin = 3;
% Cmax = n;
% 
% nrOnes = sum(x);
% if nrOnes < Cmin
%     if islogical(x)
%         idx0 = find(~x);
%     else
%         idx0 = find(x == 0);
%     end
%     idx = idx0(randperm(numel(idx0), Cmin-nrOnes));
%     x(idx) = true;
% elseif nrOnes > Cmax
%     if islogical(x)
%         idx1 = find(x);
%     else
%         idx1 = find(x == 1);
%     end
%     idx = idx1(randperm(numel(idx1), nrOnes-Cmax));
%     x(idx) = false;
% end
% 
% return

% function [x] = Repair(x, varargin)
% 
% n = numel(x);
% Cmin = 3;
% Cmax = n;
% 
% nrOnes = sum(x);
% if nrOnes < Cmin
%     if islogical(x)
%         idx0 = find(~x);
%     else
%         idx0 = find(x == 0);
%     end
%     idx = idx0( randperm0(numel(idx0), Cmin-nrOnes) );
%     x(idx) = true;
% elseif nrOnes > Cmax
%     if islogical(x)
%         idx1 = find(x);
%     else
%         idx1 = find(x == 1);
%     end
%     idx = idx1( randperm0(numel(idx1), nrOnes-Cmax) );
%     x(idx) = false;
% end
% 
% return
% 
% function [x] = randperm0(l, n, varargin)
% 
% x = randperm(l);
% x = x(1:n);
% 
% return