function [a,u,v] = fcn_edgets_v1(ts,zflag)
% This is an adapted version of the orignal scirpt written by Faskowitz, J. (https://github.com/brain-networks/edge-centric_demo/blob/main/fcn/fcn_edgets.m).
% Ziwei Zhang changed the z-score function to drop NaN before calculating the mean and
% standard deviation.

% fcn_edgets    
%
%   [a,u,v] = fcn_edgets(ts,zflag) 
%
%   Edge time series are a time-resolved account of the co-fluctuation
%   between nodes. These time series are an intermediate step in the 
%   calculation of Pearson correlation. 
%
%   Inputs:
%       ts, 
%           time series, size: (time)x(node)
%       zflag, 
%           do z-score? (time series should be z-scored to obtain ets)
%
%   Outputs:
%       a, 
%           edge time series (ets), size: (time)x(edge)
%       u,v, 
%           variables that are helpful for indexing the resultant ets
% 

if nargin == 1
    zflag = true;
end
[~,n] = size(ts);               % number samples/nodes
if zflag
    %z = zscore(ts);                 % z-score
    zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan')); % z-score dropping NaN
    z = zscor_xnan(ts);
else
    z = ts;
end
[u,v] = find(triu(ones(n),1));  % get edges
a = z(:,u).*z(:,v);             % edge ts products
