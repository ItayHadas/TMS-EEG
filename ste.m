function ster=ste(data,dim)%% return standard error of given data% data: a row or a colum vector of data% [n,m]=size(~isnan(data));% if n == 1        % case of a row vector%    data = data';%    n = m;% endif nargin==0    error('Not enough input arguments.');endswitch nargin    case 1        ster=nanstd(data)./sqrt(sum(~isnan(data)));    case 2        ster=nanstd(data,0,dim)./sqrt(sum(~isnan(data),dim));end