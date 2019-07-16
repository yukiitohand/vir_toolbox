function [ yq ] = interpGaussConv_robust( x,y,xq,fwhm,varargin )
% [ yq ] = interpGaussian( x,y,xq,fwhm )
%   interpolation using Gaussian convolution
%   Input Parameters
%       x    : [L x 1], samples for original data
%       y    : [L x N], data for original samples, N: #of data
%       xq   : [Lq x 1], queried samples.
%       fwhm : [Lq x 1] or scalar, full width half maximum of each queried samples
%   Output Parameters
%       yq   : [Lq x N]

if ~isvector(x) || ~isnumeric(x)
    error('x must be a numeric vector');
end

if ~isvector(xq) || ~isnumeric(xq)
    error('xq must be a numeric vector');
end

if ~isvector(fwhm) || ~isnumeric(fwhm)
    error('xq must be a numeric vector');
end

L = length(x);
[Ly,N] = size(y);

if L~=Ly
    error('size of x and y does not match');
end

Lq = length(xq); Lqfwhm = length(fwhm);
if Lqfwhm==1
    fwhm = ones([Lq,1])*fwhm;
else
    if Lq~=Lqfwhm
         error('size of xq and fwhm does not match');
    end   
end
x = x(:); xq = xq(:); fwhm = fwhm(:);

retainratio = 0.9;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'RETAINRATIO'
                retainratio = varargin{i+1};
            otherwise
                error('Undefined option:%s.',upper(varargin{i}));
        end
    end
end

% covert fwhm to sigma
sigma = fwhm2sigma(fwhm);

yq = nan([Lq,N]);
% for i =1:Lq
%     c = normpdf(x,xq(i),sigma);
%     c_sum = sum(c);
%     yq(i,:) = sum(bsxfun(@times,c,y),1) /c_sum;
% end


% compute approximate width of the wvspc
x_extend = zeros([L+2,1]);
x_extend(2:end-1) = x;
x_extend(1) = 2*x(1)-x(2);
x_extend(end)=2*x(end)-x(end-1);

% detect blank region
dx = diff(x);
dx_med = medfilt1(dx,3);
dx_med(dx_med==0) = nan; % sometimes the median value get zero
bndry = find((dx ./ dx_med) > 10);
N_bndry = length(bndry);

% compute approximate width of the wvspc
blank_regions = zeros(N_bndry+2,2);
blank_regions(1,1) = 0;
blank_regions(N_bndry+2,2) = inf;
x_bd = zeros(L,1);
for n=1:(N_bndry+1)
    if n==1
        if isempty(bndry)
            iList = 1:L;
        else
            iList = 1:bndry(n);
        end
    elseif n <= N_bndry
        iList = (bndry(n-1)+1):bndry(n);
    else
        iList = bndry(n-1)+1:L;
    end
    Li = length(iList);
    xi = x(iList);
    x_extend = zeros([Li+2,1]);
    x_extend(2:end-1) = xi;
    x_extend(1) = 2*xi(1)-xi(2);
    x_extend(end)=2*xi(end)-xi(end-1);
    x_between = (x_extend(2:end) + x_extend(1:end-1))/2;
    x_bd(iList) = x_between(2:end) - x_between(1:end-1);
    
    
    blank_regions(n,2) = xi(1);
    blank_regions(n+1,1) = xi(end);
    
end


for i = 1:Lq
    % assess whether or not the center of position is in the blank region.
    xqi_is_invalid = any(and(xq(i)>blank_regions(:,1),xq(i)<blank_regions(:,2)));
    if ~xqi_is_invalid
        c = normpdf(x,xq(i),sigma(i)) .* x_bd;
        Z = sum(c);
        valid_idx = Z > retainratio;
        yq(i,valid_idx) = nansum(c.*y(:,valid_idx),1) ./ Z(1,valid_idx);
    end
end

end