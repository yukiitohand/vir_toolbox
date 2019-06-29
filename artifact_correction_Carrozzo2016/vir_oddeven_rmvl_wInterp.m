function [xs] = vir_oddeven_rmvl_wInterp(x,wv)
% [xs] = vir_oddeven_rmvl_wInterp(x,wv)
%  
%  INPUTS
%   x:  [L x N]
%   wv: wavelength 
%  OUTPUTS
%   xs: [L x N]

% [L,N] = size(x);

wv_between = (wv(2:end)+wv(1:end-1)) / 2;

x_between = (x(2:end,:)+x(1:end-1,:)) / 2;

xs = interp1(wv_between,x_between,wv,'linear');

end


