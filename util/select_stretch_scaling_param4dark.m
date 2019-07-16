function [s_opt] = select_stretch_scaling_param4dark(y,dark,varargin)
% [s_opt] = select_stretch_scaling_param4dark(y,dark,varargin)
%  get the best scaling parameter by bisection method. Criterion is total
%  variation.
%  INPUTS
%    y: signal to be evaluated [L x N]
%    dark: dark signal [L x 1]
%  OUPUTS
%    s_opt: scalar, optimal scaling parameter
% 

[L,Ny] = size(y);

tol = 0.001;
s0  = 0.9;
s1  = 5;
denIdx = [];
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'TOL'
                tol = varargin{i+1};
            case 'S0'
                s0 = varargin{i+1};
            case 'S1'
                s1 = varargin{i+1};
            case 'DEN_IDX'
                denIdx = varargin{i+1};
            otherwise
                % Something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end


while abs(s1-s0)>tol
    yd0 = y-dark*s0;
    yd1 = y-dark*s1;
    if ~isempty(denIdx)
        if ischar(denIdx) && strcmpi(denIdx,'mean')
            yd0 = yd0 ./ nanmean(yd0,2);
            yd1 = yd1 ./ nanmean(yd1,2);
        elseif isscalar(denIdx)
            yd0 = yd0 ./ yd0(:,denIdx);
            yd1 = yd1 ./ yd1(:,denIdx);
        end
    end
    
    fs0 = total_variation(yd0,1);
    fs1 = total_variation(yd1,1);
    % fsm = total_variation(y-dark*sm);
    
    sm = (s0+s1)/2;
    if fs0>fs1
        s0 = sm;
    else
        s1 = sm;
    end
    
end

s_opt = sm;

end

