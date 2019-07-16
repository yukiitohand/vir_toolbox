function [tval] = total_variation(x,dim)
% [tval] = total_variation(x,dim)
%  evaluate total variation


if dim>1
    nd = ndims(x);
    dim_complement = setdiff(1:nd,dim);
    x = permute(x,[dim dim_complement]);
end


adff = abs(diff(x));

tval = sum(adff,'all');

end

