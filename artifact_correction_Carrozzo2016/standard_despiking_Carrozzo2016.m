function [im_in_outrmvd,outliers] = standard_despiking_Carrozzo2016(im_in,wv)
% [im_in_outrmvd,outliers] = standard_despiking_Carrozzo2016(im_in)
%  standard despiking procedures, detecting and replacing spikes
% INPUTS
%   im_in: [B x S] each column is a signal
%   wv   : [B x 1] wavelength
% OUTPUTS
%  im_in_out_rmvd: [B x S] spike replaced signals
%  outliers      : boolean matrix [B x S], true if outliers

[B,S] = size(im_in);

h = ones(3,1)/3;
im_in_fil = imfilter(im_in,h);

r = im_in./im_in_fil;
r_sgm = nanmean(((r-1).^2),'all');
r_std = sqrt(r_sgm);

outliers = abs(r-1) > 3*r_std;

% interpolate outliers by a polynomial fit with 20 neighboring points.
im_in_outrmvd = im_in;
for c=1:S
    for b=1:B
        if outliers(b,c)
            bdx_tested = max(1,(b-10)):min(B,(b+10));
            x = wv(bdx_tested);
            y = im_in(bdx_tested,c);
            good_tmp = ~outliers(bdx_tested);
            p = polyfit(x(good_tmp),y(good_tmp),2);
            im_in_outrmvd(b,c) = polyval(p,wv(b));
            
        end
    end
end