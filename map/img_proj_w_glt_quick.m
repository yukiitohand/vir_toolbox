function [img_proj] = img_proj_w_glt_quick(img,glt_x,glt_y)
% [img_proj] = img_proj_w_glt_quick(img,GLTdata)
%   projected image using GLT
%  Inputs
%   img: image data [L x S x B]
%   glt_x: glt_x image [Lg x Sg x Bg]
%   glt_y: glt_y image [Lg x Sg x Bg]
%  Outputs
%   img_proj: projected image, number of bands is same, and image size is
%   the size of the GLT image


[L,S,B] = size(img);

[Lg,Sg,Bg] = size(glt_x);

X_GLT = abs(glt_x);
Y_GLT = abs(glt_y);

cumIdx = Y_GLT(:)+(X_GLT(:)-1)*L;

img_proj_2d = nan(B,length(cumIdx));

cumIdx_nonzero_idx = (cumIdx>0);
cumIdx_nonzero = cumIdx(cumIdx_nonzero_idx);

flat_img_2d = reshape(img,[S*L,B])';
img_proj_2d_nozero = flat_img_2d(:,cumIdx_nonzero);

img_proj_2d(:,cumIdx_nonzero_idx) = img_proj_2d_nozero;

img_proj = reshape(img_proj_2d',[Lg,Sg,B]);

end