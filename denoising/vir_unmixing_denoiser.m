function [X,d,c] = vir_unmixing_denoiser(A,Y,wv)
% [X,d,c,C,Z,Xlib,X_BB] = vir_unmixing_denoiser(Alib,BB,Y,wv,varargin)
%   de-noising vir I/F data with unmixing approach.
% 
%   minimize || Y-diag(d)(AX+B)-c1^T ||_{1,1} + lamda_a* || X ||_{1,1}
%
%  using an alternating minimization approach.
%  Denoising is not performed.
% 
%  INPUTS
%    A: [L,N] library matrix, whose columns storing endmembers
%    Y   : [L,Ny  ] observation spectra (I/F)
%    wv  : [L,1   ] wavelength vector
%  OUPUTS
%    X   : [Nlib+Nbb,Ny] abundance matrix
%    d   : [L,1   ] multiplicative correction factor
%    c   : [L,1   ] additive correction term
%   

[L,Ny] = size(Y);
[L,N] = size(A);

lambda_a = 0.01;

% initialization: consistent distortions are initialized with zeros.
d = ones(L,1); % multiplicative component
c = zeros(L,1); % additive component

vec1L = ones(1,L);
vecNy1 = ones(Ny,1);
vec1Ny = ones(1,Ny);


% main loop
for i=1:5
    Ycd = (Y-c)./d;
    % minimize over X and B
    C = continuumDictionary(wv);
    s_c = vnorms(C,1);
    C = bsxfun(@rdivide,C,s_c);
    C = C*2;
    C(:,2:end-1) = -C(:,2:end-1);
    [X,Z,C,r,d,rho,Rhov,res_p,res_d] = huwacbl1_gadmm_a_v2(A,Ycd,wv,...
         'LAMBDA_A',1,'tol',1e-5,'maxiter',1000,'verbose','yes','ConcaveBase',C);
    
    % minimize over d and c
    Ym = A*X+C*Z;
    cvx_begin
        variable D(L,L) diagonal nonnegative
        variable c(L)
        variable U(L,Ny) nonnegative
        minimize(  vec1L*U*vecNy1 )
        subject to
            -U <= Y-D*Ym-c*vec1Ny
            U  >= Y-D*Ym-c*vec1Ny
    cvx_end
    
    d = diag(D);
    
    
end
    
    




end