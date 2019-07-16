function [X,d,c,C,Z,Xlib,X_BB,bp] = vir_unmixing_denoiser_e_r(Alib,BB,Y,wv,varargin)
% solve the optimization problem
% 
%   minimize || Y-diag(d)(AX+B)-c1^T ||_{1,1} + lamda_a* || X ||_{1,1}
% 

c0 = [];
d0 = [];
lambda_lib = 0.01;
lambda_BB  = 0.5;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'C0'
                c0 = varargin{i+1};
            case 'D0'
                d0 = varargin{i+1};
            case 'LAMBDA_A_LIB'
                lambda_lib = varargin{i+1};
            case 'LAMBDA_A_BB'
                lambda_BB = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end


[L,Ny] = size(Y);
[L,Nlib] = size(Alib);
[L,Nbb] = size(BB);

A = [Alib BB];
% A = Alib;

%lambda_a = zeros(Nlib,1);
lambda_a = zeros(Nlib+Nbb,1);
lambda_a(1:Nlib) = lambda_lib;
lambda_a((Nlib+1):(Nlib+Nbb)) = lambda_BB;

% initialization: consistent distortions are initialized with zeros.
if isempty(d0)
    d = ones(L,1); % multiplicative component
else
    d = d0;
end
if isempty(c0)
    c = zeros(L,1); % additive component
else
    c = c0;
end

vec1L = ones(1,L);
vecNy1 = ones(Ny,1);
vec1Ny = ones(1,Ny);

% main loop
for i=1:1
    Ycd = (Y-c)./d;
    % minimize over X and B
    C = continuumDictionary(wv);
    s_c = vnorms(C,1);
    C = bsxfun(@rdivide,C,s_c);
    C = C*2;
    C(:,2:end-1) = -C(:,2:end-1);
    [X,Z,C,r,dd,rho,Rhov,res_p,res_d] = huwacbl1_gadmm_a_v2_vir(A,Ycd,wv,...
         'LAMBDA_A',lambda_a,'ConcaveBase',C,'tol',1e-5,'maxiter',1000,'verbose','no');
    
    Ym = A*X+C*Z;
    if i>0
        rr = Ycd-Ym;
        rr_median = nanmedian(rr,2);
        rr_median_out = rr-rr_median;
        bp_std = robust_v3('std',rr,2,'Noutliers',5);
%         bp_by_std = bp_std > 0.0006;
        bp_singles = rr_median_out > 10*bp_std;
        bp = bp_singles;
%         bp = or(bp_by_std,bp_singles);
    
        Ycd(bp_singles) = Ym(bp_singles);
        Y = d.*Ycd + c;
    end
    
    % minimize over d and c
    
    cvx_begin quiet
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
Xlib = X(1:Nlib,:);
X_BB = X(Nlib+1:end,:);
% [X,Z,C,r,d,rho,Rhov,res_p,res_d] = huwacbl1_gadmm_a_v2(A,Ycd,wv,...
%          'LAMBDA_A',lambda_a,'tol',1e-5,'maxiter',60,'verbose','yes');