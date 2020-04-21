%% HGDP 
cd('C:\Dropbox\Projects\Sketching SVD\Experiments\Exp4-Empirical')
addpath('C:\Dropbox\Projects\Sketching SVD\Experiments\Exp4-Empirical\hgdp\chr22')
%% load the raw data processed into Matlab format
clear; 
load('chr22_geno_norm.mat');
R = c22_norm; clear('c22_norm');
%set to 0 missing entries;
R(isnan(R))=0;

rng(2);
%% three options for matrix to look at
%report results for all three
%1. full matrix
%Y = R;

%2. subsample dimension
%Y = R(:,1:20:end);

%3. subsample dimension and sample size
Y = R(1:10:end,1:20:end);
%%
[n,p] = size(Y);
Y = n^(-1/2)*Y;
%% 
tic;
[U,s,V] = svd(Y,'econ');
toc;
s =diag(s);

%%
tic;
r = floor(0.8*n);
S = orth(randn(n, r))'; % S: r * n
Y_t = S*Y;
[V_t, s_t, U_t] = svd(Y_t,'econ');
toc;
s_t =diag(s_t);

%%
%T=\frac{\lambda^{-1}\left(\sigma_1(Y),p/n\right)}{\lambda^{-1}\left(\sigma_1(SY),r/n\right)}.

%[ell,cos_right,cos_left] = standard_spiked_inverse(lambda,gamma)

ell_p = standard_spiked_inverse(s(1)^2,p/n);
ell_r = standard_spiked_inverse(s_t(1)^2,p/r);
T = ell_p/ell_r;
%generally get values between 1.2-1.35

%%
figure
hist(s.^2, 100)
figure
hist(s_t.^2, 100)


