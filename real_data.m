

%% load the raw data processed into Matlab format
msd = readtable('MSD_standardized.csv', 'headerlines', 1);
msd = msd(1:5000, 2: 91); % 100000 * 90
msd = msd{:, :};
flt = readtable('flt_standardized.csv');
flt = flt(1:2000, 2:22); % 60448 * 21
flt = flt{:, :};

%% msd
[U, s, V] = svd(msd, 'econ');
s = diag(s);
n = size(msd, 1);
p = size(msd, 2);
r = floor(0.8 * n);
S = orth(randn(n, r))'; % S: r * n
msd_t = S * msd;
[V_t, s_t, U_t] = svd(msd_t,'econ');
s_t = diag(s_t);
ell_p = standard_spiked_inverse(s(1)^2, p/n);
ell_r = standard_spiked_inverse(s_t(1)^2, p/r);
T = ell_p / ell_r;
T

%% flt
[U, s, V] = svd(flt, 'econ');
s = diag(s);
n = size(flt, 1);
p = size(flt, 2);
r = floor(0.8 * n);
S = orth(randn(n, r))'; % S: r * n
flt_t = S * flt;
[V_t, s_t, U_t] = svd(flt_t, 'econ');
s_t = diag(s_t);
ell_p = standard_spiked_inverse(s(1)^2, p/n);
ell_r = standard_spiked_inverse(s_t(1)^2, p/r);
T = ell_p / ell_r;
T


%%
figure
hist(s.^2, 100)
hold on
hist(s_t.^2, 100, 'color', 'r')


