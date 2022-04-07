%% Figure 9-10: large signal, single spike
signal_strengths = linspace(1, 20, 30);
nrep = 20;
n = 4000;
p = 800;
r = 400;
gamma = p / n;
xi = r / n;
k = 1; 
num_rep = 20; 
print_every = 10;
lambda_res = zeros(length(signal_strengths), nrep);
cos_res = zeros(length(signal_strengths), nrep);
theory_cos = zeros(length(signal_strengths), 1);
theory_lambda = zeros(length(signal_strengths), 1);

% Figure 10: use the following Toeplitz Sigma
% q = 0.9;
% Sigma_t = toeplitz(q.^linspace(0, p - 1, p));
% rho_1 = 1;
% % rho_2 = (1 + q^2) / (1 - q^2);
% rho_2 = 1 + 2 / p * (p * q^2 / (1 - q^2) - q^2 * (1 - q^(2 * p)) / (1 - q^2)^2);

% Figure 9: use the following Sigma
% generate another Sigma
O = orth(rand(p, p));
L = ones(p, 1);
L(1) = 5;
L(2:floor(p / 2)) = 2;
Sigma = O' * diag(L) * O;
% 
rho_1 = mean(L);
rho_2 = mean(L .^ 2);

for i = 1:length(signal_strengths)
    disp(i);
    d = signal_strengths(i);
    D = diag(d);
    X = randn(n, p) * sqrtm(Sigma_t) / sqrt(n);
    W = orth(randn(n, k));
    U = orth(randn(p, k));
    signal_mat = W * D * U';
    Y = signal_mat + X;
    E = U' * Sigma_t * U;
    for j = 1:num_rep
        S = orth(randn(n, r))'; % S: r * n
        Y_t = S * Y;
        [~, D_t, U_t] = svd(Y_t, 'econ');
        lambda_res(i, j) = diag(D_t(1:k, 1:k)) .^ 2;
        cos_res(i, j) = (U_t(:, 1:k)' * U(:, 1:k)) .^ 2;
    end
    % theory
    theory_lambda(i) = max(xi * d ^ 2 + xi * E + gamma * rho_1, ...
        (1 + sqrt(gamma / xi)) * (sqrt(gamma) + sqrt(xi)) * sqrt(xi));

    temp = (xi - gamma / d ^ 4 * rho_2) / (xi + gamma / d ^ 2 ...
            * (rho_1 + (rho_2 - rho_1 * E) / d ^ 2));
    theory_cos(i) = max(temp, 0);
end

csvwrite('results/cos_large_simu.csv', cos_res);
csvwrite('results/cos_large_theory.csv', theory_cos);
csvwrite('results/lambda_large_simu.csv', lambda_res);
csvwrite('results/lambda_large_theory.csv', theory_lambda);


fig_path = 'plots/';
name = sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi);
plot_result(cos_res', theory_cos, signal_strengths, n, k, gamma, xi, num_rep, strcat(fig_path, 'cos_large_orth_'), 'cos', name)
plot_result(lambda_res', theory_lambda, signal_strengths, n, k, gamma, xi, num_rep, strcat(fig_path, 'lambda_large_orth_'), 'lbd', name)

