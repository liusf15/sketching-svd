function [] = large_signal_toe(n, p, r, k, num_rep, seed, print_every)
    gamma = p / n;
    xi = r / n;
    d = linspace(30, 5, k);
    D = diag(d);
    rng(seed);
    cos_large_toe = zeros(num_rep, k, k);
    lambda_large_toe = zeros(num_rep, k);
    theory_cos_large_toe = zeros(k, k);
    theory_lambda_large_toe = zeros(k, 1);
    rng(1);
    % generate Sigma
    q = 0.9;
    Sigma_t = toeplitz(q.^linspace(0, p - 1, p));
%     Sigma_t = Sigma_t * 2; % for q = 0.5
    rho_1 = 1;
    rho_2 = (1 + q^2) / (1 - q^2);
    
    % simulation
    X = randn(n, p) * sqrtm(Sigma_t) / sqrt(n);
    W = orth(randn(n, k));
    U = orth(randn(p, k));
    signal_mat = W * D * U';
    Y = signal_mat + X;
    E = U' * Sigma_t * U;

    for j = 1:num_rep
        if mod(j, print_every) == 0
            fprintf('large signal %d \n', j);
        end
        S = orth(randn(n, r))'; % S: r * n
        Y_t = S * Y;
        [~, D_t, U_t] = svd(Y_t, 'econ');
        lambda_large_toe(j, :) = diag(D_t(1:k, 1:k)) .^ 2;
        cos_large_toe(j, :, :) = (U_t(:, 1:k)' * U(:, 1:k)) .^ 2;
    end


    % theory

    E_diag = diag(E);
    for i = 1:k
        theory_lambda_large_toe(i) = max(xi * d(i) ^ 2 + xi * E_diag(i) + gamma * rho_1, ...
            (1 + sqrt(gamma / xi)) * (sqrt(gamma) + sqrt(xi)) * sqrt(xi));
        temp = (xi - gamma / d(i) ^ 4 * rho_2) / (xi + gamma / d(i) ^ 2 ...
                * (rho_1 + (rho_2 - rho_1 * E_diag(i)) / d(i) ^ 2));
        theory_cos_large_toe(i, i) = max(temp, 0);
        for j = 1:k
            if j ~= i
                theory_cos_large_toe(i, j) = max(temp * E(i, j) ^ 2 / (d(i) ^ 2 - d(j) ^ 2) ^ 2, 0);
            end
        end
    end

    % plot cos
    me = mean(cos_large_toe, 1);
    me = reshape(me, k, k);
    sd = std(cos_large_toe, 1);
    sd = reshape(sd, k, k);

    figure;
    hold on
    errorbar(d, diag(me), diag(sd), ...
        'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9]);
    plot(d, diag(theory_cos_large_toe), 'DisplayName', 'Theory', ...
         'linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
    xlabel('$$d_i$$', 'Interpreter', 'LaTex');
    ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
    set(gca,'fontsize',20)
    title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
    legend('location','southeast');
    grid on;

    filename = sprintf('cos_toep_q_%0.1f_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', q, gamma, xi, num_rep, n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)

    % plot lambda
    figure;
    hold on;
    errorbar(d, mean(lambda_large_toe, 1), std(lambda_large_toe, 1), ...
        'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9]);
    plot(d, theory_lambda_large_toe, 'DisplayName', 'Theory', ...
        'linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
    xlabel('$$d_i$$', 'Interpreter', 'LaTex')
    ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
    set(gca,'fontsize',20)
    title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
    legend('location','southeast');
    grid on;

    
    filename = sprintf('lambda_toep_q_%0.1f_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', q, gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)








