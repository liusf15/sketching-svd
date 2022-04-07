% Simulation multi spikes
% Figures: 11, 12

n = 20000;
p = 2500;
r = 2000;
gamma = p / n;
xi = r / n;
k = 5;
rng(1);
d = linspace(20, 3, k);
D = diag(d);

num_rep = 20;
names = ["orthogonal projection", "Gaussian projection", "uniform sampling", "Hadamard", "countSketch", "countSketch-normalized", "leverage", "osnap"];
cos = zeros(length(names), k, num_rep);
lambda = zeros(length(names), k, num_rep);
%% simulation
for j = 1:num_rep
    disp(j);
%     X = randn(n, p) / sqrt(n);  % Gaussian
%     X = (rand(n, p) * 2 - 1) * sqrt(3) / sqrt(n); % uniform data
    X = (binornd(1, 0.5, n, p) * 2 - 1) / sqrt(n);  % Rademacher data
    W = orth(randn(n, k));
    U = orth(randn(p, k));
    signal_mat = W * D * U';
    Y = signal_mat + X;
    [cos(1, :, j), lambda(1, :, j)] = sketchingMethods.orth(Y, U, r);
    [cos(2, :, j), lambda(2, :, j)] = sketchingMethods.gaus(Y, U, r);
    [cos(3, :, j), lambda(3, :, j)] = sketchingMethods.unif(Y, U, r);
    [cos(4, :, j), lambda(4, :, j)] = sketchingMethods.hada(Y, U, r);  % Y has n equal to pow2(floor(log2(n)))
    [cos(5, :, j), lambda(5, :, j)] = sketchingMethods.coun(Y, U, r, false);
    [cos(6, :, j), lambda(6, :, j)] = sketchingMethods.coun(Y, U, r, true);
    [cos(7, :, j), lambda(7, :, j)] = sketchingMethods.leve(Y, U, r);
    [cos(8, :, j), lambda(8, :, j)] = sketchingMethods.osna(Y, U, r, 40);
end
if ~exist('results/', 'dir')
       mkdir('results/')
end
if ~exist('plots/', 'dir')
       mkdir('plots/')
end
for i = 1:8
    filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
    csvwrite(strcat('results/cos_rade_', names(i), filename), squeeze(cos(i, :, :)))
    csvwrite(strcat('results/lambda_rade_', names(i), filename), squeeze(lambda(i, :, :)))
end

%% plot
% 1. orthogonal
i = 1;
filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
cos = csvread(strcat('results/cos_', names(i), filename));
lambda = csvread(strcat('results/lambda_', names(i), filename));
[theory_cos, theory_lambda] = sketchingMethods.theory_orth(gamma, xi, d);
plot_result(cos', theory_cos, d, n, k, gamma, xi, num_rep, 'plots/multi_cos_orth_', 'cos', 'Orthogonal projection', 35)
plot_result(lambda', theory_lambda, d, n, k, gamma, xi, num_rep, 'plots/multi_lambda_orth_', 'lbd', 'Orthogonal projection', 35)


% 2. gaussian
i = 2;
filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
cos = csvread(strcat('results/cos_', names(i), filename));
lambda = csvread(strcat('results/lambda_', names(i), filename));
[theory_cos, theory_lambda, theory_cos_approx, theory_lambda_approx] = sketchingMethods.theory_gaus(gamma, xi, d);

%%% plot cos
figure, hold on
a = {'-','--','-.',':'};
me = mean(cos', 1);
re = std(cos', 1);
h1 = errorbar(d, me, re, 'lineWidth', 3, 'DisplayName', ...
    'Simulation', 'color', [.1, .4, .9]);
h2 = plot(d, theory_cos, 'DisplayName', 'Theory', ...
    'linewidth', 4, 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(d, theory_cos_approx, 'LineStyle','--', ...
    'linewidth', 4, 'color', [1, .5, 0], 'DisplayName', 'Theory approx');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize', 35)
title('Gaussian projection', 'fontsize', 35, 'Interpreter', 'LaTex');
legend('location','southeast');
grid on;
filename = sprintf('plots/multi_cos_gaus_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', gamma, xi,num_rep,n, k);
saveTightFigure(gcf, filename);
fprintf(['Saved Results to ' filename '\n']);
close(gcf)
%%% plot lambda
figure, hold on
me = mean(lambda', 1);
re = std(lambda', 1);
h1 = errorbar(d, me, re, 'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = plot(d, theory_lambda, 'DisplayName', 'Theory','linewidth', 4, 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(d, theory_lambda_approx, 'DisplayName', 'Theory approx','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex')
ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize', 35)
title('Gaussian projection', 'FontSize', 35, 'Interpreter', 'LaTex');
legend('location','northwest');
grid on;
filename = sprintf('plots/multi_lambda_gaus_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', gamma, xi,num_rep,n, k);
saveTightFigure(gcf, filename);
fprintf(['Saved Results to ' filename '\n']);
close(gcf)

% 3. uniform
i = 3;
filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
cos = csvread(strcat('results/cos_', names(i), filename));
lambda = csvread(strcat('results/lambda_', names(i), filename));
[theory_cos, theory_lambda] = sketchingMethods.theory_unif(gamma, xi, d);
plot_result(cos', theory_cos, d, n, k, gamma, xi, num_rep, 'plots/multi_cos_unif_', 'cos', 'Uniform sampling', 35)
plot_result(lambda', theory_lambda, d, n, k, gamma, xi, num_rep, 'plots/multi_lambda_unif_', 'lbd', 'Uniform sampling', 35)

% 4. Hadamard
i = 4;
filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
cos = csvread(strcat('results/cos_', names(i), filename));
lambda = csvread(strcat('results/lambda_', names(i), filename));
[theory_cos, theory_lambda] = sketchingMethods.theory_hada(gamma, xi, d);
plot_result(cos', theory_cos, d, n, k, gamma, xi, num_rep, 'plots/multi_cos_hada_', 'cos', 'SRHT', 35)
plot_result(lambda', theory_lambda, d, n, k, gamma, xi, num_rep, 'plots/multi_lambda_hada_', 'lbd', 'SRHT', 35)

% 5. Count sketch
filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
cos = csvread(strcat('results/cos_', names(5), filename))';
lambda = csvread(strcat('results/lambda_', names(5), filename))';
cos_norm = csvread(strcat('results/cos_', names(6), filename))';
lambda_norm = csvread(strcat('results/lambda_', names(6), filename))';
[theory_cos, theory_lambda] = sketchingMethods.theory_coun(gamma, xi, d);
%%% plot cos
figure, hold on
h1 = errorbar(d, mean(cos, 1), std(cos, 1),'lineWidth', 3, 'DisplayName', 'CS', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = errorbar(d, mean(cos_norm, 1), std(cos_norm, 1),'lineWidth', 3, 'DisplayName', 'CS Normalized', 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(d, theory_cos, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize', 35)
title('Count sketch', 'FontSize', 35, 'Interpreter', 'LaTeX')
legend('location','southeast');
grid on;
filename = sprintf('plots/multi_cos_count_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', gamma, xi,num_rep,n, k);
saveTightFigure(gcf, filename);
fprintf(['Saved Results to ' filename '\n']);
close(gcf)

%%% plot lambda
figure, hold on
h1 = errorbar(d, mean(lambda * xi, 1), std(lambda * xi, 1),'lineWidth', 3, 'DisplayName', 'CS', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = errorbar(d, mean(lambda_norm, 1), std(lambda_norm, 1),'lineWidth', 3, 'DisplayName', 'CS Normalized', 'color', [0, .8, 0], 'linestyle', '--');
h3 = plot(d, theory_lambda, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', ':');
xlabel('$$d_i$$', 'Interpreter', 'LaTex')
ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize', 35)
title('Count sketch', 'FontSize', 35, 'Interpreter', 'LaTeX')
legend('location','northwest');
grid on;
filename = sprintf('plots/multi_lambda_count_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', gamma, xi,num_rep,n, k);
saveTightFigure(gcf, filename);
fprintf(['Saved Results to ' filename '\n']);
close(gcf)




