% Simulation with single spike
% Figures: 3, 4, 5, 6, 7, 8, 13

n = 4000;
p = 800;
r = 200;
gamma = p / n;
xi = r / n;
k = 1;
rng(1);
signal_strengths = linspace(1, 20, 10);

num_rep = 20;
names = ["orthogonal projection", "Gaussian projection", "uniform sampling", "Hadamard", "countSketch", "countSketch-normalized", "leverage", "osnap"];
cos = zeros(length(names), 10, num_rep);
lambda = zeros(length(names), 10, num_rep);
%% simulation with all sketching methods considered in the paper
for j = 1:num_rep
    disp(j);
    for i = 1:10
        d = signal_strengths(i);
%         X = randn(n, p) / sqrt(n);  % Gaussian
        X = (rand(n, p) * 2 - 1) * sqrt(3) / sqrt(n); % uniform data
%         X = (binornd(1, 0.5, n, p) * 2 - 1) / sqrt(n);  % Rademacher data
        W = orth(randn(n, k));
        U = orth(randn(p, k));
        signal_mat = W * d * U';
        Y = signal_mat + X;
        [cos(1, i, j), lambda(1, i, j)] = sketchingMethods.orth(Y, U, r);  % orthogonal projection
        [cos(2, i, j), lambda(2, i, j)] = sketchingMethods.gaus(Y, U, r);  % Gaussian projection
        [cos(3, i, j), lambda(3, i, j)] = sketchingMethods.unif(Y, U, r);  % uniform sampling
        [cos(4, i, j), lambda(4, i, j)] = sketchingMethods.hada(Y, U, r);  % SRHT (Y has n equal to pow2(floor(log2(n))))
        [cos(5, i, j), lambda(5, i, j)] = sketchingMethods.coun(Y, U, r, false);  % countSketch (unnormalized)
        [cos(6, i, j), lambda(6, i, j)] = sketchingMethods.coun(Y, U, r, true);  % countSketch (normalized)
        [cos(7, i, j), lambda(7, i, j)] = sketchingMethods.leve(Y, U, r);  % leverage score sampling
        [cos(8, i, j), lambda(8, i, j)] = sketchingMethods.osna(Y, U, r, 40);  % osnap

    end
end
% save all simulation results for later use
if ~exist('results/', 'dir')
       mkdir('results/')
end
for i = 1:8
    filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
    csvwrite(strcat('results/uniX_cos_', names(i), filename), squeeze(cos(i, :, :)))
    csvwrite(strcat('results/uniX_lambda_', names(i), filename), squeeze(lambda(i, :, :)))
end

%% Make Figure 13
r = 200;
idx_sort = [7, 1, 3, 4, 6, 5, 2, 8];
for i = 1:8
    filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
    cos(i, :, :) = csvread(strcat('results/uniX_cos_', names(i), filename));
    lambda(i, :, :) = csvread(strcat('results/uniX_lambda_', names(i), filename));
end
% cos
figure, hold on;
mark = {':', '-', ':', '-', '--', '-.', '-.', '--', '-.'};
for j = 1:8
    i = idx_sort(j);
    errorbar(signal_strengths, mean(cos(i, :, :), 3), std(cos(i, :, :), 0, 3),'lineWidth', 2, 'DisplayName', names(i), 'linestyle', mark{j});
end
legend('location','northeast');
legend('off')
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
ylim([0, 1]);
set(gca,'fontsize',20)
legend('location','southeast');
grid on;
ax1 = gca;
axes(ax1);
magnifyOnFigure;

filename = sprintf('plots/final/magnify_uniX_cos_all_n_%d_p_%d_r_%d_k_%d_nrep_%d.pdf', n, p, r, k, num_rep);
saveTightFigure(gcf, filename);
close(gcf);

% lambda
figure, hold on;
for j = 1:8
    i = idx_sort(j);
    errorbar(signal_strengths, mean(lambda(i, :, :), 3), std(lambda(i, :, :), 0, 3),'lineWidth', 3, 'DisplayName', names(i), 'linestyle', ':');
end
plot([0 20], [0, 20].^2, 'lineWidth', 3, 'linestyle', "-", 'DisplayName', "y=x")
legend('location','northwest');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize', 14)
grid on;
filename = sprintf('plots/compare_all/uniX_lambda_all_n_%d_p_%d_r_%d_k_%d_nrep_%d.png', n, p, r, k, num_rep);
saveas(gcf, filename);
close(gcf);

%% Make Figures 3, 4, 5, 6, 7, 8
fig_path = 'plots/single_uniX_';

% Figure 3: orthogonal projection
i = 1;
filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
cos = csvread(strcat('results/uniX_cos_', names(i), filename));
lambda = csvread(strcat('results/uniX_lambda_', names(i), filename));
theory_cos = zeros(10, 1);
theory_lambda = zeros(10, 1);
for j = 1:10
    [theory_cos(j), theory_lambda(j)] = sketchingMethods.theory_orth(gamma, xi, signal_strengths(j));
end
plot_result(cos', theory_cos, signal_strengths, n, k, gamma, xi, num_rep, strcat(fig_path, 'cos_orth_'), 'cos')
plot_result(lambda', theory_lambda, signal_strengths, n, k, gamma, xi, num_rep, strcat(fig_path, 'lambda_orth_'), 'lbd')


% Figure 4: iid projection
i = 2;
filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
cos = csvread(strcat('results/uniX_cos_', names(i), filename));
lambda = csvread(strcat('results/uniX_lambda_', names(i), filename));
theory_cos = zeros(10, 1);
theory_lambda = zeros(10, 1);
theory_cos_approx = zeros(10, 1);
theory_lambda_approx = zeros(10, 1);
for j = 1:10
    [theory_cos(j), theory_lambda(j), theory_cos_approx(j), theory_lambda_approx(j)] = sketchingMethods.theory_gaus(gamma, xi, signal_strengths(j));
end
%%% plot cos
figure, hold on
a = {'-','--','-.',':'};
me = mean(cos', 1);
re = std(cos', 1);
h1 = errorbar(signal_strengths, me, re, 'lineWidth', 3, 'DisplayName', ...
    'Simulation', 'color', [.1, .4, .9]);
h2 = plot(signal_strengths, theory_cos, 'DisplayName', 'Theory', ...
    'linewidth', 4, 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(signal_strengths, theory_cos_approx, 'LineStyle','--', ...
    'linewidth', 4, 'color', [1, .5, 0], 'DisplayName', 'Theory approx');
xlabel('$$d_1$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_1,\tilde\xi_1\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
legend('location','southeast');
grid on;
filename = sprintf('plots/single_uniX_cos_gaus_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', gamma, xi,num_rep,n, k);
saveTightFigure(gcf, filename);
fprintf(['Saved Results to ' filename '\n']);
close(gcf)
%%% plot lambda
figure, hold on
me = mean(lambda', 1);
re = std(lambda', 1);
h1 = errorbar(signal_strengths, me, re, 'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = plot(signal_strengths, theory_lambda, 'DisplayName', 'Theory','linewidth', 4, 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(signal_strengths, theory_lambda_approx, 'DisplayName', 'Theory approx','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_1$$', 'Interpreter', 'LaTex')
ylabel('$$\tilde\lambda_1$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
legend('location','southeast');
grid on;
filename = sprintf('plots/single_uniX_lambda_gaus_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', gamma, xi,num_rep,n, k);
saveTightFigure(gcf, filename);
fprintf(['Saved Results to ' filename '\n']);
close(gcf)

% Figure 5: uniform sampling
i = 3;
filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
cos = csvread(strcat('results/uniX_cos_', names(i), filename));
lambda = csvread(strcat('results/uniX_lambda_', names(i), filename));
theory_cos = zeros(10, 1);
theory_lambda = zeros(10, 1);
for j = 1:10
    [theory_cos(j), theory_lambda(j)] = sketchingMethods.theory_unif(gamma, xi, signal_strengths(j));
end
plot_result(cos', theory_cos, signal_strengths, n, k, gamma, xi, num_rep, strcat(fig_path, 'cos_unif_'), 'cos')
plot_result(lambda', theory_lambda, signal_strengths, n, k, gamma, xi, num_rep, strcat(fig_path, 'lambda_unif_'), 'lbd')

% Figure 6: SRHT
i = 4;
filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
cos = csvread(strcat('results/uniX_cos_', names(i), filename));
lambda = csvread(strcat('results/uniX_lambda_', names(i), filename));
theory_cos = zeros(10, 1);
theory_lambda = zeros(10, 1);
for j = 1:10
    [theory_cos(j), theory_lambda(j)] = sketchingMethods.theory_hada(gamma, xi, signal_strengths(j));
end
plot_result(cos', theory_cos, signal_strengths, n, k, gamma, xi, num_rep, strcat(fig_path, 'cos_hada_'), 'cos')
plot_result(lambda', theory_lambda, signal_strengths, n, k, gamma, xi, num_rep, strcat(fig_path, 'lambda_hada_'), 'lbd')

% Figure 7, 8: countSketch
n = 500; %20,50,100,500
p = 500;
xi = 0.2;
r = floor(n * xi);
gamma = p / n;

%% simulation
cos_all = zeros(2, 10, num_rep);
lambda_all = zeros(2, 10, num_rep);
for j = 1:num_rep
    disp(j);
    for i = 1:10
        d = signal_strengths(i);
        X = (rand(n, p) * 2 - 1) * sqrt(3) / sqrt(n); 
        W = orth(randn(n, k));
        U = orth(randn(p, k));
        signal_mat = W * d * U';
        Y = signal_mat + X;
        [cos_all(1, i, j), lambda_all(1, i, j)] = sketchingMethods.coun(Y, U, r, false);
        [cos_all(2, i, j), lambda_all(2, i, j)] = sketchingMethods.coun(Y, U, r, true);
        
    end
end

% filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
% cos = csvread(strcat('results/uniX_cos_', names(5), filename))';
% lambda = csvread(strcat('results/uniX_lambda_', names(5), filename))';
% cos_norm = csvread(strcat('results/uniX_cos_', names(6), filename))';
% lambda_norm = csvread(strcat('results/uniX_lambda_', names(6), filename))';
theory_cos = zeros(10, 1);
theory_lambda = zeros(10, 1);
for j = 1:10
    [theory_cos(j), theory_lambda(j)] = sketchingMethods.theory_coun(gamma, xi, signal_strengths(j));
end
%%% plot cos
cos = squeeze(cos_all(1, :, :));
cos_norm = squeeze(cos_all(2, :, :));
figure, hold on
h1 = errorbar(signal_strengths, mean(cos, 2), std(cos'),'lineWidth', 3, 'DisplayName', 'CS', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = errorbar(signal_strengths, mean(cos_norm, 2), std(cos_norm'),'lineWidth', 3, 'DisplayName', 'CS Normalized', 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(signal_strengths, theory_cos, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_1$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_1,\tilde\xi_1\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',30)
legend('location','southeast');
grid on;
filename = sprintf('plots/single_uniX_cos_count_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', gamma, xi,num_rep,n, k);
saveTightFigure(gcf, filename);
fprintf(['Saved Results to ' filename '\n']);
close(gcf)

%%% plot lambda
lambda = squeeze(lambda_all(1, :, :));
lambda_norm = squeeze(lambda_all(2, :, :));
figure, hold on
h1 = errorbar(signal_strengths, mean(lambda * xi, 2), std(lambda' * xi),'lineWidth', 3, 'DisplayName', 'CS', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = errorbar(signal_strengths, mean(lambda_norm, 2), std(lambda_norm'),'lineWidth', 3, 'DisplayName', 'CS Normalized', 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(signal_strengths, theory_lambda, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_1$$', 'Interpreter', 'LaTex')
ylabel('$$\tilde\lambda_1$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize', 30)
legend('location','northwest');
grid on;
filename = sprintf('plots/single_uniX_lambda_count_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.pdf', gamma, xi,num_rep,n, k);
saveTightFigure(gcf, filename);
fprintf(['Saved Results to ' filename '\n']);
close(gcf)


