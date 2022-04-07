% compare SRHT, uniform sampling, CountSketch, Gaussian projection
% Figure 2

n = 2500;
p = 800;
r = 2000;
gamma = p / n;
xi = r / n;
k = 7;
rng(1);
d = linspace(20, 3, k);
D = diag(d);

num_rep = 20;
names = ["Hadamard", "uniform sampling", "countSketch", "countSketch-normalized", "Gaussian projection"];
cos = zeros(length(names), k, num_rep);
lambda = zeros(length(names), k, num_rep);
%% simulation
for j = 1:num_rep
    disp(j);
%     X = randn(n, p) / sqrt(n);
    X = (rand(n, p) * 2 - 1) * sqrt(3) / sqrt(n);
    W = orth(randn(n, k));
    U = orth(randn(p, k));
    signal_mat = W * D * U';
    Y = signal_mat + X;
    [cos(1, :, j), lambda(1, :, j)] = sketchingMethods.hada(Y, U, r);
    [cos(2, :, j), lambda(2, :, j)] = sketchingMethods.unif(Y, U, r);
    [cos(3, :, j), lambda(3, :, j)] = sketchingMethods.coun(Y, U, r, false);
    [cos(4, :, j), lambda(4, :, j)] = sketchingMethods.coun(Y, U, r, true);
    [cos(5, :, j), lambda(5, :, j)] = sketchingMethods.gaus(Y, U, r);
end
if ~exist('results/', 'dir')
       mkdir('results/')
end
if ~exist('plots/', 'dir')
       mkdir('plots/')
end
% write results
for i = 1:5
    filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
    csvwrite(strcat('results/uniX2_cos_', names(i), filename), squeeze(cos(i, :, :)))
    csvwrite(strcat('results/uniX2_lambda_', names(i), filename), squeeze(lambda(i, :, :)))
end
% read results
for i = 1:5
    filename = sprintf('_n_%d_p_%d_r_%d_k_%d_nrep_%d.csv', n, p, r, k, num_rep);
    cos(i, :, :) = csvread(strcat('results/uniX2_cos_', names(i), filename));
    lambda(i, :, :) = csvread(strcat('results/uniX2_lambda_', names(i), filename));
end

%% plot
% cos
figure, hold on;
mark = {':', '-', ':', '-', '--', '-.', '-.', '--', '-.'};
for i = 1:5
    errorbar(d, mean(cos(i, :, :), 3), std(cos(i, :, :), 0, 3),'lineWidth', 2, 'DisplayName', names(i), 'linestyle', mark{i});
end
legend('location','southeast');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
ylim([0.9, 1]);
set(gca,'fontsize',25)
grid on;
filename = sprintf('plots/final/cos_all_n_%d_p_%d_r_%d_k_%d_nrep_%d.pdf', n, p, r, k, num_rep);
saveTightFigure(gcf, filename);
close(gcf);

% lambda
figure, hold on;
for i = 1:8
    errorbar(d .^ 2, mean(lambda(i, :, :), 3), std(lambda(i, :, :), 0, 3),'lineWidth', 3, 'DisplayName', names(i), 'linestyle', ':');
end
plot([0 250], [0, 250], 'lineWidth', 3, 'linestyle', "-", 'DisplayName', "y=x")
legend('location','northwest');
xlabel('$$d_i^2$$', 'Interpreter', 'LaTex');
ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize', 14)
grid on;
filename = sprintf('lambda_all_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
saveas(gcf, filename);


