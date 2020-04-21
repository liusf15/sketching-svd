% compare all sketching methods
n = 4096;
p = 3200;
r = 2500;
gamma = p / n;
xi = r / n;
k = 10;
rng(1);
d = sort(rand(1, k) * 20 + 1, 'descend');
D = diag(d);

num_rep = 5;
names = ["orthogonal projection", "Gaussian projection", "uniform sampling", "Hadamard", "countSketch", "countSketch-normalized", "leverage", "osnap"];
cos = zeros(length(names), k, num_rep);
lambda = zeros(length(names), k, num_rep);
%% simulation
for j = 1:num_rep
    disp(j);
    X = randn(n, p) / sqrt(n);
    W = orth(randn(n, k));
    U = orth(randn(p, k));
    [cos(1, :, j), lambda(1, :, j)] = sketchingMethods.orth(W, D, U, X, r);
    [cos(2, :, j), lambda(2, :, j)] = sketchingMethods.gaus(W, D, U, X, r);
    [cos(3, :, j), lambda(3, :, j)] = sketchingMethods.unif(W, D, U, X, r);
    [cos(4, :, j), lambda(4, :, j)] = sketchingMethods.hada(W, D, U, X, r);
    [cos(5, :, j), lambda(5, :, j)] = sketchingMethods.coun(W, D, U, X, r, false);
    [cos(6, :, j), lambda(6, :, j)] = sketchingMethods.coun(W, D, U, X, r, true);
    [cos(7, :, j), lambda(7, :, j)] = sketchingMethods.leve(W, D, U, X, r);
    [cos(8, :, j), lambda(8, :, j)] = sketchingMethods.osna(W, D, U, X, r, 500);
end

csvwrite("compare_all_cos.csv", cos)
csvwrite("compare_all_lambda.csv", lambda)

%% plot
% cos
figure, hold on;
mark = {':', '-', ':', '-', '--', '-.', '-.', '--', '-.'};
for i = 1:8
    errorbar(d, mean(cos(i, :, :), 3), std(cos(i, :, :), 0, 3),'lineWidth', 2, 'DisplayName', names(i), 'linestyle', mark{i});
end
legend('location','southeast');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
ylim([0, 1]);
set(gca,'fontsize',14)
% title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
grid on;
filename = sprintf('cos_all_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
saveas(gcf, filename);

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
% title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
grid on;
filename = sprintf('lambda_all_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
saveas(gcf, filename);


