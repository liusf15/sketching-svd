
% multi spikes

cd('/Users/sifanliu/Dropbox/Projects/Sketching SVD (1)/Experiments/Exp1-Check_formulas');
%% parameters
n = 4000;
p = 3200;
r = 2500;
gamma = p / n;
xi = r / n;

%% 3.1 Orthogonal projection
num_rep = 5;
k = 10;
res = zeros(num_rep, k);
lambda = zeros(num_rep, k);
theory_cos = zeros(k, 1);
theory_lambda = zeros(k, 1);
rng(1);
d = sort(rand(1, k) * 20 + 1, 'descend');
D = diag(d);

for j = 1:num_rep
    disp(j);

    X = randn(n, p) / sqrt(n);
    S = orth(randn(n, r))'; % S: r * n

    W = orth(randn(n, k));
    U = orth(randn(p, k));

    signal_mat = W * D * U';
    Y = signal_mat + X;
    Y_t = S * Y;
    [V_t, D_t, U_t] = svd(Y_t,'econ');

    for i = 1:k
        res(j, i) = (U_t(:,i)'*U(:,i))^2;
        lambda(j, i) = D_t(i, i) ^ 2;
    end
end

% theory
for i = 1:k
    theory_cos(i) = max((xi - gamma / d(i)^4) / (xi + gamma / d(i)^2),0);
    theory_lambda(i) = max( (1 + d(i)^2) * ( xi + gamma / d(i)^2), ...
        (1+sqrt(gamma/xi))*(sqrt(gamma)+sqrt(xi))*sqrt(xi) );
end


%%% plot cos
savefigs = 1; 
a = {'-', '--', '-.', ':'};
figure, hold on

me = mean(res, 1);
re = std(res, 1);
h1 = errorbar(d, me, re,'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = plot(d, theory_cos, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;
if savefigs==1
    filename = sprintf('multispike_cos_orth_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi, num_rep, n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%%% plot lambda
savefigs =1;
figure, hold on
me = mean(lambda, 1);
re = std(lambda, 1);
h1 = errorbar(d, me, re, 'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
set(h1,'LineStyle',a{1});
h2 = plot(d, theory_lambda, 'DisplayName', 'Theory', 'linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
set(h2,'LineStyle',a{2});
xlabel('$$d_i$$', 'Interpreter', 'LaTex')
ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;

if savefigs==1
    filename = sprintf('multispike_lambda_orth_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%% 3.3 uniform sampling
n = 4000;
p = 3200;
r = 2500;
gamma = p / n;
xi = r / n;
num_rep = 5;
k = 10;
rng(1);
d = sort(rand(1, k) * 20 + 1, 'descend');
D = diag(d);
cos_unif = zeros(num_rep, k);
lambda_unif = zeros(num_rep, k);
theory_cos_unif = zeros(k, 1);
theory_lambda_unif = zeros(k, 1);
for j = 1:num_rep
    disp(j);

    X = randn(n, p) / sqrt(n);
    S = orth(randn(n, r))'; % S: r * n

    W = orth(randn(n, k));
    U = orth(randn(p, k));

    signal_mat = W * D * U';
    Y = signal_mat + X;
    sampled_ind = binornd(1,r/n,n,1);
    Y_t = Y(sampled_ind==1,:);
    [V_t, D_t, U_t] = svd(Y_t,'econ');

    for i = 1:k
        cos_unif(j, i) = (U_t(:,i)'*U(:,i))^2;
        lambda_unif(j, i) = D_t(i, i) ^ 2;
    end
end

for i = 1:k
    theory_cos_unif(i) = max((xi - gamma / d(i)^4) / (xi + gamma / d(i)^2),0);
    theory_lambda_unif(i) = max( (1 + d(i)^2) * ( xi + gamma / d(i)^2), (1+sqrt(gamma/xi))*(sqrt(gamma)+sqrt(xi))*sqrt(xi) );
end

%%% plot cos
savefigs =1; a = {'-','--','-.',':'};
figure, hold on

me = mean(cos_unif, 1);
re = std(cos_unif, 1);
h1 = errorbar(d, me, re,'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = plot(d, theory_cos_unif, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;

if savefigs==1
    filename = sprintf('multispike_cos_unif_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%%% plot lambda
savefigs =1; a = {'-','--','-.',':'};
figure, hold on

me = mean(lambda_unif, 1);
re = std(lambda_unif, 1);
h1 = errorbar(d, me, re, 'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = plot(d, theory_lambda_unif, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex')
ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;

if savefigs==1
    filename = sprintf('multispike_lambda_unif_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%% 3.2 Gaussian projection
n = 4000;
p = 3200;
r = 2500;
gamma = p / n;
xi = r / n;
k = 10;
num_rep = 5;
rng(1);
d = sort(rand(1, k) * 20 + 1, 'descend');
D = diag(d);
cos_gaus = zeros(num_rep, k);
lambda_gaus = zeros(num_rep, k);
theory_cos_gaus = zeros(k, 1);
theory_cos_gaus_approx = zeros(k, 1);
theory_lambda_gaus = zeros(k, 1);
theory_lambda_gaus_approx = zeros(k, 1);

for j = 1:num_rep
    disp(j);

    X = randn(n, p) / sqrt(n);
    S = randn(r, n) / sqrt(n); % S: r * n

    W = orth(randn(n, k));
    U = orth(randn(p, k));

    signal_mat = W * D * U';
    Y = signal_mat + X;
    Y_t = S*Y;
    [V_t, D_t, U_t] = svd(Y_t,'econ');

    for i = 1:k
        cos_gaus(j, i) = (U_t(:,i)'*U(:,i))^2;
        lambda_gaus(j, i) = D_t(i, i) ^ 2;
    end
end

for i = 1:k
%     theory_cos_gaus(i) = max((xi - gamma / d(i)^4) / (xi + gamma / d(i)^2), 0); %% TBC
    theory_cos_gaus_approx(i) = max(gauss_cos_approx(d(i), gamma, xi), 0);
    theory_cos_gaus(i) = max(gauss_cos(gamma, xi, d(i)), 0);
    b = gamma / d(i)^2;
    m = - b/((1+b)*(xi+b));
    C1 = m^3;
    C2 = -(1+xi-2*gamma)*m^2-m;
    C3 = -(1-gamma)*(gamma-xi)*m-gamma;
    pol = [C1 C2 C3];
    root = roots(pol);
    theory_lambda_gaus(i) = min(real(root));
    theory_lambda_gaus_approx(i) = xi * d(i) ^ 2 + (xi * gamma + xi + gamma) + ...
        (gamma + xi + 1) * gamma / d(i) ^ 2;
end


%%% plot cos
savefigs =1;
figure, hold on
a = {'-','--','-.',':'};
me = mean(cos_gaus, 1);
re = std(cos_gaus, 1);
h1 = errorbar(d, me, re, 'lineWidth', 3, 'DisplayName', ...
    'Simulation', 'color', [.1, .4, .9]);
% h1.Marker = '*';
% h1.MarkerSize = 10;

h2 = plot(d, theory_cos_gaus, 'DisplayName', 'Theory', ...
    'linewidth', 4, 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(d, theory_cos_gaus_approx, 'LineStyle','--', ...
    'linewidth', 4, 'color', [1, .5, 0], 'DisplayName', 'Theory approx');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;

if savefigs==1
    filename = sprintf('multispike_cos_gaus_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%%% plot lambda
savefigs =1; 
figure, hold on

me = mean(lambda_gaus, 1);
re = std(lambda_gaus, 1);
h1 = errorbar(d, me, re, 'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
% h1.Marker = '*';
% h1.MarkerSize = 10;
h2 = plot(d, theory_lambda_gaus, 'DisplayName', 'Theory','linewidth', 4, 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(d, theory_lambda_gaus_approx, 'DisplayName', 'Theory approx','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex')
ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;

if savefigs==1
    filename = sprintf('multispike_lambda_gaus_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%% 3.2 Gaussian projection one spike
n = 1000;
p = 800;
r = 500;
gamma = p / n;
xi = r / n;
k = 10;
res = zeros(num_rep, k);
lam = zeros(num_rep, k);
theory_cos = zeros(k, 1);
theory_lam = zeros(k, 1);
d = linspace(1,k,k);

for i=1:length(d)
    i
    D = d(i);
    for j = 1:num_rep
        X = randn(n, p) / sqrt(n);
        S = randn(r,n)/sqrt(n); % S: r * n
        
        W = orth(randn(n, 1));
        U = orth(randn(p, 1));
        
        signal_mat = W * D * U';
        Y = signal_mat + X;
        Y_t = S*Y;
        [V_t, D_t, U_t] = svd(Y_t,'econ');
        
        
        res(j, i) = (U_t(:,1)'*U(:,1))^2;
        lam(j,i) = D_t(1,1); %top sketched spiked sval
    end
end
%
for i = 1:k
    %theory_cos(i) = max((xi - gamma / d(i)^4) / (xi + gamma / d(i)^2),0);
    theory_cos(i) = gauss_cos_approx(d(i), gamma, xi);
    b = gamma / d(i)^2;
    m = - b/((1+b)*(xi+b));
    C1 = m^3;
    C2 = -(1+xi-2*gamma)*m^2-m;
    C3 = -(1-gamma)*(gamma-xi)*m-gamma;
    pol = [C1 C2 C3];
    root = roots(pol);
    
    theory_lam(i) = sqrt(min(real(root)));
end

%%% plot cos
savefigs =1;
figure, hold on

me = mean(res, 1);
re = std(res, 1);
h1 = errorbar(linspace(1, k, k), me, re,'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = plot(linspace(1, k, k), theory_cos, 'DisplayName', 'Theory', 'linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('Signal Strengh');
ylabel('$$|\langle u_1,\tilde\xi_1\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;

if savefigs==1
    filename = sprintf('singlespike_cos_gauss_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d.png', gamma, xi,num_rep,n);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%%% plot lambda
savefigs =1; 
figure, hold on

me = mean(lam, 1);
re = std(lam, 1);
h1 = errorbar(linspace(1, k, k), me, re, 'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = plot(linspace(1, k, k), theory_lam, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('Signal strength')
ylabel('$$\tilde\lambda_1$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;

if savefigs==1
    filename = sprintf('singlespike_lambda_gauss_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d.png', gamma, xi,num_rep,n);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end



%% 3.4 Hadamard projection
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

cos_hada = zeros(num_rep, k);
lambda_hada = zeros(num_rep, k);
theory_cos_hada = zeros(k, 1);
theory_lambda_hada = zeros(k, 1);

for j = 1:num_rep
    disp(j);

    X = randn(n, p) / sqrt(n);
    S = randn(r,n)/sqrt(n); % S: r * n

    W = orth(randn(n, k));
    U = orth(randn(p, k));

    signal_mat = W * D * U';
    Y = signal_mat + X;
    %Hadamard
    Di = 2*binornd(1,1/2,n,1)-1;
    Y_t = (Di*ones(1,p)).*Y;
    H = hadamard(n);
    Y_t = H*Y_t/sqrt(n);
    sampled_ind = binornd(1,r/n,n,1);
    Y_t = Y_t(sampled_ind==1,:);
    [V_t, D_t, U_t] = svd(Y_t,'econ');
    
    for i = 1:k
        cos_hada(j, i) = (U_t(:,i)'*U(:,i))^2;
        lambda_hada(j, i) = D_t(i, i) ^ 2;
    end
end

for i = 1:k
    theory_cos_hada(i) = max((xi - gamma / d(i)^4) / (xi + gamma / d(i)^2),0);
    theory_lambda_hada(i) = max( (1 + d(i)^2) * ( xi + gamma / d(i)^2), (1+sqrt(gamma/xi))*(sqrt(gamma)+sqrt(xi))*sqrt(xi) );
end

%%% plot cos
savefigs =1;
figure, hold on

me = mean(cos_hada, 1);
re = std(cos_hada, 1);
h1 = errorbar(d, me, re,'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = plot(d, theory_cos_hada, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;

if savefigs==1
    filename = sprintf('multispike_cos_hada_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%%% plot lambda
savefigs =1; 
figure, hold on

me = mean(lambda_hada, 1);
re = std(lambda_hada, 1);
h1 = errorbar(d, me, re, 'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = plot(d, theory_lambda_hada, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex')
ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;


if savefigs==1
    filename = sprintf('multispike_lambda_hada_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end


%% 3.5 count sketch
n = 4000;
p = 3200;
r = 2500;
gamma = p / n;
xi = r / n;
xi_eff = xi*(1-exp(-1/xi));
k = 10;
num_rep = 5;
rng(1);
d = sort(rand(1, k) * 20 + 1, 'descend');
D = diag(d);
cos_count = zeros(num_rep, k);
lambda_count = zeros(num_rep, k);
cos_count_norm = zeros(num_rep, k);
lambda_count_norm = zeros(num_rep, k);
theory_cos_count = zeros(k, 1);
theory_lambda_count = zeros(k, 1);

for j = 1:num_rep
    disp(j);
    X = randn(n, p) / sqrt(n);

    W = orth(randn(n, k));
    U = orth(randn(p, k));

    signal_mat = W * D * U';
    Y = signal_mat + X;

    ra = randi([1 r],n,1);
    count = zeros(r,1);
    for ii=1:r
        count(ii) = sum(ra==ii);
    end

    sig = 2*randi([0 1],n,1)-1;
    Y_t = zeros(r,p);

    %original count-sketch
    for ii = 1:n
        Y_t(ra(ii),:) = Y_t(ra(ii),:)+sig(ii)*Y(ii,:);
    end

    [~, D_t, U_t] = svd(Y_t,'econ');
    
    for i = 1:k
        cos_count(j, i) = (U_t(:,i)'*U(:,i))^2;
        lambda_count(j, i) = D_t(i, i) ^ 2;
    end

    %normalized count-sketch
    ind_nonempty = (count>0);   
    Y_t2 = Y_t(ind_nonempty,:);
    count = count(ind_nonempty);
    normalizing_mx = (sqrt(count)*ones(1,p));
    %normalizing_mx = n/r*(ones(length(count),1)*ones(1,p));
    Y_t2 = Y_t2./normalizing_mx; 
    [~, D_t2, U_t2] = svd(Y_t2,'econ');
    for i = 1:k
        cos_count_norm(j, i) = (U_t2(:,i)'*U(:,i))^2;
        lambda_count_norm(j, i) = D_t2(i, i) ^ 2;
    end
end

for i = 1:k
    theory_cos_count(i) = max((xi_eff - gamma / d(i)^4) / (xi_eff + gamma / d(i)^2),0);
    theory_lam_count(i) = max( (1 + d(i)^2) * ( xi_eff + gamma / d(i)^2), (1+sqrt(gamma/xi_eff))*(sqrt(gamma)+sqrt(xi_eff))*sqrt(xi_eff) );
end

%%% plot cos
savefigs =1;
figure, hold on

h1 = errorbar(d, mean(cos_count, 1), std(cos_count, 1),'lineWidth', 3, 'DisplayName', 'CS', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = errorbar(d, mean(cos_count_norm, 1), std(cos_count_norm, 1),'lineWidth', 3, 'DisplayName', 'CS Normalized', 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(d, theory_cos_count, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;

if savefigs==1
    filename = sprintf('multispike_cos_count_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%%% plot lambda
savefigs =1; 
figure, hold on
h1 = errorbar(d, mean(lambda_count, 1), std(lambda_count, 1),'lineWidth', 3, 'DisplayName', 'CS', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = errorbar(d, mean(lambda_count_norm, 1), std(lambda_count_norm, 1),'lineWidth', 3, 'DisplayName', 'CS Normalized', 'color', [0, .8, 0], 'linestyle', '-');
h3 = plot(d, theory_lam_count, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex')
ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;


if savefigs==1
    filename = sprintf('multispike_lambda_count_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end


%% 3.6 large signals
n = 4000;
p = 3200;
r = 2500;
gamma = p / n;
xi = r / n;
num_rep = 3;
k = 10;
cos_large = zeros(num_rep, k);
lambda_large = zeros(num_rep, k);
theory_cos_large = zeros(k, 1);
theory_lambda_large = zeros(k, 1);
rng(1);
% d = sort(rand(1, k) * 20 + 5, 'descend') * 20;
d = sqrt(linspace(2, 100, k));
d = sort(d, 'descend');
D = diag(d);
Sigma = eye(p);

% simulation
for j = 1:num_rep
    disp(j);

    X = randn(n, p) * sqrtm(Sigma) / sqrt(n);
    S = orth(randn(n, r))'; % S: r * n

    W = orth(randn(n, k));
    U = orth(randn(p, k));

    signal_mat = W * D * U';
    Y = signal_mat + X;
    Y_t = S * Y;
    [V_t, D_t, U_t] = svd(Y_t, 'econ');

    for i = 1:k
        cos_large(j, i) = (U_t(:,i)'*U(:,i))^2;
        lambda_large(j, i) = D_t(i, i) ^ 2;
    end
end


% theory
E_diag = ones(k, 1);
rho_1 = 1;
rho_2 = 1;
for i = 1:k
    theory_cos_large(i) = (xi - gamma / d(i) ^ 4 * rho_2) / (xi + gamma / d(i) ^ 2 ...
    * (rho_1 + (rho_2 - rho_1 * E_diag(i)) / d(i) ^ 2));
    theory_cos_large(i) = max(theory_cos_large(i), 0);
%     (1 + d(i)^2) * ( xi + gamma / d(i)^2)
%     (1 + sqrt(gamma / xi)) * (sqrt(gamma) + sqrt(xi)) * sqrt(xi)
    theory_lambda_large(i) = xi * d(i) ^ 2 + xi * E_diag(i) + gamma * rho_1;
end


%%% plot cos
savefigs = 1; 
a = {'-', '--', '-.', ':'};
figure, hold on

me = mean(cos_large, 1);
re = std(cos_large, 1);
h1 = errorbar(d, me, re,'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
h2 = plot(d, theory_cos_large, 'DisplayName', 'Theory','linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
xlabel('$$d_i$$', 'Interpreter', 'LaTex');
ylabel('$$|\langle u_i,\tilde\xi_i\rangle|^2$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;
if savefigs==1
    filename = sprintf('multispike_cos_large_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi, num_rep, n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end

%%% plot lambda
savefigs =1;
figure, hold on
me = mean(lambda_large, 1);
re = std(lambda_large, 1);
h1 = errorbar(d, me, re, 'lineWidth', 3, 'DisplayName', 'Simulation', 'color', [.1, .4, .9], 'linestyle', '-');
set(h1,'LineStyle',a{1});
h2 = plot(d, theory_lambda_large, 'DisplayName', 'Theory', 'linewidth', 4, 'color', [1, .5, 0], 'linestyle', '--');
set(h2,'LineStyle',a{2});
xlabel('$$d_i$$', 'Interpreter', 'LaTex')
ylabel('$$\tilde\lambda_i$$', 'Interpreter', 'LaTeX');
set(gca,'fontsize',20)
title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
legend('location','southeast');
grid on;

if savefigs==1
    filename = sprintf('multispike_lambda_large_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end


%% large signals: toeplitz sigma
n = 4000;
p = 3200;
r = 2500;
gamma = p / n;
xi = r / n;
num_rep = 5;
k = 10;
cos_large_toe = zeros(num_rep, k, k);
lambda_large_toe = zeros(num_rep, k);
theory_cos_large_toe = zeros(k, k);
theory_lambda_large_toe = zeros(k, 1);
rng(1);
% d = sort(rand(1, k) * 20 + 5, 'descend') * 20;
d = sqrt(linspace(2, 100, k));
d = sort(d, 'descend');
D = diag(d);
% generate Sigma
q = 0.5;
Sigma_t = toeplitz(q.^linspace(0, p - 1, p));
Sigma_t = Sigma_t * 2; % for q = 0.5
rho_1 = 1;
rho_2 = 1 / (1 - q);
rho_1 = rho_1 * 2;
rho_2 = rho_2 * 2;

% simulation
X = randn(n, p) * sqrtm(Sigma_t) / sqrt(n);
W = orth(randn(n, k));
U = orth(randn(p, k));
signal_mat = W * D * U';
Y = signal_mat + X;
E = U' * Sigma_t * U;

for j = 1:num_rep
    j
    S = orth(randn(n, r))'; % S: r * n
    Y_t = S * Y;
    [V_t, D_t, U_t] = svd(Y_t, 'econ');
    lambda_large_toe(j, :) = diag(D_t(1:k, 1:k)) .^ 2;
    cos_large_toe(j, :, :) = (U_t(:, 1:k)' * U(:, 1:k)) .^ 2;
%     for i = 1:k
%         cos_large_toe(j, i) = (U_t(:,i)'*U(:,i))^2;
%         lambda_large_toe(j, i) = D_t(i, i);
%     end
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
if savefigs==1
    filename = sprintf('cos_toep_q_%0.1f_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', q, gamma, xi, num_rep, n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
end

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

if savefigs==1
    filename = sprintf('lambda_toep_q_%0.1f_gamma_%0.1f_xi_%0.2f_nrep_%d_n_%d_k_%d.png', q, gamma, xi,num_rep,n, k);
    saveas(gcf, filename);
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end









