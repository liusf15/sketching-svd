% sketching methods
classdef sketchingMethods
    methods(Static)
        % orthogonal projection
        function [cos, lambda] = orth(Y, U, r)
            [n, ~] = size(Y);
            [~, k] = size(U);
            S = orth(randn(n, r))'; % S: r * n
            Y_t = S * Y;
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1: k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end
        % formula for orthogonal sketching
        function [theory_cos, theory_lambda] = theory_orth(gamma, xi, d)
            k = length(d);
            theory_cos = zeros(k, 1);
            theory_lambda = zeros(k, 1);
            for i = 1:k
                theory_cos(i) = max((xi - gamma / d(i)^4) / (xi + gamma / d(i)^2),0);
                theory_lambda(i) = max( (1 + d(i)^2) * ( xi + gamma / d(i)^2), ...
                    (1+sqrt(gamma/xi))*(sqrt(gamma)+sqrt(xi))*sqrt(xi) );
            end
        end
        % Gaussian projection
        function [cos, lambda] = gaus(Y, U, r)
            [n, ~] = size(Y);
            [~, k] = size(U);
            S = randn(r, n) / sqrt(n); % S: r * n
            Y_t = S * Y;
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1 : k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end
        % formula for Gaussian projection
        function [theory_cos_gaus, theory_lambda_gaus, theory_cos_gaus_approx, theory_lambda_gaus_approx] = theory_gaus(gamma, xi, d)
            k = length(d);
            theory_cos_gaus = zeros(k, 1);
            theory_lambda_gaus = zeros(k, 1);
            theory_cos_gaus_approx = zeros(k, 1);
            theory_lambda_gaus_approx = zeros(k, 1);
            for i = 1:k
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
        end
        % uniform sampling
        function [cos, lambda] = unif(Y, U, r)
            [n, ~] = size(Y);
            [~, k] = size(U);
            sampled_ind = binornd(1, r / n, n, 1);
            Y_t = Y(sampled_ind == 1, :);
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1 : k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end
        % formula for uniform sampling
        function [theory_cos_unif, theory_lambda_unif] = theory_unif(gamma, xi, d)
            k = length(d);
            theory_cos_unif = zeros(k, 1);
            theory_lambda_unif = zeros(k, 1);
            for i = 1:k
                theory_cos_unif(i) = max((xi - gamma / d(i)^4) / (xi + gamma / d(i)^2),0);
                theory_lambda_unif(i) = max( (1 + d(i)^2) * ( xi + gamma / d(i)^2), (1+sqrt(gamma/xi))*(sqrt(gamma)+sqrt(xi))*sqrt(xi) );
            end
        end
        % SRHT
        function [cos, lambda] = hada(Y, U, r)
            [n, ~] = size(Y);
            n = pow2(floor(log2(n)));
            Y = Y(1:n, :);
            [~, k] = size(U);
            Di = 2 * binornd(1, 1 / 2, n, 1) - 1;
            Y_t = diag(Di) * Y;
            H = hadamard(n);
            Y_t = H * Y_t / sqrt(n);
            sampled_ind = binornd(1, r / n, n, 1);
            Y_t = Y_t(sampled_ind == 1, :);
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1 : k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end
        % formula for SRHT
        function [theory_cos_hada, theory_lambda_hada] = theory_hada(gamma, xi, d)
            k = length(d);
            theory_cos_hada = zeros(k, 1);
            theory_lambda_hada = zeros(k, 1);
            for i = 1:k
                theory_cos_hada(i) = max((xi - gamma / d(i)^4) / (xi + gamma / d(i)^2),0);
                theory_lambda_hada(i) = max( (1 + d(i)^2) * ( xi + gamma / d(i)^2), (1+sqrt(gamma/xi))*(sqrt(gamma)+sqrt(xi))*sqrt(xi) );
            end
        end
        % Count sketch
        function [cos, lambda] = coun(Y, U, r, normalized)
            [n, p] = size(Y);
            [~, k] = size(U);
            ra = randi([1 r], n, 1);
            count = zeros(r, 1);
            for ii = 1 : r
                count(ii) = sum(ra == ii);
            end

            sig = 2 * randi([0 1] , n, 1) - 1;
            Y_t = zeros(r, p); 

            for ii = 1:n
                Y_t(ra(ii), :) = Y_t(ra(ii), : ) + sig(ii) * Y(ii, : );
            end
            %normalized count-sketch
            if normalized
                ind_nonempty = (count > 0);   
                Y_t = Y_t(ind_nonempty, :);
                count = count(ind_nonempty);
                normalizing_mx = (sqrt(count) * ones(1, p));
                %normalizing_mx = n/r*(ones(length(count),1)*ones(1,p));
                Y_t = Y_t ./ normalizing_mx; 
            end
            
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1 : k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end        
        % formula for CountSketch
        function [theory_cos_coun, theory_lambda_coun] = theory_coun(gamma, xi, d)
            k = length(d);
            xi_eff = xi*(1-exp(-1/xi));
            theory_cos_coun = zeros(k, 1);
            theory_lambda_coun = zeros(k, 1);
            for i = 1:k
                theory_cos_coun(i) = max((xi_eff - gamma / d(i)^4) / (xi_eff + gamma / d(i)^2),0);
                theory_lambda_coun(i) = max( (1 + d(i)^2) * ( xi_eff + gamma / d(i)^2), (1+sqrt(gamma/xi_eff))*(sqrt(gamma)+sqrt(xi_eff))*sqrt(xi_eff) );
            end
        end
        % leverage sampling
        function [cos, lambda] = leve(Y, U, r)
            [n, ~] = size(Y);
            [~, k] = size(U);
            Hat = Y * ((Y' * Y) \ Y');
            leverage_scores = diag(Hat);
            leverage_prob = min(leverage_scores ./ sum(leverage_scores) * r, 1);
            sampled_ind = binornd(1, leverage_prob, n, 1);
            Y_t = Y(sampled_ind == 1, :);
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1 : k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end
        % osnap
        function [cos, lambda] = osna(Y, U, r, s)
            [n, ~] = size(Y);
            [~, k] = size(U);
            S = zeros(r, n);
            for i = 1:r
                ind = randsample(r, s);
                signs = 2 * binornd(1, 0.5, s, 1) - 1;
                S(i, ind) = signs ./ sqrt(s);
            end
            Y_t = S * Y;
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1 : k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end
    end
end