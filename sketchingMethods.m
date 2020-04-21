% sketching methods
classdef sketchingMethods
    methods(Static)
        % orthogonal projection
        function [cos, lambda] = orth(W, D, U, X, r)
            [n, k] = size(W);
            signal_mat = W * D * U';
            Y = signal_mat + X;
            S = orth(randn(n, r))'; % S: r * n
            Y_t = S * Y;
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1: k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end
        % Gaussian projection
        function [cos, lambda] = gaus(W, D, U, X, r)
            [n, k] = size(W);
            signal_mat = W * D * U';
            Y = signal_mat + X;
            S = randn(r, n) / sqrt(n); % S: r * n
            Y_t = S * Y;
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1 : k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end
        % uniform sampling
        function [cos, lambda] = unif(W, D, U, X, r)
            [n, k] = size(W);
            signal_mat = W * D * U';
            Y = signal_mat + X;
            sampled_ind = binornd(1, r / n, n, 1);
            Y_t = Y(sampled_ind == 1, :);
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1 : k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end
        % Hadamard projection
        function [cos, lambda] = hada(W, D, U, X, r)
            [n, k] = size(W);
            signal_mat = W * D * U';
            Y = signal_mat + X;
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
        % Count sketch
        function [cos, lambda] = coun(W, D, U, X, r, normalized)
            [n, k] = size(W);
            [p, ~] = size(U);
            signal_mat = W * D * U';
            Y = signal_mat + X;
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
        % leverage sampling
        function [cos, lambda] = leve(W, D, U, X, r)
            [n, k] = size(W);
            signal_mat = W * D * U';
            Y = signal_mat + X;
            Hat = Y * inv(Y' * Y) * Y';
            leverage_scores = diag(Hat);
            leverage_prob = min(leverage_scores ./ sum(leverage_scores) * r, 1);
            sampled_ind = binornd(1, leverage_prob, n, 1);
            Y_t = Y(sampled_ind == 1, :);
            [~, D_t, U_t] = svd(Y_t, 'econ');
            cos = diag((U_t(:, 1 : k))' * U) .^ 2;
            lambda = diag(D_t(1 : k, 1 : k)) .^ 2;
        end
        % osnap
        function [cos, lambda] = osna(W, D, U, X, r, s)
            [n, k] = size(W);
            signal_mat = W * D * U';
            Y = signal_mat + X;
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