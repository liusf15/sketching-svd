function [y] = gauss_lambda(gamma, xi, d)
    m = - (gamma / d ^ 2) / (1 + gamma / d ^ 2) / (xi + gamma / d ^ 2);
    y = g_1c(m, gamma, xi);
end