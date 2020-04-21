function [y] = gauss_cos_approx(d, gamma, xi)
    y = (xi - (1 + xi) * gamma / d ^ 4) / (xi + (1 + xi) * gamma / d ^ 2);
end