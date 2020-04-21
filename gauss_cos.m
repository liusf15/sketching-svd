function [y] = gauss_cos(gamma, xi, d)
    m = -(gamma / d ^ 2) / (1 + gamma / d ^ 2) / (xi + gamma / d ^ 2);
%     C1 = m^3;
%     C2 = -(1+xi-2*gamma)*m^2-m;
%     C3 = -(1-gamma)*(gamma-xi)*m-gamma;
%     pol = [C1 C2 C3];
%     root = roots(pol);
%     theta = sqrt(min(real(root)));
    alpha = m;
    y = alpha ^ 4 / d ^ 2 * g_1c_prime(alpha, gamma, xi) / (m_s_2c_prime(- 1 / alpha, xi) - alpha ^ 2 * (1 + gamma / d ^ 2));
end
    