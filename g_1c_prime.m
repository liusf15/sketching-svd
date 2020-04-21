function [y] = g_1c_prime(m, gamma, xi)
    y = (gamma - xi) / m ^ 2 + 2 * xi / m ^ 3 * m_s_1c(-1 / m, xi) - xi / m ^ 4 * m_s_1c_prime(-1 / m, xi);
end
    