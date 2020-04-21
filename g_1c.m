function [y] = g_1c(m, gamma, xi)
    y = - gamma / m + xi / m * (1 - m_s_1c(-1 / m, xi) / m);
end