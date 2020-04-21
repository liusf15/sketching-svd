function [y] = m_s_2c(z, xi)
    lambda_p = (1 + sqrt(xi)) ^ 2;
    lambda_n = (1 - sqrt(xi)) ^ 2;
    y = (-z -1 + xi + sqrt((z - lambda_p) * (z - lambda_n))) / (2 * z);
end