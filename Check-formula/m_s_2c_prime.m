function [y] = m_s_2c_prime(z, xi)
    lambda_p = (1 + sqrt(xi)) ^ 2;
    lambda_n = (1 - sqrt(xi)) ^ 2;
    y1 = (-1 + (z - xi - 1) / sqrt((z - lambda_p) * (z - lambda_n))) / (2 * z);
    y2 = (z + 1 - xi - sqrt((z - lambda_p) * (z - lambda_n))) / (2 * z ^ 2);
    y = y1 + y2;
end