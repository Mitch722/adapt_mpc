function [mat_k, mat_k_s, mat_k_ss] = matern_kernel52(x_data, x_pred, sigF, sigN, l)


[d_kern_sq, dstar_sq, dstarstar_sq] = kernels(x_data, x_pred);

d_k = d_kern_sq.^0.5;
d_star = dstar_sq.^0.5;
d_starstar = dstarstar_sq.^0.5;

mat_k = sigF^2 * (1 + 5^0.5 / l * d_k + 5/3/l^2 * d_k.^2) .* exp(-5^0.5 / l * d_k) + sigN^2*eye(size(d_kern_sq));

mat_k_s = sigF^2 * (1 + 5^0.5 / l * d_star + 5/3/l^2 * d_star.^2) .* exp(-5^0.5 / l * d_star);

mat_k_ss = sigF^2 * (1 + 5^0.5 / l * d_starstar + 5/3/l^2 * d_starstar.^2) .* exp(-5^0.5 / l * d_starstar) + sigN^2*eye(size(d_starstar));

end