function [rbf_k, rbf_kstar, rbf_kstarstar] = rbf_kernel(x_data, x_pred, sigF, sigN, l)


[k_kern, kstar, kstarstar] = kernels(x_data, x_pred);

rbf_k = sigF^2*exp(-k_kern/(2*l^2)) + sigN^2*eye(size(k_kern));

rbf_kstar = sigF^2*exp(-kstar/(2*l^2)) + sigN^2*eye(size(kstar));

rbf_kstarstar = sigF^2*exp(-kstarstar/(2*l^2)) + sigN^2*eye(size(kstarstar));

end