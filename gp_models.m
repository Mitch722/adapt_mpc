function hyper_params = gp_models(mu, y_data, x_data, params_set, no_models)

npred = 10*length(x_data(:, 1));
x_pred = zeros(npred, 4);
mu_pred = mu*ones(length(x_pred(:, 1)), 1);

mu = mu*ones(length(x_data(:, 1)), 1);

[~, dim] = size(x_data);

for j = 1 : dim

    x_pred(:, j) = params_set(j, 1) + (params_set(j, 2) - params_set(j, 1))*rand(npred, 1);
    % x_pred(:, j) = linspace(params_set(j, 1), params_set(j, 2), npred);
end

[rbf_k, rbf_kstar, rbf_kstarstar] = rbf_kernel(x_data, x_pred, 0.1, 0.1, 1);

inv_K = inv(rbf_k + 0.00001*eye(size(rbf_k)));

y_pred = mu_pred + rbf_kstar * inv_K *(y_data - mu);

sigStar = rbf_kstarstar - rbf_kstar * inv_K * rbf_kstar';

% find the minima of the y_pred

hyper_params = zeros(no_models, 4);

[~, min_index] = min(y_pred);
y_pred_update = y_pred;

[~, max_index] = max(diag(sigStar));
sigStar_up = sigStar;


for j = 1 : round(no_models/2)
    
    hyper_params(j, :) = x_pred(min_index, :);
    
    y_pred_update(min_index, :) = [];
    [~, min_index] = min(y_pred_update);
    
end

for j = round(no_models/2)+1 : no_models
    
    hyper_params(j, :) = x_pred(max_index, :);
    
    sigStar_up(max_index, :) = [];
    [~, max_index] = max(diag(sigStar_up));
    
end

end