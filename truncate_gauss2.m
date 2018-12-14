function hyper_params = truncate_gauss2(mu, sigma, params_set, no_models)

hyper_params = mvnrnd(mu, sigma, no_models);

rep_params_low = repmat(params_set(:, 1)', no_models, 1);

index_low = hyper_params < rep_params_low;

hyper_params(index_low) = rep_params_low(index_low);


rep_params_high = repmat(params_set(:, 2)', no_models, 1);

index_high = hyper_params > rep_params_high;

hyper_params(index_high) = rep_params_high(index_high);


