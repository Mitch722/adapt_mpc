function hyper_params = truncate_gauss3(mu, sigma, params_set, no_models)

% begin the hyper parameters
hyper_params = mvnrnd(mu, sigma, no_models);
% low 
rep_mu = repmat(mu, no_models, 1);

rep_params_low = repmat(params_set(:, 1)', no_models, 1);

index_low = hyper_params < rep_params_low;

hyper_params(index_low) = rep_mu(index_low);


rep_params_high = repmat(params_set(:, 2)', no_models, 1);

index_high = hyper_params > rep_params_high;

hyper_params(index_high) = rep_mu(index_high);


