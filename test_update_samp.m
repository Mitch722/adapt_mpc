mu_arr = [0.5, 0.6, 0.8, 0.1];
mu_arr = [mu_arr, mu_arr];

std_arr = ones(1, length(mu_arr));

pd = makedist('Normal');



