
m = 1;
M = 4;
l = 1;
g = 9.81;

Ts = 0.01;

a11 = m*g/M;
a22 = (M+m)*g/(M*l);

b11 = 1/M;
b22 = 1/(M*l);

var = 0.4;

a11_set = [a11-var*a11, a11+var*a11];
a22_set = [a22-var*a22, a22+var*a22];

b11_set = [b11-var*b11, b11+var*b11];
b22_set = [b22-var*b22, b22+var*b22];

params_set = [a11_set; a22_set; b11_set; b22_set];

load('advs_s_100.mat')

ndp = 100;

y_data = 100*rand(ndp, 1);
mu = mean(y_data);
x_data = a_mat(1:ndp, :);

no_models = 10;

hyper_params_gp = gp_models(mu, y_data, x_data, params_set, no_models);


