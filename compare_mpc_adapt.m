
Time_out = 25;
seed = 100;
adv_true = 1;
%% Adaptive MPC using random sampling

[y_adapt, x_adapt, Ck_adapt, error, yhat_adapt, xhat_adapt, a_init, true_params,...
    u_adpt, a_mat, params_mat] = adapt_ctrl(Time_out, seed, adv_true);
%% Adaptive GP
[y_agp, x_agp, Ck_agp, errorgp, yhat_gp, xhat_gp, a_init_gp, true_params_gp,...
    u_agp, a_mat_gp, params_mat_gp] = adapt_gp(Time_out, seed, adv_true);

%%
[ympc, xmpc, Ck, yhat_mpc, xhat_mpc, true_params_mpc, u_mpc] = mpc_ctrl(Time_out, seed, adv_true);

%%
base_line_error = sum((true_params - a_init).^2);

figure
plot(linspace(0, Time_out, length(error)), error, 'b')
hold on 
plot([0, Time_out], [sum((true_params - a_init).^2), sum((true_params - a_init).^2)], 'r')

grid on

error_adapt = sum(sum(y_adapt-yhat_adapt, 2).^2, 1);
error_mpc = sum(sum(ympc-yhat_mpc, 2).^2, 1);
title('Model Parameter Error')

%% State Estimation Error

error_Xadapt = sum(sum((x_adapt-xhat_adapt).^2, 2), 1);
error_Xmpc = sum(sum((xmpc-xhat_mpc).^2, 2), 1);

cum_error_Xadapt_sq = sum(cumsum((x_adapt-xhat_adapt).^2, 2), 1);
cum_error_Xmpc_sq = sum(cumsum((xmpc-xhat_mpc).^2, 2), 1);
cum_error_Xagp_sq = sum(cumsum((x_agp-xhat_gp).^2, 2), 1);

cum_time = [1:1:length(cum_error_Xadapt_sq)-1];

MSE_adapt = cum_error_Xadapt_sq(2:end)./cum_time;
MSE_mpc = cum_error_Xmpc_sq(2:end)./cum_time;
MSE_agp = cum_error_Xagp_sq(2:end)./cum_time;


cum_error_Xadapt = sum(cumsum((x_adapt-xhat_adapt), 2), 1);
cum_error_Xmpc = sum(cumsum((xmpc-xhat_mpc), 2), 1);

c_theta_adpt = sum(cumsum((x_adapt(2, :)-xhat_adapt(2, :)).^2, 2), 1);
c_theta_mpc = sum(cumsum((xmpc(2, :)-xhat_mpc(2, :)).^2, 2), 1);

c_xpos_adpt = sum(cumsum((x_adapt(1, :)-xhat_adapt(1, :)).^2, 2), 1);
c_xpos_mpc = sum(cumsum((xmpc(1, :)-xhat_mpc(1, :)).^2, 2), 1);

% do a polyfit

x = linspace(0, length(cum_error_Xadapt_sq), length(cum_error_Xadapt_sq));
coefs_adpt = polyfit(x, cum_error_Xadapt_sq, 1);
coefs_mpc = polyfit(x, cum_error_Xmpc_sq, 1);

x_inputs = linspace(0, length(cum_error_Xadapt_sq), 10);


%% Optimality

J_adapt = cumsum(u_adpt.^2) + cumsum(sum(x_adapt(1:2, :).^2, 1));
J_mpc = cumsum(u_mpc.^2) + cumsum(sum(xmpc(1:2, :).^2, 1));
J_gp = cumsum(u_agp.^2) + cumsum(sum(x_agp(1:2, :).^2, 1));


%% cum state error

figure
plot(cum_error_Xadapt_sq, 'b')
hold on 
grid on
plot(cum_error_Xmpc_sq, 'r')
hold on
plot(cum_error_Xagp_sq, 'k')
plot(x_inputs, coefs_adpt*[x_inputs; ones(1, 10)], 'k+')
plot(x_inputs, coefs_mpc*[x_inputs; ones(1, 10)], 'k+')

title('Cumlative State Error Squared adapt vs mpc')
legend('Adaptive Ctrl', 'MPC', 'Adaptive GP')
xlabel('Time Step: 100Hz Controller')

figure
plot(MSE_adapt, 'b')
hold on 
grid on
plot(MSE_mpc, 'r')
hold on
plot(MSE_agp, 'k')

title('Mean-Squared State Error adapt vs mpc')
legend('Adaptive Ctrl', 'MPC', 'Adaptive GP')
xlabel('Time Step: 100Hz Controller')

figure
plot(J_adapt, 'b')
hold on 
grid on
plot(J_mpc, 'r')
plot(J_gp, 'k')
title('Cumulative Input Optimality')
legend('Adaptive Ctrl', 'MPC', 'Adaptive GP')
xlabel('Time Step: 100Hz Controller')

%%

figure
plot(cum_error_Xadapt, 'b')
hold on 
grid on
plot(cum_error_Xmpc, 'r')

title('Cumlative State Error adapt vs mpc')
legend('Adaptive Ctrl', 'MPC')

%%
figure
plot(c_theta_adpt, 'b')
hold on 
grid on
plot(c_theta_mpc, 'r')

title('Cumlative State Error Squared \theta adapt vs mpc')
legend('Adaptive Ctrl', 'MPC')

figure
plot(c_xpos_adpt, 'b')
hold on 
grid on
plot(c_xpos_mpc, 'r')

title('Cumlative State Error Squared Cart Position adapt vs mpc')
legend('Adaptive Ctrl', 'MPC')

%%

figure
loglog(cum_error_Xadapt_sq, 'b')
grid on
hold on
loglog(cum_error_Xmpc_sq, 'r')

title('Log Log Cumlative State Error Squared')
legend('Adaptive Ctrl', 'MPC')

% %%
% figure
% 
% hold on
% % plot(sum((ympc-yhat_mpc).^2, 1), 'r');
% plot(sum((ympc(2,:)-yhat_mpc(2, :)).^2, 1)-sum((y_adapt(2, :)-yhat_adapt(2, :)).^2, 1), 'b');
% 
% grid on

%%

figure
title('How the Model Parameters change')

subplot(2, 2, 1)
plot(a_mat(:, 1), 'b+-')
grid on

subplot(2, 2, 2)
plot(a_mat(:, 2), 'r+-')
grid on

subplot(2, 2, 3)
plot(a_mat(:, 3), 'm+-')
grid on

subplot(2, 2, 4)
plot(a_mat(:, 4), 'k+-')

grid on

figure
title('How the Model Parameters change GP')

subplot(2, 2, 1)
plot(a_mat_gp(:, 1), 'b+-')
grid on

subplot(2, 2, 2)
plot(a_mat_gp(:, 2), 'r+-')
grid on

subplot(2, 2, 3)
plot(a_mat_gp(:, 3), 'm+-')
grid on

subplot(2, 2, 4)
plot(a_mat_gp(:, 4), 'k+-')

grid on

figure
title('How the Model Parameters change GP')

subplot(2, 2, 1)
plot(a_mat_gp(:, 1), 'b')
hold on
plot(a_mat(1:length(a_mat(:, 1))/length(a_mat_gp(:, 1)):end, 1), 'r')
grid on

subplot(2, 2, 2)
plot(a_mat_gp(:, 2), 'b')
hold on
plot(a_mat(1:length(a_mat(:, 1))/length(a_mat_gp(:, 1)):end, 2), 'r')
grid on

subplot(2, 2, 3)
plot(a_mat_gp(:, 3), 'b')
hold on
plot(a_mat(1:length(a_mat(:, 1))/length(a_mat_gp(:, 1)):end, 3), 'r')
grid on

subplot(2, 2, 4)
plot(a_mat_gp(:, 4), 'b')
hold on
plot(a_mat(1:length(a_mat(:, 1))/length(a_mat_gp(:, 1)):end, 4), 'r')

grid on

%% 

figure

plot(y_adapt(1, :), 'b')
hold on
plot(ympc(1, :), 'r')
plot(y_agp(1, :), 'k')
grid on


title('Outputs: Displacements of Cart and Angle of Pendulum')
legend('Adaptive Ctrl', 'MPC', 'Adaptive Bayesian Optimisation')


