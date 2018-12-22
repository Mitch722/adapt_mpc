
Time_out = 30;
seed = 100;
adv_true = 1;

[y_adapt, x_adapt, Ck_adapt, error, yhat_adapt, xhat_adapt, a_init, true_params,...
    u_adpt] = adapt_ctrl(Time_out, seed, adv_true);
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

cum_error_Xadapt = sum(cumsum((x_adapt-xhat_adapt), 2), 1);
cum_error_Xmpc = sum(cumsum((xmpc-xhat_mpc), 2), 1);

c_theta_adpt = sum(cumsum((x_adapt(2, :)-xhat_adapt(2, :)).^2, 2), 1);
c_theta_mpc = sum(cumsum((xmpc(2, :)-xhat_mpc(2, :)).^2, 2), 1);

c_xpos_adpt = sum(cumsum((x_adapt(1, :)-xhat_adapt(1, :)).^2, 2), 1);
c_xpos_mpc = sum(cumsum((xmpc(1, :)-xhat_mpc(1, :)).^2, 2), 1);


%% Optimality

J_adapt = cumsum(u_adpt.^2) + cumsum(sum(x_adapt(1:2, :).^2, 1));
J_mpc = cumsum(u_mpc.^2) + cumsum(sum(xmpc(1:2, :).^2, 1));


%%

figure
plot(cum_error_Xadapt_sq, 'b')
hold on 
grid on
plot(cum_error_Xmpc_sq, 'r')

title('Cumlative State Error Squared adapt vs mpc')
legend('Adaptive Ctrl', 'MPC')

figure
plot(J_adapt, 'b')
hold on 
grid on
plot(J_mpc, 'r')
title('Cumulative Input Optimality')
legend('Adaptive Ctrl', 'MPC')

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




