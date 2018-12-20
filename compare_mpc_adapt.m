
Time_out = 20;
seed = 5000;

[y_adapt, x_adapt, Ck_adapt, error, yhat_adapt, xhat_adapt, a_init, true_params] = adapt_ctrl(Time_out, seed);
[ympc, xmpc, Ck, yhat_mpc, xhat_mpc] = mpc_ctrl(Time_out, seed);

%%
base_line_error = sum((true_params - a_init).^2);

figure
plot(linspace(0, Time_out, length(error)), error, 'b')
hold on 
plot([0, Time_out], [sum((true_params - a_init).^2), sum((true_params - a_init).^2)], 'r')

grid on

error_adapt = sum(sum(y_adapt-yhat_adapt, 2).^2, 1);
error_mpc = sum(sum(ympc-yhat_mpc, 2).^2, 1);
%%

error_Xadapt = sum(sum((x_adapt-xhat_adapt).^2, 2), 1);
error_Xmpc = sum(sum((xmpc-xhat_mpc).^2, 2), 1);

cum_error_Xadapt = sum(cumsum((x_adapt-xhat_adapt).^2, 2), 1);
cum_error_Xmpc = sum(cumsum((xmpc-xhat_mpc).^2, 2), 1);

%%

figure
plot(cum_error_Xadapt, 'b')
hold on 
grid on
plot(cum_error_Xmpc, 'r')

title('Cumlative State Error adapt vs mpc')
legend('Adaptive Ctrl', 'MPC')



