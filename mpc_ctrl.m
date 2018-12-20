function [y, x, Ck, yhat, xhat] = mpc_ctrl(Time_out, seed)

%% initial estimated values
% set the seed
rng(seed);

m = 1;
M = 4;
l = 1;
g = 9.81;

Ts = 0.01;
%%
a11 = m*g/M;
a22 = (M+m)*g/(M*l);

b11 = 1/M;
b22 = 1/(M*l);
%%

sysd = makesysd_a(a11, a22, b11, b22, Ts);

Q = sysd.C'*sysd.C;
Q(1, 1) = 0.1;
R = 1;

[P, L_eig, Kopt] = dare(sysd.A, sysd.B, Q, R);
%%

Phi = sysd.A - sysd.B*Kopt;

multiple = 0.999;
% Put poles near orgin of unit circle in z-domain
%obvs_poles = multiple * real(L_eig);
obvs_poles = [0.0001, 0.001, 0.002, 0.0002]';

L = place(sysd.A', sysd.C', obvs_poles);
L = L';

%% State Estimator Dynamics

Aobv = sysd.A - sysd.B*Kopt;
    
Bobv = sysd.B;
         
Cobv = sysd.C;

Dobv = [0; 0];


%%
[~, no_states] = size(Aobv);
[no_outputs, ~] = size(Cobv);

%% set of m, M, l
var = 0.4;

a11_set = [a11-var*a11, a11+var*a11];
a22_set = [a22-var*a22, a22+var*a22];

b11_set = [b11-var*b11, b11+var*b11];
b22_set = [b22-var*b22, b22+var*b22];

params_set = [a11_set; a22_set; b11_set; b22_set];
%% 
jitter = 0.1;
true_params = [a11_set(1)+jitter*a11_set(1), a22_set(1)+jitter*a22_set(1)... 
                ,b11_set(1)+jitter*b11_set(1), b22_set(1)+jitter*b22_set(1)];

sys_low = makesysd_a(true_params(1), true_params(2), true_params(3), true_params(4),Ts);
sys_hih = makesysd_a(a11_set(2), a22_set(2), b11_set(2), b22_set(2), Ts);

% Built the System model
A = sys_low.A;
B = sys_low.B;
C = sys_low.C;

%% check feasibility 
value = 0;
Phi_low = sys_low.A - sys_low.B*Kopt;
Phi_hih = sys_hih.A - sys_hih.B*Kopt;

if value == 1
   cvx_begin sdp
   
   variable P(4, 4) symmetric
   
   minimize(trace(P))
   Phi_hih*P*Phi_hih' - P <= -eye(4)
   P>=0
   cvx_end 
   
   cvx_begin sdp
   
   variable P(4, 4) symmetric
   
   minimize(trace(P))
   Phi_low*P*Phi_low' - P <= -eye(4)
   P>=0
   cvx_end 
   
end

%%
% Time_out = 15;
mu = [a11, a22, b11, b22];

sigma = 1*diag(mu);
% time horizon window
p = 15;
maxF = 100;
main_bounds = [0.8, 0.2, 1]';
[H, f, Ac, Ax, b1, lb, ub, options] = MPC_vars(Aobv, Bobv, Cobv, Kopt, R, p, main_bounds, maxF);
%%

x = zeros(4, Time_out/Ts);
y = zeros(no_outputs, Time_out/Ts);

xhat = zeros(4, Time_out/Ts);
yhat = zeros(2, Time_out/Ts);

X = xhat(:, 1);
Ck = x(1, :);
%%

for k = 1 : Time_out/Ts - 1
    
    % Optimization part
    
    b = b1 + Ax*X;
    ck = quadprog(H, f, -Ac, -b, [], [], lb, ub, [], options);
    if isempty(ck)
        c = 0;
    else
        c = ck(1);
    end
    Ck(1, k) = c;
    
    
    v = 0.008*randn(no_states, 1);
    v(2, 1) = 0.01*v(2, 1);
    w = 0.01*randn(no_outputs, 1);
    
    
    
    %% System dynamics
    uk = -Kopt*xhat(:, k) + c;
    
    x(:, k+1) = A*x(:, k) + B*uk + v;
    y(:, k) = C*x(:, k) + w;
    
    
    %% Observer dynamics
    yhat(:, k) = Cobv*xhat(:, k);
    xhat(:, k+1) = Aobv*xhat(:, k) + Bobv*c + L*(y(:, k) - yhat(:, k));
    
    
    X = xhat(:, k+1); 
 
    
end

%%

figure
plot(y(1, :))
grid on
title('Cart Position MPC')

figure
plot(y(2, :))
grid on
title('Angle of Pendulum phi')

figure
stairs(Ck)

grid on
title('Reference Input MPC')




