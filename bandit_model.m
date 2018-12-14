%% initial estimated values
% set the seed
rng('default');

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
[no_outputs, ~] = size(C);

%% set of m, M, l
var = 0.4;

a11_set = [a11-var*a11, a11+var*a11];
a22_set = [a22-var*a22, a22+var*a22];

b11_set = [b11-var*b11, b11+var*b11];
b22_set = [b22-var*b22, b22+var*b22];

params_set = [a11_set; a22_set; b11_set; b22_set];
%% 

sys_low = makesysd_a(a11_set(1), a22_set(1), b11_set(1), b22_set(1), Ts);
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
Time_out = 1;
mu = [a11, a22, b11, b22];

sigma = 1*diag(mu);
% time horizon window
p = 20;
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

%% Test Models
no_models = 1000;
len_test = 5;
counter = 0;
    
a_mat = zeros(round(Time_out/Ts/len_test), 4);

y_buffer = zeros(no_outputs*no_models, len_test);
x_buffer = zeros(4*no_models, len_test);

a_init = [a11; a22; b11; b22]';
a_prev = a_init;

hyper_params = truncate_gauss2(mu, sigma, params_set, no_models);
%%
for k = 1 : Time_out/Ts - 1
    
 
    % Generate models
    [blkA_sparse, blkB_sparse, blkC_sparse, kalm_gain] = makeSparseBlkdiag(hyper_params, 4, Ts);
    
    
    %% Optimization part
    
    b = b1 + Ax*X;
    ck = quadprog(H, f, -Ac, -b, [], [], lb, ub, [], options);
    if isempty(ck)
        c = 0;
    else
        c = ck(1);
    end
    Ck(1, k) = c;
    
    
    v = 0.01*randn(no_states, 1);
    v(2, 1) = 0.1*v(2, 1);
    w = 0.001*randn(no_outputs, 1);
    
    
    
    %% System dynamics
    uk = -Kopt*xhat(:, k) + c;
    
    x(:, k+1) = A*x(:, k) + B*uk + v;
    y(:, k) = C*x(:, k) + w;
    
    
    %% Observer dynamics
    yhat(:, k) = Cobv*xhat(:, k);
    xhat(:, k+1) = Aobv*xhat(:, k) + Bobv*c + L*(y(:, k) - yhat(:, k));
    
    
    X = xhat(:, k+1); 
    
    %% run input through a test models
    % fill the buffer
    pntr = mod(k, len_test)+1;
    
    if pntr == 1
       y_buffer = zeros(no_outputs*no_models, len_test);
       x_buffer = zeros(4*no_models, len_test); 
       
       x_buffer(:, 1) = repmat(xhat(:, k), no_models, 1);
    end
    
    y_buffer(:, pntr) = blkC_sparse*x_buffer(:, pntr);
    kal_buff_gain = kalm_gain*(repmat(y(:, k), no_models, 1) - y_buffer(:, pntr));
    
    x_buffer(:, pntr+1) = blkA_sparse*x_buffer(:, pntr) + blkB_sparse*uk + kal_buff_gain;
    
    
    % find the optimal model 
    if pntr == len_test
       
        y_buff_reshape = [y_buffer(1:2:end, 1:end-1), y_buffer(2:2:end, 1:end-1)];
        y_data = [y(1, k-len_test+2:k), y(2, k-len_test+2:k)];
        
        [index, rms_vals] = rms_est(y_buff_reshape, y_data);
        
        a_new = hyper_params(index, :);
        
        % update observers and models
        new_sys = makesysd_a(a_new(1), a_new(2), a_new(3), a_new(4), Ts);
        A = new_sys.A;
        B = new_sys.B;
        C = new_sys.C;
        
        L = place(A', C', obvs_poles);
        L = L';
        
        [H, f, Ac, Ax, b1, lb, ub, options] = MPC_vars(A-B*Kopt, B, C, Kopt, R, p, main_bounds, maxF);
        
        kernel_func = 0.1*diag((a_prev - a_new).^2) + 0.0001*eye(length(a_prev));
        sigma = kernel_func;
        mu = a_new;
        
        hyper_params = truncate_gauss2(mu, sigma, params_set, no_models);
        % hyper_params = mvnrnd(mu, sigma, no_models);
        
        counter = counter + 1;
        a_mat(counter, :) = a_new;
        
    end
    
    
    
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


