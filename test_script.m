%% initial estimated values
m = 1;
M = 4;
l = 1;
g = 9.81;

Ts = 0.001;
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
obvs_poles = multiple * real(L_eig);
obvs_poles = [0.0001, 0.001, 0.002, 0.0002]';

L = place(sysd.A', sysd.C', obvs_poles);
L = L';

% Built the Observer model

A = [sysd.A - sysd.B*Kopt, sysd.B*Kopt;
        zeros(size(sysd.A)),    sysd.A - L*sysd.C];
    
B = [ sysd.B;
         zeros(size( sysd.B ))];
         
C = [sysd.C, zeros(size(sysd.C))];

D = [0; 0];
%%
[~, no_states] = size(A);
[no_outputs, ~] = size(C);

%% set of m, M, l
var = 0.4;

a11_set = [a11-var*a11, a11+var*a11];
a22_set = [a22-var*a22, a22+var*a22];

b11_set = [b11-var*b11, b11+var*b11];
b22_set = [b22-var*b22, b22+var*b22];

%% 

sys_low = makesysd_a(a11_set(1), a22_set(1), b11_set(1), b22_set(1), Ts);
sys_hih = makesysd_a(a11_set(2), a22_set(2), b11_set(2), b22_set(2), Ts);

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

sigma = 0.05*diag(mu);
% time horizon window
p = 30;
maxF = 100;
main_bounds = [0.8, 0.2, 1]';
[H, f, Ac, Ax, b1, lb, ub, options] = MPC_vars(A, B, C, Kopt, R, p, main_bounds, maxF);
%%

x = zeros(no_states, Time_out/Ts);
y = zeros(no_outputs, Time_out/Ts);

X = x(:, 1);
Ck = x(1, :);

for k = 1 : Time_out/Ts - 1
    
    % kernel_func = 1./(a_prev - a_new).^2
    % sigma = diag(kernel_func)
    R = mvnrnd(mu,sigma,2);
    
    b = b1 + Ax*X;
    ck = quadprog(H, f, -Ac, -b, [], [], lb, ub, [], options);
    if isempty(ck)
        c = 0;
    else
        c = ck(1);
    end
    Ck(1, k) = c;
    
    
    w = 0.01*randn(no_states, 1);
    v = 0.001*randn(no_outputs, 1);
    
    x(:, k+1) = A*x(:, k) + B*c + w;
    y(:, k) = C*x(:, k) + v;
    
    X = x(:, k+1);

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


