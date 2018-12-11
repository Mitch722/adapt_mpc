function sys_var_d = makesysd(m, M, l, Ts)

g = 9.81;

a11 = m*g/M;
a22 = (M+m)*g/(M*l);

b11 = 1/M;
b22 = 1/(M*l);

A = [1, 0, Ts, 0;
     0, 1, 0, Ts;
     0, a11*Ts, 0, 0;
     0, a22*Ts, 0, 0];
 
B = [0, 0, b11*Ts, b22*Ts]';

C = [1, 0, 0, 0;
     0, 1, 0, 0];

sys_var_d.A = A;
sys_var_d.B = B;
sys_var_d.C = C;

