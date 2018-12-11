function sys_true = inv_p(Ts)

m = 1;
M = 4;
l = 1;
g = 9.81;

a11 = m*g/M;
a22 = (M+m)*g/(M*l);

b11 = 1/M;
b22 = 1/(M*l);

A = [0, 0, 1, 0;
     0, 0, 0, 1;
     0, a11, 0, 0;
     0, a22, 0, 0];
 
B = [0, 0, b11, b22]';

C = diag([1, 1, 0, 0]);

D = 0;

sys_true_c = ss(A, B, C, D);

sys_true = c2d(sys_true_c, Ts);



