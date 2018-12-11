function sys_var_d = makesysd_a(a11, a22, b11, b22, Ts)

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

