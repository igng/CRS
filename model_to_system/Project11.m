close all;
clearvars;
clc;

syms s;
S_model =   [-(s+1)^2, 0, 0, 2+s^2;
            2*s, -s, 0, 2*s - 1;
            s+1, -1, s-2, -3;
            2*s, 0, 0, -2*s - 1];
T = -S_model(1:3, 1:3);
U = S_model(1:3, 4);
V = S_model(4, 1:3);
W = S_model(4,4);

%%
detT = expand(simplify(det(T)));
r = rank(T);
n = length(coeffs(detT, 'all')) - 1;

%%
TF_model = simplify(V*(T^-1)*U + W);

%%
[L_1, R_1, T_1] = smithForm(T);
S_1_model = simplify([L_1, zeros(3, 1); zeros(1, 3), 1]*S_model*[R_1, zeros(3, 1); zeros(1, 3), 1]);

%%
S_2_model = [-eye(1), zeros(1, 4); zeros(4, 1), S_1_model];
T_2 = -S_2_model(1:n, 1:n);
U_2 = S_2_model(1:n, n+1);
V_2 = S_2_model(n+1, 1:n);
W_2 = S_2_model(n+1, n+1);

%%
A = [0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 2, 3, 0];
[L_3_bar, R_3_bar, T_2_temp] = smithForm(s*eye(n) - A);

%%
L_3 = L_3_bar^-1;
R_3 = R_3_bar^-1;
S_3_model = simplify([L_3, zeros(n, 1); zeros(1, n), 1]*S_2_model*[R_3, zeros(n, 1); zeros(1, n), 1]);
U_3 = S_3_model(1:n, n+1);
V_3 = S_3_model(n+1, 1:n);
W_3 = S_3_model(n+1, n+1);

%%
[N, D] = numden(V_3*((s*eye(n) - A)^-1));
N_3 = quorem(N, D);
S_4_model = simplify([eye(4), zeros(n, 1); N_3, 1]*S_3_model);

%%
[N, D] = numden(((s*eye(n) - A)^-1)*U_3);
N_4 = quorem(N, D);
S_5_model = simplify(S_4_model*[eye(4), N_4; zeros(1, n), 1]);

%%
S_system = S_5_model;
A = S_system(1:n, 1:n) + s*eye(n);
B = S_system(1:n, n+1);
C = S_system(n+1, 1:n);
D = S_system(n+1, n+1);

%%
syms xi x1 x2 x3 x4 u;
x = [x1; x2; x3; x4];
xi = simplify(R_1*[zeros(r, n-r), eye(r)]*R_3*(x + N_4*u));
TF_system = simplify(C*(s*eye(n) - A)^-1*B + D);
isequal(TF_model, TF_system);

%%
xi0 = simplify([eye(n-r), zeros(1, r)]*R_3*(x + N_4*u));
x_delta = A*x + B*u;
subs(xi0, s*x4, x_delta(4));