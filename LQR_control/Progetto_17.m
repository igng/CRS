clearvars;
close all;
clc;

set(0, 'DefaultFigureWindowStyle', 'docked');
%% Model description
syms s;

A = [   -0.07, -0.017, 16.62, -18.4, 0.001, -1.0, 0.02, -0.07;
        0.04, -0.65, 0.14, -1.39, -0.04, 0.07, -0.33, -0.03;
        0.01, 0.007, -2.72, -2.22, 0.0002, 0.15, -0.001, -0.04;
        0, 0, 1, 0, 0, 0, 0, 0;
        -0.007, -0.006, -0.97, 0.005, -0.14, -6.91, 22.3, 3.76;
        -0.0006, 0.003, -0.81, 0.001, -0.014, -4.56, -6.26, 0.63;
        0, 0, 0, 0, 0, 1, 0, 0;
        0.007, 0.015, -0.55, 0.0001, 0.014, -1.03, -0.92, -3.68    ];
    
B = [   -2.2, 0.54, 0, 0.0001;
        -0.01, -12.1, -314.45, 0;
        0.36, -0.003, -0.001, 0.008;
        0, 0, 0, 0;
        -0.034, -0.17, 1.81, -1.0;
        0.093, -0.098, 1.09, -0.25;
        0, 0, 0, 0;
        0.25, 0.04, 0.04, 0.73  ];
    
B = B(:, 1:3);
% B = B(:, [1 2 4]);
% B = B(:, [1 3 4]);
% B = B(:, 2:4);
    
C = [   0, 0, 0, 1, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 1;
        0, 0, 0, 0, 1, 0, 0, 0  ];

[n, ~] = size(A);
[~, p] = size(B);
[q, ~] = size(C);

D = zeros(q, p);

omega = 0.1;    

gamma1 = 0;         % reference
gamma2 = 1i*omega;  % noise
opt = stepDataOptions('StepAmplitude', -1);
S2b = ss(A, B, C, D);
%% Condizioni di esistenza della soluzione al problema | pag 422-423
% C_hat = C
% Specifica (0)
if ~D
    ben_conn = 1;
    fprintf(1, 'La matrice di trasferimento ingresso-uscita è nulla, quindi il sistema è sicuramente ben connesso\n');
end

% Specifica (a1)
g1_stab = rank([A - gamma1*eye(n) B]);
g2_stab = rank([A - gamma2*eye(n) B]);
cb_stab = (g1_stab == n && g2_stab == n);
if (cb_stab)
    fprintf(1, 'La terna (A, B, C) è Cb-stabilizzabile\n');
else
    fprintf(1, 'La terna (A, B, C) non è Cb-stabilizzabile\n');
end
g1_rilev = rank([A - gamma1*eye(n); C]);
g2_rilev = rank([A - gamma2*eye(n); C]);
cb_rilev = g1_rilev == n && g2_rilev == n;
if (cb_rilev)
    fprintf(1, 'La terna (A, B, C) è Cb-rilevabile\n');
else
    fprintf(1, 'La terna (A, B, C) non è Cb-rilevabile\n');
end
a1_cond = cb_stab && cb_rilev;
if (a1_cond)
    fprintf(1, 'La condizione a1) è soddisfatta\n');
else
    fprintf(1, 'La condizione a1) non è soddisfatta quindi non esiste soluzione al problema proposto\n');
end

% Specifica (b)
rank_g1 = (rank([A - gamma1*eye(n) B; C D]) == (n + q));
rank_g2 = (rank([A - gamma2*eye(n) B; C D]) == (n + q));
b_cond  = rank_g1 && rank_g2;
if (b_cond)
    fprintf(1, 'La condizione b) è soddisfatta\n');
else
    fprintf(1, 'La condizione b) non è soddisfatta quindi non esiste soluzione al problema proposto\n');
end
% Conclusione
if (ben_conn && a1_cond && b_cond)
    fprintf(1, 'Esiste una controllore dinamico K tale da soddisare le specifiche (0), (a1), (b) e (c1)\n');
else
    fprintf(1, 'Non esiste una controllore dinamico K tale da soddisare le specifiche (0), (a1), (b) e (c1)\n');
end
%% Modello interno | pag 430-431
syms x
phi = expand((x - gamma1)*(x - gamma2)*(x - gamma2'));
cof = coeffs(phi, 'all');
mu = length(cof) - 1;
% A_phi = [zeros(mu-1, 1) eye(mu-1); -cof(end:-1:2)];
A_phi = [0 1 0; 0 0 1; 0 -1/100 0];
B_phi = [zeros(mu-1, 1); 1];
A_M1 = blkdiag(A_phi, A_phi, A_phi);
B_M1 = blkdiag(B_phi, B_phi, B_phi);
C_M1 = [1 0 0 0 1 0 0 0 1;
        0 1 0 1 0 1 0 1 0;
        0 0 1 0 0 0 1 0 0];
% rank(obsv(A_M1, C_M1))
D_M1 = zeros(q, p);
S2a = ss(A_M1, B_M1, C_M1, D_M1);
%% Interconnessione [S2a + S2b]
S2 = series(S2a, S2b);
[A2, B2, C2, D2] = ssdata(S2);
[n2, ~] = size(A2);
[~, p2] = size(B2);
[q2, ~] = size(C2);
%% LQR (1/3)
% Retroazione dallo stato 
Q2 = eye(n2);
R2 = eye(p2);
K = lqr(A2, B2, Q2, R2);
eK = eig(A2 - B2*K);
% disp('Autovalori di A2 - B2*K');
% disp(eK);
%% LQR (2/3)
% Guadagno del filtro di Kalman
W2 = eye(n2);
V2 = eye(p2);
L = lqr(A2',C2', W2, V2)';
eL = eig(A2 - L*C2);
% disp('Autovalori di A2 - L*C2');
% disp(eL);
%% LQR (3/3)
A1 = A2 - B2*K - L*C2 + L*D2*K;
B1 = L;
C1 = -K;
D1 = zeros(size(K, 1), size(L, 2));
S1 = ss(A1, B1, C1, D1);
CL = feedback(series(S1, S2), eye(q), 1); % retroazione positiva
% [Acl,~,~,~] = ssdata(CL);
% eCl = eig(Acl);
% disp(eCl);
figure(1)
step(CL, opt);
grid on; hold on;
%% Velocizzazione dei modi tramite diversi pesi
% % Pesiamo di meno il controllo ==> possiamo avere un controllo più aggressivo
% % LQR (1/3)
% % Retroazione dallo stato 
% Q2 = eye(n2);
% R2 = 1e-3*eye(p2);
% K = lqr(A2, B2, Q2, R2);
% eK = eig(A2 - B2*K);
% % LQR (2/3)
% % Guadagno del filtro di Kalman
% W2 = eye(n2);
% V2 = 1e-3*eye(p2);
% L = lqr(A2',C2', W2, V2)';
% % eL = eig(A2 - L*C2);
% % disp('Autovalori di A2 - L*C2');
% % disp(eL);
% % LQR (3/3)
% A1 = A2 - B2*K - L*C2 + L*D2*K;
% B1 = L;
% C1 = -K;
% D1 = zeros(size(K, 1), size(L, 2));
% S1 = ss(A1, B1, C1, D1);
% CL = feedback(series(S1, S2), eye(q), 1); % retroazione positiva
% % [Acl,~,~,~] = ssdata(CL);
% % eCl = eig(Acl);
% % disp(eCl);
% step(CL, opt);
% % legend('Controllore standard', 'Controllore aggressivo');
% % grid on; hold on;
%% Velocizzazione dei modi tramite la traslazione
Q2 = eye(n2);
R2 = eye(p2);
W2 = eye(n2);
V2 = eye(p2);
figure(10)
for alfa = 0.3:0.3:1.5
    K_alfa = lqr(A2 + alfa*eye(size(A2)), B2, Q2, R2);
    L_alfa = lqr(A2' + alfa*eye(size(A2)), C2', W2, V2)';
    eK_alfa = eig(A2 - B2*K);
    eL_alfa = eig(A2 - L*C2);
    A1_alfa = A2 - B2*K_alfa - L_alfa*C2 + L_alfa*D2*K_alfa;
    B1_alfa = L_alfa;
    C1_alfa = -K_alfa;
    D1_alfa = zeros(size(K_alfa,1),size(L_alfa,2));

    S1_alfa = ss(A1_alfa, B1_alfa, C1_alfa, D1_alfa);
    CL_alfa = feedback(series(S1_alfa, S2), eye(q), 1); % retroazione positiva
%     [Acl_alfa,~,~,~] = ssdata(CL_alfa);
%     eCl_alfa = eig(Acl_alfa);
%     disp(eCl_alfa);
%     figure;
    step(CL_alfa, opt); grid on; hold on
end
%% Perturbazione additiva sulla coppia (A, B)
epsilon = 0.05;
delta_A = zeros(size(A));
delta_B = zeros(size(B));
% for i = 1:n
%     for j = 1:n
%         delta_A(i,j) = epsilon*rand*A(i,j);
%     end
%     for l = 1:p
%         delta_B(i,l) = epsilon*rand*B(i,l);
%     end
% end

for i = 1:n
    for j = 1:n
        delta_A(i,j) = (-1)^(rand > 0.5)*epsilon*A(i,j);
    end
    for l = 1:p
        delta_B(i,l) = (-1)^(rand > 0.5)*epsilon*B(i,l);
    end
end

A_pert = A + delta_A;
B_pert = B + delta_B;
S2b_pert = ss(A_pert, B_pert, C, D);
S2_pert  = series(S2a, S2b_pert);
CL = feedback(series(S1, S2_pert), eye(q), 1); % retroazione positiva
figure(100);
step(CL, opt);
grid on;