clearvars;
close all;
clc;

set(0, 'DefaultFigureWindowStyle', 'docked');
%% Model description
syms s;
sI = s*eye(8);

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

C = [   0, 0, 0, 1, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 1;
        0, 0, 0, 0, 1, 0, 0, 0  ];

[n, ~] = size(A);
[~, p] = size(B);
[q, ~] = size(C);

%% 1) Verifiche preliminari
P = ctrb(A, B);
O = obsv(A, C);

if (rank(P) == n)
    disp("Il sistema è raggiungibile")
else
    disp("Il sistema NON è raggiungibile");
end

if (rank(O) == n)
    disp("Il sistema è osservabile")
else
    disp("Il sistema NON è osservabile");
end

if (p == q)
    disp("L'impianto è quadrato");
else
    disp("L'impianto NON è quadrato");
end

%% 2) LQR (pag 1029)
Cq = eye(q, n);
Q = Cq.'*Cq;
Oq = obsv(A, Cq);

if (rank(Oq) == n)
    disp("La coppia (A, Cq) è rilevabile");
else
    disp("La coppia (A, Cq) NON è rilevabile");
end

Pq = Cq*((sI-A)\B);
McSq = mcs(Pq);
ZdTq = double(zdt(McSq));

fprintf("\nGli zeri di trasmissione della matrice Cq*(sI-A)^-1*B");

if (~ sum(ZdTq > 0))
    fprintf(" NON ");
end
fprintf("hanno parte reale maggiore o uguale a zero\n");

disp(ZdTq);

R = 1e2*eye(3);
K_inf = lqr(A, B, Q, R);

SL = K_inf*((sI-A)\B);

U = (eye(p) + SL)\(SL);
[numU, denU] = numden(U);
Uh = tf(zeros(3));
for i = 1:3
    for j = 1:3
        nn = double(coeffs(numU(i,j), 'All'));
        dd = double(coeffs(denU(i,j), 'All'));
        Uh(i,j) = tf(nn, dd);
    end
end

close(figure(1));
figure(1);
sigma(Uh, 'm'); hold on;
set(findall(gcf,'type','line'),'linewidth', 3)

%% LTR (pag 1032)
Pp = C*((sI-A)\B);
McS = mcs(Pp);
ZdTp = double(zdt(McS));

fprintf("\nGli zeri di trasmissione della matrice C*(sI-A)^-1*B");
if (~ sum(ZdTp == 0))
    fprintf(" NON ");
end

fprintf("sono sull'asse immaginario\n");

disp(ZdTp);

for sigma_v = 100.^(1:3)
    W = eye(p);
    V = sigma_v^2*(B*B.');
    L = lqr(A', C', V, W)';

    SLo = -K_inf*((sI - A + L*C + B*K_inf)\(L*C*((sI-A)\B)));

    Uo = (eye(p) + SLo)\(SLo);
    [numUo, denUo] = numden(Uo);
    Uho = tf(zeros(3));
    for i = 1:3
        for j = 1:3
            nno = double(coeffs(numUo(i,j), 'All'));
            ddo = double(coeffs(denUo(i,j), 'All'));
            Uho(i,j) = tf(nno, ddo);
        end
    end
    
    sigma(Uho);
end
grid on; hold off;
legend('0', '1e2', '1e4', '1e6');