clc, clear;
%Derivadas de las variables de estado_
beta_p = sym('beta_p');
p_p = sym('p_p');
r_p = sym('r_p');
phi_p = sym('phi_p');
psi_p = sym('psi_p');
x_p = [beta_p, p_p, r_p, phi_p, psi_p]
%Variables de estado_
beta = sym('beta');
p = sym('p');
r = sym('r');
phi = sym('phi');
psi = sym('psi');
x = [beta, p, r, phi, psi]
%Entradas del sistema
delta_a =sym('delta_a');
delta_r =sym('delta_r');
u = [delta_a, delta_r]
%valor de lambda para filtro IIR (Leaky Integrator)
lambda = 0.5;

%Matriz A
A = [-0.0558 0.08 -0.997 0.0415 0.0033;
    -3.05 -0.465 0.388 0 0;
    0.598 -0.318 -0.115 0 0;
    0 1 0 0 0;
    0 0 1 0 0]
%Matriz B
B = [0 0.00729;
    0.143 0.153;
    0.00775 -0.475;
    0 0;
    0 0]
%Matriz C
C = [0 0 0 1 0;
    0 0 0 0 1]
%Matriz D
D = [0 0;0 0]

Q = ctrb(A,B)
fprintf("El rango es %d", rank(Q))

Ob = obsv(A,C);
if rank(Ob)== 5
    fprintf("El sistema es observable")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aaug = [0 0 0 0 0;
        0 0 0 0 0]';
caug = [0 0; 0 0]';
baug = [0 0;
        0 0];
Aaug = [A aaug;
    -C caug]
Baug = [B; baug]
Caug = [C caug]
%Se hacen a las matrices diagonales porque ambas salidas están desacopladas
%%%%%%%%%%%%%%%%FUNCIONA BIEN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a_i = sqrt(0.1428);
a_i = 1;
rho_i = sqrt(0.0005);%0.001  0.0005

Qaug = [a_i^2/100^2 0 0 0 0 0 0;
    0 a_i^2/100^2 0 0 0 0 0;
    0 0 a_i^2/100^2 0 0 0 0;
    0 0 0 a_i^2/.5^2 0 0 0;  %.5
    0 0 0 0 a_i^2/.5^2 0 0;  %.5
    0 0 0 0 0 a_i^2/.1^2 0;  %.1
    0 0 0 0 0 0 a_i^2/.1^2]; %.1

Raug = [rho_i^2/1^2 0;
    0 rho_i^2/1^2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a_i = sqrt(1);
% rho_i = sqrt(1);%0.001  0.0005
% Qaug = [a_i^2/1^2 0 0 0 0 0 0;
%     0 a_i^2/1^2 0 0 0 0 0;
%     0 0 a_i^2/1^2 0 0 0 0;
%     0 0 0 a_i^2/.5^2 0 0 0;  %.5
%     0 0 0 0 a_i^2/.5^2 0 0;  %.5
%     0 0 0 0 0 a_i^2/.1^2 0;  %.1
%     0 0 0 0 0 0 a_i^2/.1^2]; %.1
% 
% Raug = [rho_i^2/2000 0;
%     0 rho_i^2/2000];


QQ = ctrb(Aaug, Baug)

rank(QQ)
%Obtención de las ganancias optimizadas
[Kaug, Saug, Paug] = lqr(Aaug,Baug,Qaug,Raug)

Aclaug = Aaug-Baug*Kaug

%Corregir Kraug
%Kraug = -inv(Caug*inv(Aclaug)*Baug)
Kp = Kaug(:,1:5);
Ki = [Kaug(:,6) Kaug(:,7)]

%LQI (LQR + Integrador)
sys = ss(A,B,C,D)
[Kai, Sai, Pai]=lqi(sys,Qaug,Raug)
KP = Kai(:,1:5)
KI = Kai(:,6:7)
%Covarianzas de ruido de proceso y de ruido de sensor
Qk_coef = 1.5;
% Qk = Qk_coef*eye(2);
Qk = [2.5e-2 0;0 1e-1]
%Qk = [1.5 0; 0 1.1];
%Qk = 1;
Rk = [0.0001 0; 0 0.0001]; %Funciona mejor que la covarianza igual a 1
%Rk = [1 0; 0 1];
[kalmf,L,~,Mx,Z] = kalman(sys,Qk,Rk)