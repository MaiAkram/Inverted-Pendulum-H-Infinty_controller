% MODEL PARAMETERS
M = 0.5; %0.5;
m = 0.3; %0.3;
b = 0.1;
J = 0.006; %0.004;
g = 9.8;
L = 0.3; %0.3;

q = (M+m)*(J+m*L^2)-(m*L)^2;
s = tf('s');

% TRANSFER FUNCTIONS
%G_x = (((J+m*L^2)/q)*s^2 - (m*g*L/q)) / (s^4 + (b*(J + m*L^2))*(s^3)/q - ((M + m)*m*g*L)*(s^2)/q - b*m*g*L*s/q);
G_x = ((J+m*L^2)/q) / (s^2 - ((M + m)*m*g*L)/q);
%G_x = ((J+m*L^2)*(s^2)-m*g*L) / (((J*(M+m)+M*m*(L^2)))*s^2 - ((M + m)*m*g*L));
G_phi = (m*L*s/q) / (s^3 + (b*(J + m*L^2))*(s^2)/q - ((M + m)*m*g*L)*s/q - b*m*g*L/q);
%G_phi = (m*L) / (((m^2)*(L^2)-(M+m)*(J+(m*(L^2))))*(s^2) + (M+m)*m*g*L); % MK TF

% PID CONTROLLER GAINS
Kp_p = 100; %-120; %75;
Ki_p = 1; %-355; %120; %150
Kd_p = 20; %-10; %9; %8.64

Kp = 100;
Ki = 1; %120; %150
Kd = 20; %8.64

Kp_c = 80;   %100
Ki_c = 145;  %50
Kd_c = 10.8; %33

% STATE SPACE AND TRANSFER FUNCTIONS
% x' = Ax + Bu
% y = Cx + Du
% x' = [x' x'' phi' phi'']
% y = [x phi]
A = [0 1 0 0;
     0 -(J+m*L^2)*b/(J*(M+m)+M*m*L^2) (m^2)*g*(L^2)/(J*(M+m)+M*m*L^2) 0;
     0 0 0 1;
     0 -(m*L*b)/(J*(M+m)+M*m*L^2) (m*g*L)*(M+m)/(J*(M+m)+M*m*L^2) 0];
B = [0;
     (J+m*L^2)/(J*(M+m)+M*m*L^2);
     0;
     (m*L)/(J*(M+m)+M*m*L^2)];
C = [1 0 0 0;
     0 0 1 0];
D = [0; 0];
% I = [1 0 0 0;
%      0 1 0 0;
%      0 0 1 0;
%      0 0 0 1];
% F = C*inv((s*I-A))*B+D;
% G_x = F(1);

% CONTROLLABILITY CHECK
r_c = rank(ctrb(A,B));

% OBSERVABILITY CHECK
r_o = rank(obsv(A,C));

% Define the sensitivity (WS) and complementary sensitivity (WT) weighting functions
s = tf('s');

WS = 1 / (0.53 * s + 1); % Sensitivity weighting function (low-pass filter)
WT = 1 / (0.53 * s + 1); % Complementary sensitivity weighting function (unity gain)
W1 = (0.0999*(s^2)+0.99*s+0.25)/(0.7006*(s^2)+1.044*s+2.36);
W2 = (2.805*s+0.1536)/(0.1104*s+4.194);
W3 = (0.103*(s^2)+0.4988*s+0.5596)/(0.0115*(s^2)+0.5047*s+1.5);
WR = (2*s/3+10)/(s+0.001);
WP = WR;%0.32*(s+3)/(s+0.1);

% Define the unstructured uncertainty
Wu = [];

% Construct the generalized plant
P_x = augtf(G_x, WS, WT, Wu);
P_phi = augtf(G_phi, WS, WT, Wu);

% Synthesize the H-infinity controller
K_x = hinfsyn(P_x, 1, 1);
Kx_tf = tf(K_x);
%Kx_tf = ((s^4)+(s^3)+53*(s^2)-0.314*s-105)/((s^4) + 5.406*(s^3) + 70.39*(s^2) + 61.66*s);
K_phi = hinfsyn(P_phi, 1, 1);
Kphi_tf = tf(K_phi);
