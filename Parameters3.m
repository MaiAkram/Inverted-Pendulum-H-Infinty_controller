% MODEL PARAMETERS
M = 2;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
L = 0.5;
q = (M+m)*(I+m*L^2)-(m*L)^2;
s = tf('s');

% TRANSFER FUNCTIONS
G_x = (((I+m*L^2)/q)*s^2 - (m*g*L/q))/(s^4 + (b*(I + m*L^2))*s^3/q - ((M + m)*m*g*L)*s^2/q - b*m*g*L*s/q);
G_phi = (m*L*s/q)/(s^3 + (b*(I + m*L^2))*s^2/q - ((M + m)*m*g*L)*s/q - b*m*g*L/q);

% STATE SPACE AND TRANSFER FUNCTIONS
% x' = Ax + Bu
% y = Cx + Du
% x' = [x' x'' phi' phi'']
% y = [x phi]
A = [0 1 0 0;
     0 -(I+m*L^2)*b/(I*(M+m)+M*m*L^2) (m^2)*g*(L^2)/(I*(M+m)+M*m*L^2) 0;
     0 0 0 1;
     0 -(m*L*b)/(I*(M+m)+M*m*L^2) (m*g*L)*(M+m)/(I*(M+m)+M*m*L^2) 0];
B = [0;
     (I+m*L^2)/(I*(M+m)+M*m*L^2);
     0;
     (m*L)/(I*(M+m)+M*m*L^2)];
C = [1 0 0 0;
     0 0 1 0];
D = [0; 0];

% PID CONTROLLER GAINS
Kp_p = -146.230455344679;
Ki_p = -407.797870671116;
Kd_p = -12.8757934020323;

Kp_c = -1;
Ki_c = 0;
Kd_c = -3;

% Define the sensitivity (WS) and complementary sensitivity (WT) weighting functions
s = tf('s');

%WS = 1 / (1 + 0.5 * s); % Sensitivity weighting function (low-pass filter)
WS = 17*(s+9) / (s+3000);
%WT = 1 / (1 + 0.5 * s); % Complementary sensitivity weighting function (unity gain)
WT = 17*(s+120) / (s+5000);
W1 = (0.0999*(s^2)+0.99*s+0.25)/(0.7006*(s^2)+1.044*s+2.36);
W2 = (2.805*s+0.1536)/(0.1104*s+4.194);
W3 = (0.103*(s^2)+0.4988*s+0.5596)/(0.0115*(s^2)+0.5047*s+1.5);

% Define the unstructured uncertainty
Wu = [];

% Construct the generalized plant
%P_x = augtf(G_x, WS, WT, Wu);
P_x = augtf(G_x, W1, W2, W3);
P_phi = augtf(G_phi, WS, WT, Wu);

% Synthesize the H-infinity controller
K_x = hinfsyn(P_x, 1, 1);
Kx_tf = tf(K_x);
%Kx_tf = ((s^4)+(s^3)+53*(s^2)-0.314*s-105)/((s^4) + 5.406*(s^3) + 70.39*(s^2) + 61.66*s);
K_phi = hinfsyn(P_phi, 1, 1);
Kphi_tf = tf(K_phi);
