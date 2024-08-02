%% CLEAR

clear all 

clc

%% COFF

m = 0.2;

M = 0.5;

L = 0.3;

J = 0.006;

b = 0.1;

g = 9.81; 

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

s = tf('s'); 
I = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
F = C*inv((s*I-A))*B+D;

G_x = F(1);
G_phi = F(2);


%% H-Infinty

% Convert the transfer function G to state-space representation
%[A, B, C, D] = tf2ss(num, den);

% Define the sensitivity (WS) and complementary sensitivity (WT) weighting functions
%s = tf('s');

WS = 1 / (1 + 0.5 * s); % Sensitivity weighting function (low-pass filter)
WT = 1 / (1 + 0.5 * s); % Complementary sensitivity weighting function (unity gain)

% Define the unstructured uncertainty
Wu = [];

% Construct the generalized plant
P_x = augtf(G_x, WS, WT, Wu);
P_phi = augtf(G_phi, WS, WT, Wu);

% Synthesize the H-infinity controller
K_x = hinfsyn(P_x, 1, 1);
Kx_tf = tf(K_x);
K_phi = hinfsyn(P_phi, 1, 1);
Kphi_tf = tf(K_phi);

% Display the H-infinity controller
%disp('H-infinity Controller:');
%disp(K_tf);

% Plot the step response of the closed-loop system
%figure;
%step(feedback(G * K, 1));
%title('Closed-Loop Step Response');
