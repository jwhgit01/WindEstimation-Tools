close all
clear
clc
addpath ../src

% Mean
mu_V = 18;
mu_alpha = 5*pi/180;
mu_betaf = 0;
mu_vi = 20*randn(3,1);
mu_eul = [0;0;0];
mu_omega = zeros(3,1);
mu = [mu_V;mu_alpha;mu_betaf;mu_vi;mu_eul;mu_omega];

% Covariance
s_V = 0.1;
s_alpha = pi/180;
s_betaf = pi/180;
s_vi = 0.01*ones(3,1);
s_eul = [0.001;0.001;0.01];
s_omega = 0.0001*ones(3,1);
Cx = diag([s_V;s_alpha;s_betaf;s_vi;s_eul;s_omega].^2);

% ADU position
r_ADU = [1;0;0];

% Direct reconstruction
w = windReconstruction(mu_V,mu_alpha,mu_betaf,mu_vi,mu_eul,mu_omega,r_ADU)

% Mean and covariance of reconstructed wind
[muw,Cw] = airDataUnitUncertainty(mu,Cx,r_ADU)
