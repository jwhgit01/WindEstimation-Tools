function w = windReconstruction(V,alph,betaf,vi,Theta,omega,r_ADU)
%windReconstruction
%
% Copyright (c) 2024 Jeremy W. Hopwood. All rights reserved.
%
% This function reconstructs wind velocity using a vaned air data unit
% along with estimates/measurements of inertial velocity, attitude, and
% angular velocity. The reconstructed wind, w, satisfies
%
%            w = vi - R_IB*(v_ADU - cross(omega,r_ADU))
%
% where v_ADU = R_BW*e1*V. Here, vi is the NED inertial velocity, R_IB is
% the rotation matrix from the body frame to the inertial frame, R_BW is
% the rotation matrix from the wind frame to the body frame, e1 = [1;0;0],
% V is the airspeed, omega is the angular velocity of the body frame, and
% r_ADU is the position of the geometric center of the vanes in the body
% frame.
%
% Inputs:
%
%  
% Outputs:
%
%   w        The 3x1 vector of the reconstructed wind velocity.
%

% Basis vectors
e1 = [1;0;0];

% Rotation matrix from body to inertial
phi = Theta(1);
theta = Theta(2);
psi = Theta(3);
R1_phi = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
R2_theta = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
R3_psi = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];
R_IB = R3_psi*R2_theta*R1_phi;

% Rotation matrix from wind to body
beta = atan(tan(betaf)*cos(alph));
R2_alpha = [cos(-alph),0,sin(-alph);0,1,0;-sin(-alph),0,cos(-alph)];
R3_betaf = [cos(beta),-sin(beta),0;sin(beta),cos(beta),0;0,0,1];
R_BW = R2_alpha*R3_betaf;

% Reconstructed wind
w = vi - R_IB*(R_BW*e1*V - cross(omega,r_ADU));

end