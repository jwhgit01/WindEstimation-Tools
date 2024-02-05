function [muw,Cw] = airDataUnitUncertainty(mu,Cx,r_ADU,N)
%airDataUnitUncertainty
%
% Copyright (c) 2024 Jeremy W. Hopwood. All rights reserved.
%
% This function computes the variance of reconstructed wind velocity using
% a vaned air data unit and estimates/measurements of inertial velocity,
% attitude, and angular velocity. The reconstructed wind, w, satisfies
%
%            w = vi - R_IB*(v_ADU - cross(omega,r_ADU))
%
% where v_ADU = R_BW*e1*V. Here, vi is the NED inertial velocity, R_IB is
% the rotation matrix from the body frame to the inertial frame, R_BW is
% the rotation matrix from the wind frame to the body frame, e1 = [1;0;0],
% V is the airspeed, omega is the angular velocity of the body frame, and
% r_ADU is the position of the geometric center of the vanes in the body
% frame. This computaion is performed under the assumption that the random
% errors in these measurements are mutually uncorrelated and Gasussian.
%
% Inputs:
%
%   mu      The 12x1 vector containing the expected value of
%           x = [V;alpha;betaf;vi;Theta;omega]
%
%   Cx      The 12x12 positive definite covariance matrix of
%           x = [V;alpha;betaf;vi;Theta;omega]
%
%   r_ADU   The position of the ADU in the body frame
%
%   N       The number of monte-carlo samples (optional). If not given or
%           given as an empty array, a default value of 10000 is used.
%  
% Outputs:
%
%   muw     The 3x1 mean value of the reconstructed wind
%
%   Cw      The 3x3 covariance matrix of the reconstructed wind
%

% Number of Monte-Carlo samples
if nargin < 4 || isempty(N)
    N = 10000;
end

% Latin hypercube sampling
X = lhsnorm(mu,Cx,N);

% Monte-Carlo
w = zeros(N,3);
for ii = 1:N
    V = X(ii,1);
    alph = X(ii,2);
    betaf = X(ii,3);
    vi = X(ii,4:6).';
    Theta = X(ii,7:9).';
    omega = X(ii,10:12).';
    w(ii,:) = windReconstruction(V,alph,betaf,vi,Theta,omega,r_ADU).';
end

% Mean
muw = mean(w).';

% Variance
Cw = cov(w); 

end % airDataUnitUncertainty