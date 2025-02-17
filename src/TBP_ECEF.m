function [dXdt] = TBP_ECEF(t,X, mu)
%   By: Jorge Chavarín
% Equations of motion rewrited as a system of first order ordinary differential equation
%   
% INPUTS:
%   t       Time
%   X       State vector [position, velocity] or [r,v]
%   mu       Gravitational Parameter G(m1+m2)
%
% OUTPUT:
%   dXdt    Time derivate of state vector [vel., acc.]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transport Theorem problem equations of motion:
%
% Equation rewrited from I frame
% ṙ = dr / dt   = v              % v = vel.  --> Vi
% r̈ = dv / dt   = -(mu/r^3)*r    %r^3 -> scalar in ECI Frame   --> Ai
% W = Wf/i  = Angular velocity of the ECEF frame as seen from the ECI frame
% W = 7.2921x10-5 Z [rad/s] 
%
% At ECEF 
% i(ṙ) = i(Vf) = Vi - W x r           % Vf is vel. at ECEF seen from the ECI 
% ṙ = Vf = [FI] i(Vf)                 % Vf seen from ECEF
% [FI] --> Identity matrix as both frames are parallel at Xo means  ϴ = 0
% r̈ = Af = Ai - 2(W x Vf) - W x W x r
% r̈ = Af = Ai - 2(W x (Vi - W x r)) - W x W x r
% r̈ = -(mu/r^3)*r - 2(W x (Vi - W x r) - W x W x r
% X = [r, Vf]
% Ẋ = [Vf, Af]^T

% r scalar = [ (xPos)^2 + (yPos)^2 + (zPos)^2 } ^ (1/2) 
r = sqrt((X(1)^2) + (X(2)^2) + (X(3)^2) );

w = 7.2921e-5;               % Angular vel. of the ECEF frame. [rad/s]
W = [0;0;w];
posf = [X(1);X(2); X(3)];
velf = [X(4); X(5); X(6)]; 

% Constant value used at r̈    
Ai = -(mu/(r^3))*X(1:3);
FI = eye(3);
% Vector of zeros
dXdt = zeros(6,1);

% Vf = [EF] (Vi - W x r )
% Vf = Vi - W x r
%dXdt(1:3,1) =  [FI]* [X(4:6) - cross(W, X(1:3))];
dXdt(1:3, 1) = velf;
%  r̈ = Af = [EF] * [Ai - 2(W x (Vi - W x r)) - W x W x r]
%  r̈ = Af = [EF] * [Ai - 2(W x Vf) - W x W x r]
%  r̈ = Af = [EF] * [Ai - Coriolis - Centrifugal]

% Coriolis --> 2(W x (Vi - W x r))
coriolis = 2 * ( cross( W, velf ) );
% Centrifugal --> W x W x r
centrifugal = cross(W, cross(W,posf ));

dXdt(4:6, 1) = [FI]* [Ai - coriolis - centrifugal];

end