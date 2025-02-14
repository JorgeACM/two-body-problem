function [dXdt] = TBP_ECI(t, X, u)
%   By: Jorge Chavarín
% Equations of motion rewrited as a system of first order ordinary differential equation
%   
% INPUTS:
%   t       Time
%   X       State vector [position, velocity] or [r,v]
%   u       Gravitational Parameter G(m1+m2)
%
% OUTPUT:
%   dXdt    Time derivate of state vector [vel., acc.]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two-body problem equations of motion:
% r̈ = d^2r / dt = -(mu/r^3)*r.   % r^3 -> scalar.
%
% Equation rewrited.
% ṙ = dr / dt   = v              % v = vel.
% r̈ = dv / dt   = -(mu/r^3)*r    %r^3 -> scalar
% 
% X = [r, v]
% Ẋ = [v, -(mu/r^3)*r]^T  = [rdot, vdot]

% r c = [ (xPos)^2 + (yPos)^2 + (zPos)^2 } ^ (1/2) 
r = sqrt((X(1)^2) + (X(2)^2) + (X(3)^2) );

% Constant value used at r̈    
cte = -(u/(r^3));

% Vector of zeros
dXdt = zeros(6,1);

dXdt(1) = X(4);
dXdt(2) = X(5);
dXdt(3) = X(6);

dXdt(4) = cte * X(1);
dXdt(5) = cte * X(2);
dXdt(6) = cte * X(3);

end


