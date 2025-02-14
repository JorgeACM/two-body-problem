function [X] = COE2RV(coe,u)
% By: Jorge Chavarín
% Convert the state of a spacecraft from classic orbit elements 
% to ECI position and velocity components
% INPUTS:
%   coe  Classic Orbit Elements Vector 
%               coe         [a, e, i, RAAN, argPer, TA];
%               a           Semi-major axis
%               e           Eccentricity of the Orbit
%               i           Inclination of the Orbit
%               RAAN        Right Ascension of the Ascending Node (RAAN)
%               argPer      Argument of Periapsis
%               TA          True Anomaly
%
%   u    Gravitational Parameter G(m1+m2)
% OUTPUT:
%   X   ECI position and velocity coordinates

a = coe(1);
e = coe(2);
i = coe(3);
RAAN = coe(4);
argPer = coe(5);
TA = coe(6);


% Rotation Matrix JP
JP = [
    cos(RAAN)*cos(argPer) - sin(RAAN)*cos(i)*sin(argPer),   -cos(RAAN)*sin(argPer) - sin(RAAN)*cos(i)*cos(argPer),  sin(RAAN)*sin(i);
    sin(RAAN)*cos(argPer) + cos(RAAN)*cos(i)*sin(argPer),   -sin(RAAN)*sin(argPer) + cos(RAAN)*cos(i)*cos(argPer),  -cos(RAAN)*sin(i);
    sin(i)*sin(argPer),                                     sin(i)*cos(argPer),                                     cos(i);
    ];
% Map position and inertial velocity vectors from P'frame to ECI.
r = ( a*(1-e^2) ) / (1 + e*cos(TA));
h = sqrt(u*a*(1-e^2));

pos = JP*[r*cos(TA);r*sin(TA);0];                   % Position Coordinates
vel = JP*[-(u/h)*sin(TA); (u/h)*(cos(TA)+e); 0];    % Velocity Coorcinates

% ECI position and velocity coordinates
X = [pos; vel];
end