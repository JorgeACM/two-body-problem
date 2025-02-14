function [E,i] = Kepler(e,M,tol)
% By: Jorge Chavarín
% Kepler's equation solution using Newton's method.
%
% INPUTS:
%   e       Eccentricity
%   M       Mean Anomally
%   tol     Error tolerance
%
% OUTPUTS:
%   E       Eccentric Anomaly
%   i       Number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M = E-e*sin(E)
% f(x,e,M) = E - e*sin(E) - M


% Initial values
E = M;
for i = 1:100                         % Evaluate function 100 times max.

    fE = E - e*sin(E) - M;            % evaluate function in En
    dfE = 1-e*cos(E);                 % evaluate first derivative in En

    % check for convergence using IF statement
    if(abs(fE) <= tol)
        break;
    else
        % Newton's method
        dE = -fE/dfE;                        % delta fE_n
        E = E + dE;                          % update value of fEn
    end
end



%{
% OPTION #2 While Loop
function [E,i] = Kepler_while(e,M,tol)

i = 1;                   % Variable to store number of iterations
E= M;
fE = E - e*sin(E)-M;
dfE = E-e*sin(E) -M; 
if abs(fE) <= tol
    i = 0;
else
% While Loop to execute the Newtons method.
    while abs(fE) > tol
        dfE = 1-e*cos(E);
        fE = E - e*sin(E) -M;
        dE = -(fE/dfE);       % dX = -( f(Xn-1) / f'(Xn-1) )
        E = E + dE;           % Xn = Xn-1 + dXn-1
        i = i +1;             % Calculate number of iterations
    end
end

%}