% By: Jorge Chavarin
clear all; close all; clc;
format longG;

%% DATA
u = 398600.4418;        % [km3/s2]  Gravitational Parameter G(m1+m2)  
Re = 6378.137;          % [km]      Earth Ratio
w = 7.2921e-5;          % [rad/s]   Angular vel. Earth rotation
Ixx = 2500;             % [kg m^2]
Iyy = 5000;             % [kg m^2]
Izz = 6500;             % [kg m^2]
% m2km = 1/ 1000;
% Ixx = Ixx*(m2km^2);
% Iyy = Iyy*(m2km^2);
% Izz = Izz*(m2km^2);
I = [                   % Spacecraft inertia matrix
    Ixx, 0,     0;          % in principal axes frame.
    0,   Iyy,   0; 
    0,   0,     Izz];

%% Part 1: Two-Body Problem %%
% Keplerian orbit elements:
a = 7151.16 ;           % [km]      Semi-major axis
e = 0.0008;             % []        Eccentricity of the Orbit

% Inclination of the Orbit
i_deg = 98.39;          % [deg]
i = deg2rad(i_deg);     % [rad]

% Right Ascension of the Ascending Node (RAAN)
RAAN_deg = 10.0;        % [deg]
RAAN = deg2rad(10.0);   % [rad]

% Argument of Periapsis
argPer_deg = 233;       % [deg]
argPer = deg2rad(233);  % [rad]

Mo_deg = 127;           % [deg]     Initial Mean Anomaly
Mo = deg2rad(127);      % [rad]

tol = 1e-10;            % []        Tolerance factor

%% Week 2
% Solving Keplers equation using Newton's Method and Find Initial True
% Anomaly
fprintf('<strong>Mean Anomaly (Mo)</strong>\n Mo = <strong>%.4f [deg]</strong>.\n\n', Mo_deg);
[E,ii] = Kepler(e,Mo,tol);      % [rad] E = Eccentric Anomaly
Edeg = rad2deg(E);              % [deg] Eccentric Anomaly
fprintf('<strong>Eccentric Anomaly (E)</strong> after %d iterations:\n', ii);
fprintf('E = <strong>%.4f [rad]</strong> = <strong>%.4f [deg]</strong>.\n\n', E, Edeg);

% Calculating True Anomaly

%%%% True Anomaly Equation %%%%
theta = 2*atan2(sqrt(1+e) * tan(E/2), sqrt(1-e));
thetadeg = rad2deg(theta);                          % [deg] True Anomaly
fprintf('<strong>True Anomaly (θ)</strong>:\nθ = <strong>%.4f [rad]</strong> = <strong>%.4f [deg]</strong>.\n', theta, thetadeg);
fprintf('\nThe satellite is near <strong>apoapsis</strong>.\n\n');

%% Week 3
%%% Orbit Propagation via Analytical %%
%%% Solution of Two-Body Problem %%%%

% Orbit Period
n = sqrt(u/a^3);                        % [rad/s]   Mean  Motion (rad/s)
P = 2*pi/n;                             % [s]       Period of the Orbit(s)
fprintf('Orbit Period: %.4f (s)\n', P);


% Create time vector of 1000 equally spaced points t∈[0,P]
% assuming t0 = 0 s
t0 = 0;
t = linspace(t0, P+t0, 1000);

% Calculate Mean Anomally
MA = Mo + n*(t-t0);                     % [rad] Mean Anomaly

% Empty vectors of size P, to store Anomalies
TA = t; EA = t;                         % [rad] vectors
EAdeg = t; TAdeg = t; MAdeg = t;        % [deg] vectors

for j = 1:length(t)  
    [EA(j),ii] = Kepler(e,MA(j),tol);    % [rad] EA = Excentricity Anomaly
    TA(j) = 2*atan2(sqrt(1+e) * tan(EA(j)/2), sqrt(1-e)); % [rad] True Anomaly
    
    % Classic Orbits Elements
    % coe = [a; e; i; Ω; ω; θ]
    coe(:,j) = [a, e, i, RAAN, argPer, TA(j)];
    
    % Calculate ECI position and velocity
    X(:,j) = COE2RV(coe(:,j), u);      % [[km], [km/s]] [pos, vel]
end

for jj = 1:6
    XX(jj, :) = [X(jj, 1), X(jj, end)]; 
end

% X = [x; y; z; ẋ; ẏ; ż]
fprintf('X      Initial values       Final Values \n')
fprintf('x:  %15f%18f\n', XX(1,1), XX(1,2));
fprintf('y:  %15f%18f\n', XX(2,1), XX(2,2));
fprintf('z:  %15f%18f\n', XX(3,1), XX(3,2));
fprintf('ẋ: %15f%18f\n', XX(4,1), XX(4,2));
fprintf('ẏ: %15f%18f\n', XX(5,1), XX(5,2));
fprintf('ż: %15f%18f\n', XX(6,1), XX(6,2));

% Plot Anomalies
% Radians
figure('Name','Anomalies','units','normalized','outerposition',[0 0 1 1]); hold on; grid on;
plot(t, wrapTo2Pi(MA), 'LineWidth', 6); plot(t, wrapTo2Pi(EA), 'LineWidth', 4); plot(t, wrapTo2Pi(TA), 'LineWidth', 2);
legend('Mean','Eccentric', 'True');
title('Anomalies');
xlabel('time [s]')
ylabel('Anomaly [rad]')

hold off;

%{
% Degrees
% Convert Anomalies form radians to degrees
EAdeg = rad2deg(EA);          % [deg] Eccentric Anomaly
TAdeg = rad2deg(TA);          % [deg] True Anomaly
MAdeg = rad2deg(MA);          % [deg] Mean Anomaly
figure('Name','Anomalies','units','normalized','outerposition',[0 0 1 1]); hold on; grid on;
plot(t, MAdeg, 'LineWidth', 4); plot(t, EAdeg, 'LineWidth', 2); plot(t, TAdeg, 'LineWidth', 2);
legend('Mean','Eccentric', 'True');
title('Anomalies');
xlabel('time [s]')
ylabel('Anomaly [deg]')
hold off;
%}

% Degrees
%{
MA_pi_deg = rad2deg(MA_pi);
EA_pi_deg = rad2deg(EA_pi);
TA_pi_deg = rad2deg(TA_pi);
figure('Name','Anomalies','units','normalized','outerposition',[0 0 1 1]); hold on; grid on;
plot(t, MA_pi_deg, 'LineWidth', 4); plot(t, EA_pi_deg, 'LineWidth', 2); plot(t, TA_pi_deg, 'LineWidth', 2);
legend('Mean','Eccentric', 'True');
title('Anomalies');
xlabel('time [s]')
ylabel('Anomaly [deg]')
hold off;
%}

%3D Plot of the ECI spacecraft trajectory
% Create figure
figure('Name', 'ECI Frame','units','normalized','outerposition',[0 0 1 1])
[Xe,Ye,Ze] = sphere(50);

% Plot Earth
img = imrotate(imread('BlueMarble_square.png'), 0);
surf(Re*Xe, Re*Ye, Re*Ze, 'CData', img, 'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha','.8'); hold on; axis equal;
% surf(Re*Xe, Re*Ye, Re*Ze, 'EdgeColor', '#A2142F', 'FaceColor', 'c', 'FaceAlpha','1'); hold on; axis equal;
title('ECI Spacecraft Trajectory.');
xlabel('I, ECI (km)');
ylabel('J, ECI (km)');
zlabel('K, ECI (km)');

% Plot ECI axes
quiver3(0, 0, 0, 1, 0, 0, 1e4, 'k', 'LineWidth', 2);    % I-axis
quiver3(0, 0, 0, 0, 1, 0, 1e4, 'k', 'LineWidth', 2);    % J-axis
quiver3(0, 0, 0, 0, 0, 1, 1e4, 'k', 'LineWidth', 2);    % K-axis
text(1e4, 0, 0, 'I', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'k')
text(0, 1e4, 0, 'J', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'k')
text(0, 0, 1e4, 'K', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'k')


% Plot Trajectory
plot3(X(1,:), X(2,:), X(3,:), 'k', 'LineWidth', 2);
plot3(X(1,1), X(2,1), X(3,1), 'ok', 'MarkerFaceColor', 'b');
plot3(X(1,end), X(2,end), X(3,end), 'ok', 'MarkerFaceColor', 'r');


%%%{
% Calculate eccentricity vector, specific angular momentum vector, and
% complete the triad
r = X(1:3, 1);
v = X(4:6, 1);
h = cross(r, v);
ee = cross(v, h)/u - r/norm(r);

ie = ee/norm(ee);
ih = h/norm(h);
ip = cross(ih, ie)/norm(cross(ih, ie));

quiver3(0, 0, 0, ie(1), ie(2), ie(3), 1e4, 'r', 'LineWidth', 2);
quiver3(0, 0, 0, ip(1), ip(2), ip(3), 1e4, 'r', 'LineWidth', 2);
quiver3(0, 0, 0, ih(1), ih(2), ih(3), 1e4, 'r', 'LineWidth', 2);
text(1e4*ie(1), 1e4*ie(2), 1e4*ie(3), 'i_e', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'r')
text(1e4*ip(1), 1e4*ip(2), 1e4*ip(3), 'i_p', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'r')
text(1e4*ih(1), 1e4*ih(2), 1e4*ih(3), 'i_h', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'r')

hold off
%%%}
%% Week 4 %%%%%
ode_opt = odeset('RelTol', 3e-14, 'AbsTol', 1e-16);

[tout, Xout] = ode113(@TBP_ECI, t, X(:, end), ode_opt, u);

% Function to extract ode from main script
%[tout, Xout] = TBP_ECI_ode(t, X(:, end), u);   % Ẋ = [ṙ; ̇v] =[v; -(mu/r^3)*r]

% Position and Velocity error Vector over time t 
dX = Xout.' - X;

% Calculate Positionand Velocity Scalar error.
eV = zeros(1,length(t));
eX = zeros(1,length(t));
for tt = 1:length(t)
    eX(tt) = sqrt( dX(1,tt)^2 + dX(2,tt)^2 + dX(3,tt)^2);
    eV(tt) = sqrt( dX(4,tt)^2 + dX(5,tt)^2 + dX(6,tt)^2);
end

%Plot Position and Velocity Errors. Two-Body Problem vs Numerical Integration.
figure('Name', 'Pos & Vel Errors');
title("Two-Body Problem vs Numerical Integration.");
yyaxis left
plot(t,eX)
ylabel("Position Error [km].")
yyaxis right
plot(t,eV)
ylabel("Velocity Error [km/s].")
xlabel('t [s]');

% Week 4.3 Specific Energy 
% ζ = 1/2(v^2)- mu/r = cte
specific_energy = zeros(1, length(t)) -u/(2*a);                % [km/s]
SE = zeros(1,length(t));                    % Specific Energy ζ zeros vector
err = zeros(1,length(t));                   % Error zeros vector
for tt = 1:length(t)
    r = sqrt( Xout(tt, 1)^2 + Xout(tt, 2)^2 + Xout(tt, 3)^2);   % [km]
    v = sqrt( Xout(tt, 4)^2 + Xout(tt, 5)^2 + Xout(tt, 6)^2 );  % [km/s]
    SE(tt) =  ( (v^2) /2) - ( u/r );        % Specific Energy ζ [km^2/s^2] using actual pos and vel
    err(tt) = SE(tt) - -u/(2*a);     % Calculating error
    %specific_energy(tt) = -u/(2*a);           % [km/s]
end

% Plot Specific Energy
figure('Name', 'Specific Energy graphics'); hold on; grid on;
plot(t, specific_energy, 'LineWidth', 2); 

plot(t, SE, 'LineWidth', 2); 
% legend('−µ/2a', 'ζ');
lg =legend('$-\frac{\mu}{2a}$','$\zeta = \frac{1}{2}v^2 - \frac{\mu}{r}$');
set(lg, 'Interpreter', 'latex');
title('Specific Energy graphics');
xlabel('time [s]');
ylabel('$\frac{km^2}{s^2}$', 'Interpreter', 'latex');
% $x^2+e^{\pi i}$ 
ytickformat('%.4f');
hold off;

%% Week 5
period_earth = 86164;
% period_diff = period_earth/P;
t2 = linspace(t0, 10*P, 1000);      % Time vector of 10 Orbit Periods
FI = eye(3);                        % I3 Identity matrix.   
W = [0;0;w];                        % Vector Angular vel. of Earth rotation. [rad/s]
Xin = Xout(end, :).';               % [pos; vel]

% Calculationg Velocity from ECEF
% ṙ = Vf = [FI] i(Vf)                 % Vf seen from ECEF
Xin(4:6, 1) = [FI] * [Xin(4:6) - cross(W, Xin(1:3))];   % [vx; vy; vz] [km/s]

ode_opt = odeset('RelTol', 3e-14, 'AbsTol', 1e-16);
[tout, Xecef] = ode113(@TBP_ECEF, t2, Xin, ode_opt, u);

%[tout, Xecef] = TBP_ECEF_ode(t2, Xin, u);
Xecef = Xecef.';

% Create figure
figure('Name', 'ECEF spacecraft trajectory', 'units','normalized','outerposition',[0 0 1 1]);

% Plot Earth
%img = imrotate(imread('BlueMarble_square.png'), 0);
surf(Re*Xe, Re*Ye, Re*Ze, 'CData', img, 'EdgeColor', 'none', 'FaceColor', 'texturemap', 'FaceAlpha','.9'); hold on; axis equal;
% surf(Re*Xe, Re*Ye, Re*Ze, 'EdgeColor', 'none', 'FaceColor', 'c'); hold on; axis equal;
xlabel('X, ECEF [km]');
ylabel('Y, ECEF [km]');
zlabel('Z, ECEF [km]');
title('ECEF spacecraft trajectory');
% Plot ECEF axes
quiver3(0, 0, 0, 1, 0, 0, 1e4, 'k', 'LineWidth', 2);    % I-axis
quiver3(0, 0, 0, 0, 1, 0, 1e4, 'k', 'LineWidth', 2);    % J-axis
quiver3(0, 0, 0, 0, 0, 1, 1e4, 'k', 'LineWidth', 2);    % K-axis
text(1e4, 0, 0, 'X', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'k')
text(0, 1e4, 0, 'Y', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'k')
text(0, 0, 1e4, 'Z', 'FontSize', 20, 'Interpreter', 'tex', 'Color', 'k')


% Plot Trajectory
plot3(Xecef(1,:), Xecef(2,:), Xecef(3,:), 'k', 'LineWidth', 2);
plot3(Xecef(1,1), Xecef(2,1), Xecef(3,1), 'ok', 'MarkerFaceColor', 'b');
plot3(Xecef(1,end), Xecef(2,end), Xecef(3,end), 'ok', 'MarkerFaceColor', 'r');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Week 6
% Spacecraft ECI position and velocity coordinates 1 hour after periapsis.
r_w6 = [6768.27;  870.9;2153.59];     % [km] Position vector
v_w6 = [-2.0519; -1.4150; 7.0323];  % [km s^-1] Velocity vector

% RSW Orbital Frame
h_w6 = cross(r_w6, v_w6);           % h = r x ṙ

er = r_w6 / norm(r_w6);
eh = h_w6 / norm(h_w6);
e0 = cross(eh,er);

IO = [er, e0, eh];
fprintf('\n\nThe [IO] matrix is:\n<strong>[IO]</strong>=\n');
for j = 1:3
    fprintf('%.4f  %.4f  %.4f\n', IO(j,1), IO(j,2), IO(j,3));
end
OI = transpose(IO);
fprintf('\n\nThe [OI] matrix is:\n<strong>[OI]</strong>=\n');
for j = 1:3
    fprintf('%.4f  %.4f  %.4f\n', OI(j,1), OI(j,2), OI(j,3));
end

% θ =[α; β; γ]
alpha = 30;                 % [deg] First axis - roll 
beta = 20;                  % [deg] Second axis - pitch
gamma = 10;                 % [deg] First axis rot- roll

% Angles in radians
% betaRad = deg2rad(beta);
% alphaRad = deg2rad(alpha);
% gammRad = deg2rad(gamma);

% Direction Cosine Matrix [BO]
%BO2 = Euler_matrix("roll", gamma)*Euler_matrix("pitch", beta) * Euler_matrix("roll", alpha);
BO = [  cosd(beta), sind(beta)*sind(alpha), -sind(beta)*cosd(alpha);
        sind(gamma)*sind(beta), -sind(gamma)*cosd(beta)*sind(alpha)+cosd(gamma)*cosd(alpha), sind(gamma)*cosd(beta)*cosd(alpha)+cosd(gamma)*sind(alpha);
        cosd(gamma)*sind(beta), -cosd(gamma)*cosd(beta)*sind(alpha)-sind(gamma)*cosd(alpha), cosd(gamma)*cosd(beta)*cosd(alpha)-sind(gamma)*sind(alpha)
];
fprintf(['\nDirection cosine matrix <strong>[BO}</strong> from the orbital frame <strong>O</strong> ' ...
    'to the principal axes frame of the satellite <strong>B</strong>:\n\n<strong>[BO]</strong>=\n']);
for j = 1:3
    fprintf('      %.4f  %.4f  %.4f\n', BO(j,1), BO(j,2), BO(j,3));
end

% B = [BO][OI]I
% B = [BI]I
% [BI] = [BOI][OI]
BI = BO*OI;
fprintf(['\nDirection cosine matrix <strong>[BI]</strong> from the ECI Frame ' ...
    '<strong>I</strong> to the principal axes frame of the spacecraft <strong>B</strong>:\n']);
for j = 1:3
    fprintf('%.4f  %.4f  %.4f\n', BI(j,1), BI(j,2), BI(j,3));
end
% Comparing [BO][OI] with BI given from the assignment. 
% Result should be 0 o very near cero due to floating point in matlab
BI_assignment = [0.7908, 0.3705, 0.4872; -0.0474, -0.7565, 0.6523; 0.6102, -0.5389, -0.5807];
% BI_compare = BI - BI_assignment;


%% Week 7

% Calculating values of Euler's Principal Axis and Angle

% [BJ] = cos(Φ) I3 - sin(Φ)[ẽ] + ( 1-cos(Φ) )êê^T
% [BJ] = [b11, b12, b13; b21, b22, b23; b31, 32, b33]
% Φ = cos^-1( [BJ]^T -1 / 2 )
% e1 = (b23-b32) / (2sinΦ)
% e2 = (b31-b13) / (2sinΦ)
% e3 = (b12-b21) / (2sinΦ)
% Φ = phi
% E = ẽ
% EE = êê^T

phi_rad = acos( (trace(BI)-1) / 2);     % [rad] Principal Angle
phi = rad2deg(phi_rad);                 % [deg] Principal Angle
den = 2*sin(phi_rad);
e1 = (BI(2,3)-BI(3,2))/den;             % e1 Axis
e2 = (BI(3,1)-BI(1,3))/den;             % e2 Axis
e3 = (BI(1,2)-BI(2,1))/den;             % e3 Axis
rot_axis = [e1; e2; e3];                % Priincipal Axis vector

fprintf("Euler's Principal Angle & Axis <strong>(ê,Φ)</strong>:\n");
fprintf(" <strong>Φ</strong> = <strong>%.4f rad</strong>\n", phi_rad);
fprintf("<strong>e1</strong> = <strong>%.4f</strong>\n" + ...
    "<strong>e2</strong> = <strong>%.4f</strong>\n" + ...
    "<strong>e3</strong> = <strong>%.4f</strong>\n", e1, e2, e3);
% Quaternions

B0 = cosd(phi/2);
Bcte = sind(phi/2);
B1 = e1*Bcte;
B2 = e2*Bcte;
B3 = e3*Bcte;
% B^2 B/J = B0^2 = B1^2 = B2^2 = B3^2;
BB = B0^2 + B1^2 + B2^2 + B3^2; 
B = cosd(phi/2)^2 + sind(phi/2)^2*(e1^2+e2^2+e3^2);

fprintf("\nEuler's Parameters (Quaternions):\n");
fprintf("<strong>β</strong> = %.4f + %.4f<strong>i</strong> +" + ...
    " %.4f<strong>j</strong> + %.4f<strong>k</strong>\n", B0, B1, B2, B3);
 
fprintf("\nSince β0 = %.4f is <strong>nonnegative</strong> it means β corresponds to the shortest rotation.\n", B0);

% Since B0 is nonnegative it means B corresponds to the shortest rotation.
% as phi < 180 (phi =140)

% Calculate 3-2-1 Euler Angles sequence
% Wb/j = alpha_dot K + beta_dot a + gamma_dot F
% A = {a, b, c}
% D = {d, e, f}
% [BD] = R1(roll)
% [BA] = [BD][DA] = R1(roll)R2(pitch)
% [BJ] = [BD][DA][AJ] = R1(roll)R2(pitch)R3(yaw)
% [BJ] = [BD][DA][AJ] = R1(alpha)R2(beta)R3(gamma)
%^B(Wb/j) = [BJ](0;0;gamma) + [BA](0;beta;0) + [BD](alpha;0;0)
%^B(Wb/j) = R3(gamma)R2(beta)R1(alpha)(alpha;0;0) + R3(gamma)R2(beta)(0;beta;0) + R3(gamma)(0;0;gamma)
% r_B = BI*r_ECI;

% Eulers angles in radians
yaw_rad = atan2(BI_assignment(1,2),BI_assignment(1,1));
pitch_rad = asin(-BI_assignment(1,3));
roll_rad = atan2(BI_assignment(2,3),BI_assignment(3,3));
% Eulers angles in degrees
yaw = rad2deg(yaw_rad);
pitch = rad2deg(pitch_rad);
roll = rad2deg(roll_rad);
eulers_angles = [yaw_rad; pitch_rad; roll_rad];

% Week 7.3 pt1: Reporting values of 3-2-1 Eule Angles sequence values
fprintf('Euler Angles = [ψ, θ, φ]\n' );
fprintf('<strong>ψ</strong> = <strong>%.4f rad</strong>    ', yaw_rad);
fprintf('<strong>θ</strong> = <strong>%.4f rad</strong>    ', pitch_rad);
fprintf('<strong>φ</strong> = <strong>%.4f rad</strong>\n\n', roll_rad);

fprintf('<strong>ψ</strong> = <strong>%.4f deg</strong>   ', yaw);
fprintf('<strong>θ</strong> = <strong>%.4f deg</strong>   ', pitch);
fprintf('<strong>φ</strong> = <strong>%.4f deg</strong>\n\n', roll);
%{
% 3-2-1 Euler Angles sequence values matrix [BJ] compare with original [BI]
BJ = [
    cos(pitch_rad)*cos(yaw_rad), cos(pitch_rad)*sin(yaw_rad), -sin(pitch_rad);
    sin(roll_rad)*sin(pitch_rad)*cos(yaw_rad)-cos(roll_rad)*sin(yaw_rad), sin(roll_rad)*sin(pitch_rad)*sin(yaw_rad)+cos(roll_rad)*cos(yaw_rad), sin(roll_rad)*cos(pitch_rad);
    cos(roll_rad)*sin(pitch_rad)*cos(yaw_rad)+sin(roll_rad)*sin(yaw_rad), cos(roll_rad)*sin(pitch_rad)*sin(yaw_rad)-sin(roll_rad)*cos(yaw_rad), cos(roll_rad)*cos(pitch_rad)];

fprintf('Direction cosine matrix <strong>[BJ]</strong> of a 3-2-1 Euler Angles sequence(yaw, pitch, roll):\n\n');
fprintf('<strong>[BJ]</strong> =\n');
for j = 1:3
    fprintf('%3.4f  %3.4f  %2.4f\n', BJ(j,1), BJ(j,2), BJ(j,3));
end
%}
% Week 7.3 pt2 Reporting values of B(θ)
B_theta = [
    0,              sin(roll_rad),                  cos(roll_rad);
    0,              cos(pitch_rad)*cos(roll_rad),   -cos(pitch_rad)*sin(roll_rad);
    cos(pitch_rad), sin(pitch_rad)*sin(roll_rad),   sin(pitch_rad)*cos(roll_rad)];
B_cte = 1/cos(pitch_rad); 


fprintf('Kinematic Differential Equation <strong>θ˙= B(θ)ω</strong> of a 3-2-1 Euler Angles sequence(yaw, pitch, roll):\n\n');
fprintf('<strong>B(θ)</strong> =\n');
for j = 1:3
    if j == 2
        fprintf('(%.4f) * ', B_cte);
    else
        fprintf('           ');
    end
    fprintf('[%3.4f  %3.4f  %2.4f] (w%1d)\n', B_theta(j,1), B_theta(j,2), B_theta(j,3),j);
end


%% Week 8

w1 = -3.092e-4;                         % [rad/s]
w2 = 6.6161e-4;                         % [rad/s]
w3 = 7.4606e-4;                         % [rad/s]

wB = [w1; w2; w3];                      % [rad/s]

X_B = [eulers_angles; wB];              %[rad, rad/s]

t_attitude = linspace(0, 3600, 361);    %t0=0[s], tend=3600[s]

ode_opt = odeset('RelTol', 3e-12, 'AbsTol', 1e-14);
% Different solutions using ode45() and ode113()
% [tout, X_Att_Dyn] = ode45(@AttitudeDynamics, t_attitude, X_B, ode_opt, I);
[tout, X_Att_Dyn] = ode113(@AttitudeDynamics, t_attitude, X_B, ode_opt, I);
X_Att_Dyn=X_Att_Dyn.';                  % [rad, rad/s]
X_Att_Dyn(1:3,:) = wrapToPi(X_Att_Dyn(1:3,:));
X_Att_Dyn_deg = X_Att_Dyn*180/pi;       % [deg, deg/s]


% Plot euler angle and angular velocities in one 2x1 window.
figure('Name', 'Attitude State of the spacecraft.', 'units','normalized','outerposition',[0 0 1 1]);
title('Euler Properties');
tiledlayout(2,1);

nexttile; hold on; grid on;
title('Angles')
xlabel('time [s]');
ylabel('θ [rad]');
plot(tout, X_Att_Dyn(1,:), 'LineWidth', 2); plot(tout, X_Att_Dyn(2,:), 'LineWidth', 2); plot(tout, X_Att_Dyn(3,:), 'LineWidth', 2);
legend('Yaw','Pitch', 'Roll');
hold off;

nexttile; hold on; grid on;
title('Angular velocities');
xlabel('time [s]');
ylabel('w [rad/s]');
plot(tout, X_Att_Dyn(4,:), 'LineWidth', 2); plot(tout, X_Att_Dyn(5,:), 'LineWidth', 2); plot(tout, X_Att_Dyn(6,:), 'LineWidth', 2);
legend('w1','w2', 'w3');
hold off;

% Plot euler angle and angular velocities in one 2x3 window.
figure('name', 'Attitude State of the Spacecraft.', 'units','normalized','outerposition',[0 0 1 1]);
tiledlayout(2,3);

%ψ, θ, φ
nexttile
plot(tout, X_Att_Dyn(1,:), 'LineWidth', 2, 'Color', "#0072BD");
legend('Yaw');
ylabel('ψ [rad]');
xlabel('time [s]');
grid on;

nexttile
plot(tout, X_Att_Dyn(2,:), 'LineWidth', 2, 'Color', "#D95319");
legend('Pitch');
ylabel('θ [rad]');
xlabel('time [s]');
grid on;

nexttile
plot(tout, X_Att_Dyn(3,:), 'LineWidth', 2, 'Color', "#EDB120");
legend('Roll');
ylabel('φ [rad]');
xlabel('time [s]');
grid on;

nexttile
plot(tout, X_Att_Dyn(4,:), 'LineWidth', 2, 'Color', "#0072BD");
legend('w1');
ylabel('w1 [rad/s]');
xlabel('time [s]');
grid on;

nexttile
plot(tout, X_Att_Dyn(5,:), 'LineWidth', 2, 'Color', "#D95319");
legend('w2');
ylabel('w2 [rad/s]');
xlabel('time [s]');
grid on;

nexttile
plot(tout, X_Att_Dyn(6,:), 'LineWidth', 2, 'Color', "#EDB120");
legend('w3');
ylabel('w3 [rad/s]');
xlabel('time [s]');
grid on;

fprintf(['\nBecause of the 3-2-1 sequence naturality of the rotation, ' ...
    'Pitch angle (θ) almost reaches singularity (+-π).\n' ...
    '<strong>θ ≈ -85.8 deg | θ ≈ 1.49 rad</strong> at <strong>t = 1100 s</strong>\n']);

% Angular Momentum
H = [Ixx*X_Att_Dyn(4,:); Iyy*X_Att_Dyn(5,:); Izz*X_Att_Dyn(6,:)];
H_radius = sqrt(H(1,:).^2 + H(2,:).^2 + H(3,:).^2);

% Rotational Kinetic Energy

for n=1:361
    T(n) = (1/2)*(X_Att_Dyn(4,n).'*H(1,n) + X_Att_Dyn(5,n).'*H(2,n) + X_Att_Dyn(6,n).'*H(3,n));
end


% Plot euler angle and angular velocities in one 2x3 window.
figure('name', 'Angular momentum & Rotational Kinetic Energy.', 'units','normalized','outerposition',[0 0 1 1]);
tiledlayout(2,1);

nexttile
plot(tout, H_radius, 'LineWidth', 2);
legend('Angular Momentum');
ylabel('H [kg m^2 / s]');
xlabel('time [s]');
grid on;

nexttile
plot(tout, T, 'LineWidth', 2,'Color', "#D95319");
legend('Rotational Kinetic Energy');
ylabel('T [kg m^2 / s^2]');
xlabel('time [s]');
grid on;
