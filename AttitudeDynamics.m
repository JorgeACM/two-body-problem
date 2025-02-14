function [dXdt] = AttitudeDynamics(t,X, I)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

alpha = X(1);
beta = X(2);
gamma = X(3);
I1 = I(1,1); 
I2 = I(2,2);
I3 = I(3,3);

W = X(4:6);
matrix = [
    0,              sin(gamma),                cos(gamma); 
    0,              cos(beta)*cos(gamma),     -cos(beta)*sin(gamma); 
    cos(beta),     sin(beta)*sin(gamma),    sin(beta)*cos(gamma)];

B0 = (1/cos(beta)) * matrix  *W;

w_dot = [
    (I2-I3)/I1 * W(2) * W(3); 
    (I3-I1)/I2 * W(1) * W(3);
    (I1-I2)/I3 * W(1) * W(2)];

dXdt = [B0; w_dot];

end