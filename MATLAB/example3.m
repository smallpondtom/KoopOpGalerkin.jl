%% Koopman Operator Tutorial -- Van Der Pol Oscillator Example
% 
% S. Servadio, D. Arnas, and R. Linares, “A Koopman Operator Tutorial with 
% Othogonal Polynomials.” arXiv, Jul. 14, 2022. doi: 10.48550/arXiv.2111.07485.

clear; close all; clc;
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
addpath("src\")

%% Equation parameters
% Van Der Pol Oscillator 
%
% x1dot = x2
% x2dot = mu(1 - x1^2)x2 + x1
%
% where x1 and x2 are the states and mu is a parameter. 

mu = 5;

%% Other settings
nt = 1000;

x10 = 1; % Initial x1
x20 = 0.1; % Initial x2
x0 = [x10 x20];

tf = 100; % Final time
tk = linspace(0,tf,nt);

%% Differential equation
% Elements in the array f(a,b,c):
% a: dimension where the derivative is performed
% b: term in the derivative
% c = 1: coefficient of the terms in xdot = f(x)
% c != 1: variable 'c' and its exponent in the polynomial term
Nx = 2;  % number of variables
Nt = 3;  % maximum number of terms in f(x)
fx = zeros(Nx,Nt,Nx+1);

% Coefficients
fx(:,:,1) = [
   1  0   0;
  -1 mu -mu
];
% Exponents of x1
fx(:,:,2) = [
    0 0 0;
    1 0 2
];
% Exponents of x2
fx(:,:,3) = [
    1 0 0;
    0 1 1;
];

%% Koopman Solution
c = 3;
KSol = KoopOpSol(c,fx,x0,tk);
x1 = KSol(1,:);
x2 = KSol(2,:);

%% ODE45 results for comparison
[t,states] = ode45(@(t,x) vanderpol(t,x,mu), tk, x0);


%% Plotting results
fig1 = figure(1);
scatter(states(:,1),states(:,2),"+",DisplayName="analytical")
legend;

fig2 = figure(2);
scatter(x1,x2,DisplayName="Koopman")
hold on; 
scatter(states(:,1),states(:,2),"+",DisplayName="analytical")
legend;
hold off

tcut = 1:nt;
fig3 = figure(3);
plot(tk(tcut),x1(tcut),lineWidth=2,DisplayName="$x_{1,k}$")
hold on
plot(tk(tcut),x2(tcut),lineWidth=2,DisplayName="$x_{2,k}$")
plot(tk(tcut),states(tcut,1),'--',lineWidth=2,DisplayName="$x_{1}$")
plot(tk(tcut),states(tcut,2),'--',lineWidth=2,DisplayName="$x_{2}$")
legend;
hold off


%% Function
function xdot = vanderpol(t,x,mu)
    xdot = zeros(2,1);
    x1 = x(1);
    x2 = x(2);
    xdot(1) = x2;
    xdot(2) = -x1 + mu*x2 - mu*x1^2*x2;
end