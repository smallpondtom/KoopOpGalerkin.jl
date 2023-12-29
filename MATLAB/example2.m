%% Koopman Operator Tutorial -- Duffing Oscillator No. 2 Example
% 
% S. Servadio, D. Arnas, and R. Linares, “A Koopman Operator Tutorial with 
% Othogonal Polynomials.” arXiv, Jul. 14, 2022. doi: 10.48550/arXiv.2111.07485.

clear; close all; clc;
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
addpath("src\")

%% Equation parameters
% Duffing Oscillator 
%
% x1dot = x2
% x2dot = -alpha*x1 - beta*x1^3 - delta*x2 + gamma*cos(omega*t)
%
% where x1 and x2 are the states and alpha, beta, delta are some
% parameters. gamma*cos(omega*t) = gamma*u(t) is the control in which gamma
% and omega are parameters. HOWEVER, for this example we neglect the
% controls and treat it as an autonomous system.

alpha = 1;
beta = 2.0;
delta = 0.3;
% gamma = 2.5;
% omega = 2;

%% Other settings
nt = 1000;

x10 = 1; % Initial q
x20 = 0; % Initial p
x0 = [x10 x20];

tf = 20; % Final time
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
       1     0      0;
  -alpha -beta -delta
];
% Exponents of x1
fx(:,:,2) = [
    0 0 0;
    1 3 0
];
% Exponents of x2
fx(:,:,3) = [
    1 0 0;
    0 0 1;
];

%% Koopman Solution
c = 6;
KSol = KoopOpSol(c,fx,x0,tk);
x1 = KSol(1,:);
x2 = KSol(2,:);

%% ODE45 results for comparison
param.alpha = alpha;
param.beta = beta;
param.delta = delta;
[t,states] = ode45(@(t,x) duffing2noforce(t,x,param), tk, x0);


%% Plotting results
fig1 = figure(1);
scatter(x1,x2,DisplayName="Koopman")
hold on; 
scatter(states(:,1),states(:,2),"+",DisplayName="analytical")
legend;
hold off

fig2 = figure(2);
plot(tk,x1,lineWidth=2,DisplayName="$x_{1,k}$")
hold on
plot(tk,x2,lineWidth=2,DisplayName="$x_{2,k}$")
plot(tk,states(:,1),'--',lineWidth=2,DisplayName="$x_{1}$")
plot(tk,states(:,2),'--',lineWidth=2,DisplayName="$x_{2}$")
legend;
hold off


%% Function
function xdot = duffing2noforce(t,x,p)
    alpha = p.alpha; beta = p.beta; delta = p.delta;

    xdot = zeros(2,1);
    x1 = x(1);
    x2 = x(2);
    xdot(1) = x2;
    xdot(2) = -alpha*x1 - beta*x1^3 - delta*x2;
end