%% Koopman Operator Tutorial -- Lorenz Attractor Example
% 
% S. Servadio, D. Arnas, and R. Linares, “A Koopman Operator Tutorial with 
% Othogonal Polynomials.” arXiv, Jul. 14, 2022. doi: 10.48550/arXiv.2111.07485.

clear; close all; clc;
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
addpath("src\")

%% Equation parameters
% Lorenz system
%
% x1dot = -sigma*x1 + sigma*x2
% x2dot = rho*x1 - x2 - x1x3
% x3dot = x1x2 - beta*x3
%
% where x1, x2, and x3 are the states and mu is a parameter. 

% Attractor
sigma = 10;
rho = 28;
beta = 8/3;

% sigma = 0.1;
% rho = 0.1;
% beta = 0.2;

%% Other settings
nt = 10000;

x10 = 1; % Initial x1
x20 = 2; % Initial x2
x30 = 3; % Initial x3
x0 = [x10 x20 x30];

tf = 100; % Final time
tk = linspace(0,tf,nt);

%% Differential equation
% Elements in the array f(a,b,c):
% a: dimension where the derivative is performed
% b: term in the derivative
% c = 1: coefficient of the terms in xdot = f(x)
% c != 1: variable 'c' and its exponent in the polynomial term
Nx = 3;  % number of variables
Nt = 3;  % maximum number of terms in f(x)
fx = zeros(Nx,Nt,Nx+1);

% Coefficients
fx(:,:,1) = [
    -sigma  sigma   0;
       rho     -1  -1;
         1  -beta   0
];
% Exponents of x1
fx(:,:,2) = [
    1 0 0;
    1 0 1;
    1 0 0
];
% Exponents of x2
fx(:,:,3) = [
    0 1 0;
    0 1 0;
    1 0 0
];
% Exponents of x3
fx(:,:,4) = [
    0 0 0;
    0 0 1;
    0 1 0;
];

%% Koopman Solution
c = 8;
KSol = KoopOpSol(c,fx,x0,tk);
x1 = KSol(1,:);
x2 = KSol(2,:);
x3 = KSol(3,:);

%% ODE45 results for comparison
params.sigma = sigma;
params.rho = rho;
params.beta = beta;
[t,states] = ode45(@(t,x) lorenz(t,x,params), tk, x0);


%% Plotting results
fig1 = figure(1);
plot3(states(:,1),states(:,2),states(:,3),DisplayName="analytical")
legend;

fig2 = figure(2);
plot3(x1,x2,x3,DisplayName="Koopman")
hold on; 
plot3(states(:,1),states(:,2),states(:,3),DisplayName="analytical")
legend;
hold off

tcut = 1:nt;
fig3 = figure(3);
tiledlayout(3,1);
    nexttile;
    plot(tk(tcut),states(tcut,1),'--',lineWidth=2,DisplayName="$x_{1}$")
    hold on
    plot(tk(tcut),x1(tcut),lineWidth=2,DisplayName="$x_{1,k}$")
    hold off
    legend
    
    nexttile;
    plot(tk(tcut),states(tcut,2),'--',lineWidth=2,DisplayName="$x_{2}$")
    hold on
    plot(tk(tcut),x2(tcut),lineWidth=2,DisplayName="$x_{2,k}$")
    hold off
    legend

    nexttile;
    plot(tk(tcut),states(tcut,3),'--',lineWidth=2,DisplayName="$x_{3}$")
    hold on
    plot(tk(tcut),x3(tcut),lineWidth=2,DisplayName="$x_{3,k}$")
    hold off
    legend;



%% Function
function xdot = lorenz(t,x,p)
    sigma = p.sigma; rho = p.rho; beta = p.beta;
    xdot = zeros(3,1);
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    xdot(1) = -sigma*x1 + sigma*x2;
    xdot(2) = rho*x1 - x2 - x1*x3;
    xdot(3) = x1*x2 - beta*x3;
end