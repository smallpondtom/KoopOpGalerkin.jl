%% Koopman Operator Tutorial -- Duffing Oscillator No. 1 Example
% 
% S. Servadio, D. Arnas, and R. Linares, “A Koopman Operator Tutorial with 
% Othogonal Polynomials.” arXiv, Jul. 14, 2022. doi: 10.48550/arXiv.2111.07485.

clear; close all; clc;
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
addpath("src\")

%% Equation parameters
sm = 1; % Mass
sk = 1; % Spring constant
sa = 1; % Unit transformation constant
se = 0.001; % Small parameter

%% Other settings
nt = 100;

q0 = 1; % Initial q
p0 = 0.0; % Initial p
x0 = [q0; p0];

tf = 10; % Final time
tk = linspace(0,tf,nt);

%% Differential equation
% Elements in the array f(a,b,c):
% a: dimension where the derivative is performed
% b: term in the derivative
% c = 1: coefficient of the terms in xdot = f(x)
% c != 1: variable 'c' and its exponent in the polynomial term
Nx = 2;  % number of variables
Nt = 2;  % maximum number of terms in f(x)
fx = zeros(Nx,Nt,Nx+1);

% Linear system
fx(1,1,1) = 1; %x2/sm
fx(1,1,3) = 1;
fx(2,1,1) = -sk; %-x1*sk
fx(2,1,2) = 1;

% Perturbing term
%-sk*se*saˆ2*x1ˆ3
fx(2,2,1) = -sk*se*sa^2;
fx(2,2,2) = 3;

%% Koopman Solution
c = 3;
KSol = KoopOpSol(c,fx,x0,tk);
q = KSol(1,:);
p = KSol(2,:);

%% ODE45 results for comparison
param.sm = sm;
param.sk = sk;
param.se = se;
param.sa = sa;
[t,states] = ode45(@(t,x) duffing1(t,x,param), tk, x0);


%% Plotting results
fig1 = figure(1);
scatter(q,p,DisplayName="Koopman")
hold on; 
scatter(states(:,1),states(:,2),"+",DisplayName="analytical")
legend;
hold off

fig2 = figure(2);
plot(tk,q,lineWidth=2,DisplayName="$q_{k}$")
hold on
plot(tk,p,lineWidth=2,DisplayName="$p_{k}$")
plot(tk,states(:,1),'--',lineWidth=2,DisplayName="$q$")
plot(tk,states(:,2),'--',lineWidth=2,DisplayName="$p$")
legend;
hold off


%% Function
function dxdt = duffing1(t,x,p)
    dxdt = zeros(2,1);
    x1 = x(1);
    x2 = x(2);
    dxdt(1) = x2/p.sm;
    dxdt(2) = -p.sk*x1 - p.sk*p.se*p.sa^2*x1^3;
end