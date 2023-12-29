%% Koopman Operator Tutorial
% 
% S. Servadio, D. Arnas, and R. Linares, “A Koopman Operator Tutorial with 
% Othogonal Polynomials.” arXiv, Jul. 14, 2022. doi: 10.48550/arXiv.2111.07485.

clear; close all; clc;
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

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

% alpha = 1;
% beta = 5;
% delta = 0.02;
% gamma = 2.5;
% omega = 2;


%% Order of the basis functions
c = 3;

%% Equation parameters
sm = 1; % Mass
sk = 1; % Spring constant
sa = 1; % Unit transformation constant
se = 0.001; % Small parameter

%% Number of points in the figure
nt = 100;

q0 = 1; % Initial q
p0 = 0.0; % Initial p
tf = 10; % Final time
tk = linspace(0,tf,nt);



%% Number of points in the figure
% nt = 100;
% 
% x10 = 1; % Initial q
% x20 = 0; % Initial p
% tf = 10; % Final time
% tk = linspace(0,tf,nt);

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

% % Coefficients
% fx(:,:,1) = [
%        1     0      0;
%   -alpha -beta -delta
% ];
% % Exponents of x1
% fx(:,:,2) = [
%     0 0 0;
%     1 3 0
% ];
% % Exponents of x2
% fx(:,:,3) = [
%     1 0 0;
%     0 0 1;
% ];

%% Order of the basis functions
% c = 3;
% assert(c >= 3, 'Order of basis should be greater than 3 for this application.');

%% Number of basis functions
% all combinations of basis depending on the number of variables in f(x)
% and the chosen order of the basis.
ns = (c + 1)*(c + 2)/2;   

%% Legendre basis indexes
ind = zeros(ns,Nx);
s = 1;
for ord = 1:c
    for i2 = 0:c
        for i1 = 0:c
            if (i1+i2) == ord
                s = s + 1;
                ind(s,:) = [i1 i2];
            end
        end
    end
end


%% Definition of the Legendre polynomials
% order 0, 1, 2, ..., c
% so (c+1)
LPC = zeros(c+1,c+1);
LPC(1,1) = 1;
LPC(2,2) = 1;
for i = 3:c+1
    for j = 1:i-1
        LPC(i,j+1) = LPC(i,j+1) + (2*(i-2) + 1)/(i-1)*LPC(i-1,j);
        LPC(i,j) = LPC(i,j) - (i-2)/(i-1)*LPC(i-2,j);
    end
end


%% Legendre polynomials in multiple dimensions
% Pre-multiplication for the normalized constant
NLPC = zeros(c+1,c+1);
for i = 1:c+1
    for j = 1:c+1
        NLPC(i,j) = sqrt((2*(i-1)+1)/2)*LPC(i,j);
    end
end

% ORIGINAL CODE
% Multiplication of one dimensional Legendre polynomials
% MLP = zeros(ns,ns);
% for i = 1:ns
%     for j = 1:ns
%         MLP(i,j) = NLPC(ind(i,1)+1,ind(j,1)+1)*...
%                         NLPC(ind(i,2)+1,ind(j,2)+1);
%     end
% end

% GENERALIZED CODE: For Nx variables
% Multiplication of one dimensional Legendre polynomials
MLP = ones(ns,ns);
for i = 1:ns
    for j = 1:ns
        for k = 1:Nx
            MLP(i,j) = MLP(i,j) * NLPC(ind(i,k)+1,ind(j,k)+1);
        end
    end
end


%% Derivative of Normalized Legendre polynomials
DLPC = zeros(c+1,c+1);
for i = 2:c+1
    for j = 1:c
        DLPC(i,j) = j*NLPC(i,j+1);
    end
end

%% Operator Matrix
% GENERALIZED CODE: For Nx variables
K = zeros(ns,ns);
MDLP = zeros(Nx,ns);
DB = zeros(ns*Nx,Nx+1);
par = zeros(1,Nx);

% i represents the selected Legendre polynomials we are evaluating 
% the derivatives of
for i = 1:ns

    % Partial Derivatives of Legendre polynomials in multiple dimensions
    for j = 1:ns
        % Partial derivatives are chain rule so just multiply 
        % e.g., dL1dt * L2 * L3 * ... * L_{Nx}
        % or,   L1 * dL2dt * L3 * ... * L_{Nx}
        for k = 1:Nx
            MDLP(k,j) = DLPC(ind(i,k)+1,ind(j,k)+1);
            for l = 1:Nx
                if k ~= l
                    MDLP(k,j) = MDLP(k,j) * NLPC(ind(i,l)+1,ind(j,l)+1);
                end
            end
        end
    end

    % Total Derivative = fx * grad(Legendre Polynomials)
    DB = zeros(ns*2,Nx+1);
    s = 0;
    for dim = 1:Nx
        for ifx = 1:Nt  % loop thru all the terms in f(x)
            if (fx(dim,ifx,1) == 0)
                break;
            else
                for j = 1:ns
                    if (MDLP(dim,j) ~= 0)
                        s = s + 1;
                        DB(s,1) = fx(dim,ifx,1)*MDLP(dim,j);
                        % For each variable update the database
                        for k = 2:Nx+1
                            DB(s,k) = fx(dim,ifx,k) + ind(j,k-1);
                        end
                    end
                end
            end
        end
    end

    % Matix integration
    for k = 1:ns
        for j = 1:ns
            if (MLP(k,j) ~= 0)
                for ifx = 1:s
                    flag = 1;
                    for in = 1:Nx
                        % First part of integration where we get the
                        % (exponent+1) with which we divide the integrand
                        par(in) = round(ind(j,in) + DB(ifx,in+1) + 1);
                        % Check if exponent is "even" or "odd"
                        if (mod(par(in),2) == 0)
                            flag = 0;
                            break;
                        end
                    end
                    % If odd, perform the second part of the integration
                    if (flag == 1)
                        K(i,k) = K(i,k) + (2^Nx)*MLP(k,j)*DB(ifx,1) / (prod(par));
                    end
                end
            end
        end
    end
end



%% Observables
H = zeros(Nx,ns);
Obs = zeros(1,Nx+1);
for i = 1:Nx
    % Observable Polynomials
    Obs = zeros(1,Nx+1);
    Obs(1,1) = 1;
    Obs(1,i+1) = 1;

    % Matix integration
    for k = 1:ns
        for j = 1:ns
            if (MLP(k,j) ~= 0)
                flag = 1;
                for in = 1:Nx
                    par(in) = round(ind(j,in) + Obs(1,in+1) + 1);
                    if (mod(par(in),2) == 0)
                        flag = 0;
                        break;
                    end
                end
                if (flag == 1)
                    H(i,k) = H(i,k) + (2^Nx)*MLP(k,j)*Obs(1,1) / (prod(par));
                end
            end
        end
    end
end


%% Eigenvalue decomposition
[V,D] = eig(K,'nobalance');
iV = V \ eye(ns);


%% Generation of the functional space of solutions
g0 = zeros(ns,1);
arrExp = @(arr,exp) prod(arr .^ exp);
for i = 1:ns
    for j = 1:ns
        g0(i) = g0(i) + MLP(i,j) * arrExp([q0,p0],ind(j,:)); %    q0^ind(j,1)*p0^ind(j,2);
    end
end

%% Computaion of the solution as a function of time
Sol = zeros(Nx,nt);
for k = 1:nt
    Sol(:,k) = real(H * V * expm(D * tk(k)) * iV * g0);
end
q = Sol(1,:);
p = Sol(2,:);

%% ODE45 results for comparison
param.sm = sm;
param.sk = sk;
param.se = se;
param.sa = sa;
[t,states] = ode45(@(t,x) duffing(t,x,param), tk, [q0,p0]);


%% Plotting results
fig1 = figure(1);
scatter(p,q)
hold on; 
scatter(states(:,1),states(:,2),"+")
hold off

fig2 = figure(2);
plot(tk,p,lineWidth=2)
hold on
plot(tk,q,lineWidth=2)
plot(tk,states(:,1),'--',lineWidth=2)
plot(tk,states(:,2),'--',lineWidth=2)
hold off


%% Function
function dxdt = duffing(t,x,p)
    dxdt = zeros(2,1);
    x1 = x(1);
    x2 = x(2);
    dxdt(1) = x2/p.sm;
    dxdt(2) = -p.sk*x1 - p.sk*p.se*p.sa^2*x1^3;
end