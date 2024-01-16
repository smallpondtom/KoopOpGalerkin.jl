## Example 6: Lorenz System

using DifferentialEquations
using KoopOpGalerkin
using GLMakie
# Makie.inline!(true)  # uncomment to show in navigation pane

const KOG = KoopOpGalerkin
const REGULARIZE = false
const DMD_DATASET_NUM = 10

## Equation parameters for the Lorenz System
sigma = 10;
rho = 28;
beta = 8/3;

# Other settings
nt = 10000
x10 = 1.0  # Initial x1
x20 = 2.1  # Initial x2
x30 = 3.0  # Initial x3
x0 = [x10, x20, x30]
tf = 100.0  # Final time
tk = range(0, tf, length = nt)
dt = tk[2] - tk[1]

## Differential equation
Nx = 3  # Number of variables
Nt = 3  # Maximum number of terms in f(x)
fx = zeros(Float64, Nx, Nt, Nx + 1)

# Coefficients and exponents
fx[:,:,1] = [
    -sigma  sigma   0;
       rho     -1  -1;
         1  -beta   0
];
# Exponents of x1
fx[:,:,2] = [
    1 0 0;
    1 0 1;
    1 0 0
];
# Exponents of x2
fx[:,:,3] = [
    0 1 0;
    0 1 0;
    1 0 0
];
# Exponents of x3
fx[:,:,4] = [
    0 0 0;
    0 0 1;
    0 1 0;
];


## Koopman Solution
c = 8
KSol = KOG.KoopmanSolve(c, fx, x0, tk)

## ODE45 results for comparison
function lorenz!(dxdt, x, p, t)
    sigma, rho, beta = p
    x1 = x[1]
    x2 = x[2]
    x3 = x[3]
    dxdt[1] = sigma * (x2 - x1)
    dxdt[2] = x1 * (rho - x3) - x2
    dxdt[3] = x1 * x2 - beta * x3
    return dxdt
end


prob = ODEProblem(lorenz!, x0, (0.0, tf), [sigma, rho, beta])
sol_tmp = solve(prob, Tsit5(), saveat = tk)
sol = reduce(hcat, sol_tmp.u)
sol_save = deepcopy(sol)
dmd_data = sol
dsol = sol_tmp(range(tk.step.hi, tf+tk.step.hi, length=nt))
dmd_data_dot = reduce(hcat, dsol.u)

# More data for DMD
for i = 1:DMD_DATASET_NUM
    prob = ODEProblem(lorenz!, 10*rand(3).-5, (0.0, tf), [sigma, rho, beta])
    sol_new = solve(prob, Tsit5(), saveat = tk)
    dmd_data = hcat(dmd_data, reduce(hcat, sol_new.u))
    dsol_new = sol_new(range(tk.step.hi, tf+tk.step.hi, length=nt))
    dmd_data_dot = hcat(dmd_data_dot, reduce(hcat, dsol_new.u))
end

## EDMD for comparison
# Normalize the states and derivatives
lb_states = [minimum(dmd_data[i,:]) for i = 1:Nx]
ub_states = [maximum(dmd_data[i,:]) for i = 1:Nx]
states_n = KOG.normalize(dmd_data[:,:], lb_states, ub_states, -1, 1)

lb_dstates = [minimum(dmd_data_dot[i,:]) for i = 1:Nx]
ub_dstates = [maximum(dmd_data_dot[i,:]) for i = 1:Nx]
dstates_n = KOG.normalize(dmd_data_dot[:,:], lb_dstates, ub_dstates, -1, 1)

## Lift the states and derivatives
n_KO = 50
Xm_lift, Xp_lift, Ms = KOG.lift_data_RBF(states_n, dstates_n, n_KO; verbose=true)
if REGULARIZE
    Φ, Ω, Ã, ro = KOG.TREDMD(Xm_lift, Xp_lift, dt, 100, 5)
else
    Φ, Ω, Ã, ro = KOG.EDMD(Xm_lift, Xp_lift, dt, 5)
end

## Simulate the EDMD linearized system
x0_lift = KOG.lift_state_RBF(x0, n_KO, Ms)
sol_EDMD = zeros(ComplexF64, ro, nt)
b = Φ \ x0_lift
for i = 1:nt
    sol_EDMD[:, i] = exp.(Ω * tk[i]) .* b
end
sol_EDMD = Φ * sol_EDMD
sol_EDMD = real.(sol_EDMD)  # remove imaginary part


## Plotting results
fig = Figure(size=(1200, 600))
f1 = fig[1, 1] = GridLayout()
f2 = fig[2, 1] = GridLayout()
ax1 = Axis3(f1[1,1], xlabel = "x1", ylabel="x2", zlabel="x3")
scatter!(ax1, sol_save[1, :], sol_save[2, :], sol_save[3, :], color=:transparent, strokewidth=1, 
        strokecolor=:red, label="Analytical", marker=:circle)
scatter!(ax1, KSol[1, :], KSol[2, :], KSol[3, :], color=:blue, label="Koopman", marker=:cross)
scatter!(ax1, sol_EDMD[1, :], sol_EDMD[2, :], sol_EDMD[3, :], color=:transparent,
        strokecolor=:green, strokewidth=1, label="EDMD", marker=:diamond)
axislegend(ax1)
ax1.title = "Phase Space"

# Create line plot for time series
ax21 = Axis(f2[1,1:5])
lines!(ax21, tk, sol_save[1, :], color = :red, label = "x1")
lines!(ax21, tk, KSol[1, :], color = :blue, label = "x1_k", linestyle=:dash)
lines!(ax21, tk, sol_EDMD[1, :], color = :black, label = "x1_edmd", linestyle=:dash)
axislegend(ax21)

ax22 = Axis(f2[1,6:10])
lines!(ax22, tk, sol_save[2, :], color = :red, label = "x2")
lines!(ax22, tk, KSol[2, :], color = :blue, label = "x2_k", linestyle=:dash)
lines!(ax22, tk, sol_EDMD[2, :], color = :black, label = "x2_edmd", linestyle=:dash)
axislegend(ax22)

ax23 = Axis(f2[1,11:15])
lines!(ax23, tk, sol_save[3, :], color = :red, label = "x3")
lines!(ax23, tk, KSol[3, :], color = :blue, label = "x3_k", linestyle=:dash)
lines!(ax23, tk, sol_EDMD[3, :], color = :black, label = "x3_edmd", linestyle=:dash)
axislegend(ax23)

linkxaxes!(ax21, ax22, ax23)
linkyaxes!(ax21, ax22, ax23)

ax21.ylabel = "x"
Label(f2[2,:], "t", justification=:center, width=1.5)

# Display the figure
fig
