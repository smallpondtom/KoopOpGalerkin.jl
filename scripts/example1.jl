## Example 1: Duffing Oscillator with Perturbation

using DifferentialEquations
using KoopOpGalerkin
using GLMakie
# Makie.inline!(true)  # uncomment to show in navigation pane

const KOG = KoopOpGalerkin
const NORMALIZE = false
const REGULARIZE = true

## Equation parameters
sm = 1.0  # Mass
sk = 0.3  # Spring constant
sa = 1.0  # Unit transformation constant
se = 0.001 # Small parameter

## Other settings
nt = 10000
q0 = 1.0  # Initial q
p0 = 0.0  # Initial p
x0 = [q0, p0]
tf = 30.0  # Final time
tk = range(0, tf, length = nt)

## Run ODE45 system to obtain the L∞ norm of x1 and x2
if NORMALIZE
    function duffing1!(dxdt, x, p, t)
        sm, sk, se, sa = p
        x1 = x[1]
        x2 = x[2]
        dxdt[1] = x2 / sm
        dxdt[2] = -sk * x1 - sk * se * sa^2 * x1^3
    end
    p_vector = [sm, sk, se, sa]
    prob = ODEProblem(duffing1!, x0, (0.0, tf), p_vector)
    sol = solve(prob, Tsit5(), saveat = tk)
    x1max = maximum(abs.(sol[1, :]))
    x2max = maximum(abs.(sol[2, :]))
    x0 = x0 ./ [x1max, x2max]
else
    x1max = 1
    x2max = 1
end

## Differential equation
Nx = 2  # Number of variables
Nt = 2  # Maximum number of terms in f(x)
fx = zeros(Float64, Nx, Nt, Nx + 1)

# Linear system
fx[1, 1, 1] = 1.0
fx[1, 1, 3] = 1.0
fx[2, 1, 1] = -sk
fx[2, 1, 2] = 1.0

# Perturbing term
fx[2, 2, 1] = -sk * se * sa^2
fx[2, 2, 2] = 3.0

# Normalize
fx[1, :, 1] ./= x1max
fx[2, :, 1] ./= x2max

## Koopman Solution
c = 3
KSol = KOG.KoopmanSolve(c, fx, x0, tk)
q = KSol[1, :]
p = KSol[2, :]

## Normalized system: ODE45 results for comparison
function duffing1!(dxdt, x, p, t)
    sm, sk, se, sa, x1max, x2max = p
    x1 = x[1]
    x2 = x[2]
    dxdt[1] = (x2 / sm) / x1max
    dxdt[2] = (-sk * x1 - sk * se * sa^2 * x1^3) / x2max
end
p_vector = [sm, sk, se, sa, x1max, x2max]
prob = ODEProblem(duffing1!, x0, (0.0, tf), p_vector)
sol = solve(prob, Tsit5(), saveat = tk)

## EDMD for comparison
# TODO: The EDMD results are incredibly off. Probably something wrong with the implementation.
# dsol = sol(tk, Val{1})  # generate derivative data
dsol = sol(range(tk.step.hi, tf+tk.step.hi, length=nt))

# Normalize the states and derivatives
lb_states = [minimum(sol[i,:]) for i = 1:Nx]
ub_states = [maximum(sol[i,:]) for i = 1:Nx]
states_n = KOG.normalize(sol[:,:], lb_states, ub_states, -1, 1)

lb_dstates = [minimum(dsol[i,:]) for i = 1:Nx]
ub_dstates = [maximum(dsol[i,:]) for i = 1:Nx]
dstates_n = KOG.normalize(dsol[:,:], lb_dstates, ub_dstates, -1, 1)

## Lift the states and derivatives
n_KO = 20
Xm_lift, Xp_lift, Ms = KOG.lift_data_RBF(states_n, dstates_n, n_KO; verbose=true)
if REGULARIZE
    Ã, ro = KOG.TREDMD(Xm_lift, Xp_lift, 10, 5)
else
    Ã, ro = KOG.EDMD(Xm_lift, Xp_lift, 0)
end

## Simulate the EDMD linearized system
function duffing1_EDMD!(dxdt, x, Ã, t)
    dxdt .= Ã * x
end
x0_lift = KOG.lift_state_RBF(x0, n_KO, Ms)
prob = ODEProblem(duffing1_EDMD!, x0_lift[1:ro], (0.0, tf), Ã)
sol_EDMD = solve(prob, Tsit5(), saveat = tk)

## Plotting results
# Create scatter plot for phase space
fig = Figure(size=(800, 600))
ax1 = Axis(fig[1, 1], xlabel = "x1", ylabel="x2")
scatter!(ax1, sol[1, :], sol[2, :], color=:transparent, strokewidth=1, 
        strokecolor=:red, label="Analytical", marker=:circle)
scatter!(ax1, KSol[1, :], KSol[2, :], color=:blue, label="Koopman", marker=:cross)
scatter!(ax1, sol_EDMD[1, :], sol_EDMD[2, :], color=:transparent,
        strokecolor=:green, strokewidth=1, label="EDMD", marker=:diamond)
axislegend(ax1, position=:lt)
ax1.title = "Phase Space"

# Create line plot for time series
ax2 = Axis(fig[2, 1])
lines!(ax2, tk, sol[1, :], color = :red, label = "x1")
lines!(ax2, tk, sol[2, :], color = :purple, label = "x2")
lines!(ax2, tk, KSol[1, :], color = :blue, label = "x1_k", linestyle=:dash)
lines!(ax2, tk, KSol[2, :], color = :green, label = "x2_k", linestyle=:dash)
lines!(ax2, tk, sol_EDMD[1, :], color = :orange, label = "x1_EDMD", linestyle=:dashdot, linewidth=2)
lines!(ax2, tk, sol_EDMD[2, :], color = :pink, label = "x2_EDMD", linestyle=:dashdot, linewidth=2)
axislegend(ax2)
ax2.title = "Time Series"
ax2.xlabel = "t"
ax2.ylabel = "x"

# Display the figure
fig