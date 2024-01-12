## Example 2: Duffing Oscillator (general form)

using DifferentialEquations
using KoopOpGalerkin
using GLMakie
# Makie.inline!(true)  # uncomment to show in navigation pane

const KOG = KoopOpGalerkin
const NORMALIZE = false
const REGULARIZE = false

## Equation parameters for the Duffing Oscillator
alpha = 1.0
beta = 2.0
delta = 0.4

# Other settings
nt = 1000
x10 = 1.0  # Initial x1
x20 = 0.0  # Initial x2
x0 = [x10, x20]
tf = 20.0  # Final time
tk = range(0, tf, length = nt)
dt = tk[2] - tk[1]

## Run ODE45 system to obtain the L∞ norm of x1 and x2
if NORMALIZE
    function duffing2noforce!(xdot, x, p, t)
        alpha, beta, delta = p
        x1 = x[1]
        x2 = x[2]
        xdot[1] = x2
        xdot[2] = -alpha * x1 - beta * x1^3 - delta * x2
        return xdot
    end
    param = [alpha, beta, delta]
    prob = ODEProblem(duffing2noforce!, x0, (0.0, tf), param)
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
Nt = 3  # Maximum number of terms in f(x)
fx = zeros(Float64, Nx, Nt, Nx + 1)

# Coefficients and exponents
fx[:, :, 1] = [1.0 0.0 0.0; -alpha -beta -delta]
fx[:, :, 2] = [0.0 0.0 0.0; 1.0 3.0 0.0]
fx[:, :, 3] = [1.0 0.0 0.0; 0.0 0.0 1.0]

# Normalize
fx[1, :, 1] ./= x1max
fx[2, :, 1] ./= x2max

## Koopman Solution
c = 8
KSol = KOG.KoopmanSolve(c, fx, x0, tk)
x1 = KSol[1, :]
x2 = KSol[2, :]

## Normalized System: ODE45 results for comparison
function duffing2noforce!(xdot, x, p, t)
    alpha, beta, delta, x1max, x2max = p
    x1 = x[1]
    x2 = x[2]
    xdot[1] = (x2) / x1max
    xdot[2] = (-alpha * x1 - beta * x1^3 - delta * x2) / x2max
    return xdot
end
param = [alpha, beta, delta, x1max, x2max]
prob = ODEProblem(duffing2noforce!, x0, (0.0, tf), param)
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
    Φ, Ω, Ã, ro = KOG.TREDMD(Xm_lift, Xp_lift, dt, 0.1, 5)
else
    Φ, Ω, Ã, ro = KOG.EDMD(Xm_lift, Xp_lift, dt, 3)
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

# Simulate the Ã matrix linear system from EDMD
# function duffing2noforce_EDMD!(dxdt, x, Ã, t)
#     dxdt .= Ã * x
# end
# x0_lift = KOG.lift_state_RBF(x0, n_KO, Ms)
# prob = ODEProblem(duffing2noforce_EDMD!, x0_lift[1:ro], (0.0, tf), Ã)
# sol_EDMD = solve(prob, Tsit5(), saveat = tk)

## Plotting results
fig = Figure(size=(800, 600))
ax1 = Axis(fig[1, 1], xlabel = "x1", ylabel="x2")
scatter!(ax1, sol[1, :], sol[2, :], color=:transparent, strokewidth=1, 
        strokecolor=:red, label="Analytical", marker=:circle)
scatter!(ax1, KSol[1, :], KSol[2, :], color=:blue, label="Koopman", marker=:cross)
scatter!(ax1, sol_EDMD[1, :], sol_EDMD[2, :], color=:transparent,
        strokecolor=:green, strokewidth=1, label="EDMD", marker=:diamond)
axislegend(ax1)
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