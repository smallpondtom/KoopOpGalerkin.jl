## Example 1: Duffing Oscillator with Perturbation

using DifferentialEquations
using KoopOpGalerkin
using GLMakie
# Makie.inline!(true)  # uncomment to show in navigation pane

const KOG = KoopOpGalerkin

## Equation parameters
sm = 1.0  # Mass
sk = 1.0  # Spring constant
sa = 1.0  # Unit transformation constant
se = 0.001  # Small parameter

## Other settings
nt = 100
q0 = 1.0  # Initial q
p0 = 0.0  # Initial p
x0 = [q0, p0]
tf = 10.0  # Final time
tk = range(0, tf, length = nt)

## Differential equation
Nx = 2  # Number of variables
Nt = 2  # Maximum number of terms in f(x)
fx = zeros(Float64, Nx, Nt, Nx + 1)

## Linear system
fx[1, 1, 1] = 1.0
fx[1, 1, 3] = 1.0
fx[2, 1, 1] = -sk
fx[2, 1, 2] = 1.0

## Perturbing term
fx[2, 2, 1] = -sk * se * sa^2
fx[2, 2, 2] = 3.0

## Koopman Solution
c = 3
KSol = KOG.KoopmanSolve(c, fx, x0, tk)
q = KSol[1, :]
p = KSol[2, :]

## ODE45 results for comparison
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

## Plotting results
# Create scatter plot for phase space
fig = Figure(size=(800, 600))
ax1 = Axis(fig[1, 1], xlabel = "x1", ylabel="x2")
scatter!(ax1, sol[1, :], sol[2, :], color=:transparent, strokewidth=1, 
        strokecolor=:red, label="Analytical", marker=:circle)
scatter!(ax1, KSol[1, :], KSol[2, :], color=:blue, label="Koopman", marker=:cross)
axislegend(ax1)
ax1.title = "Phase Space"

# Create line plot for time series
ax2 = Axis(fig[2, 1])
lines!(ax2, tk, sol[1, :], color = :red, label = "x1")
lines!(ax2, tk, sol[2, :], color = :purple, label = "x2")
lines!(ax2, tk, KSol[1, :], color = :blue, label = "x1_k", linestyle=:dash)
lines!(ax2, tk, KSol[2, :], color = :green, label = "x2_k", linestyle=:dash)
axislegend(ax2)
ax2.title = "Time Series"
ax2.xlabel = "t"
ax2.ylabel = "x"

# Display the figure
fig