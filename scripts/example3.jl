## Example 3: Van Der Pol Oscillator

using DifferentialEquations
using KoopOpGalerkin
using GLMakie
# Makie.inline!(true)  # uncomment to show in navigation pane

const KOG = KoopOpGalerkin

## Equation parameters for the Van Der Pol Oscillator
mu = 0.01

# Other settings
nt = 1000
x10 = 1.0  # Initial x1
x20 = 0.1  # Initial x2
x0 = [x10, x20]
tf = 100.0  # Final time
tk = range(0, tf, length = nt)

## Differential equation
Nx = 2  # Number of variables
Nt = 3  # Maximum number of terms in f(x)
fx = zeros(Float64, Nx, Nt, Nx + 1)

# Coefficients and exponents
fx[:, :, 1] = [1.0  0.0 0.0; -1.0 mu -mu]
fx[:, :, 2] = [0.0 0.0 0.0; 1.0 0.0 2.0]
fx[:, :, 3] = [1.0 0.0 0.0; 0.0 1.0 1.0]

## Koopman Solution
c = 3
KSol = KOG.KoopmanSolve(c, fx, x0, tk)
x1 = KSol[1, :]
x2 = KSol[2, :]

## ODE45 results for comparison
function vanderpol!(xdot, x, mu, t)
    x1 = x[1]
    x2 = x[2]
    xdot[1] = x2
    xdot[2] = -x1 + mu * x2 - mu * x1^2 * x2
    return xdot
end

prob = ODEProblem(vanderpol!, x0, (0.0, tf), mu)
sol = solve(prob, Tsit5(), saveat = tk)

## Plotting results
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
