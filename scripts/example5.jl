## Example 5: Van Der Pol-Duffing 3D Oscillator

using DifferentialEquations
using KoopOpGalerkin
using GLMakie
# Makie.inline!(true)  # uncomment to show in navigation pane

const KOG = KoopOpGalerkin

## Equation parameters for the Van Der Pol Oscillator
β = 0.01
ϵ = 0.05
k = 0.4

# Other settings
nt = 1000
x10 = 1.0  # Initial x1
x20 = 0.1  # Initial x2
x30 = 0.0  # Initial x3
x0 = [x10, x20, x30]
tf = 100.0  # Final time
tk = range(0, tf, length = nt)

## Differential equation
Nx = 3  # Number of variables
Nt = 5  # Maximum number of terms in f(x)
fx = zeros(Float64, Nx, Nt, Nx + 1)

# Coefficients and exponents
fx[:, :, 1] = [
    1.0  0.0 0.0 0.0 0.0;
    -1.0 ϵ -k -ϵ β;
    1.0 -1.0 0.0 0.0 0.0
]
fx[:, :, 2] = [
    0.0 0.0 0.0 0.0 0.0;
    1.0 0.0 0.0 0.0 3.0;
    0.0 0.0 0.0 0.0 0.0
]
fx[:, :, 3] = [
    1.0 0.0 0.0 0.0 0.0;
    0.0 1.0 0.0 1.0 0.0;
    1.0 0.0 0.0 0.0 0.0
]
fx[:, :, 4] = [
    0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 1.0 0.0 0.0;
    0.0 1.0 0.0 0.0 0.0
]

## Koopman Solution
c = 6
KSol = KOG.KoopmanSolve(c, fx, x0, tk)
x1 = KSol[1, :]
x2 = KSol[2, :]

## ODE45 results for comparison
function vanderpolDuffing!(xdot, x, p, t)
    β, ϵ, k = p
    x1 = x[1]
    x2 = x[2]
    x3 = x[3]
    xdot[1] = x2
    xdot[2] = -x1 + ϵ * x2 - k * x3 - ϵ * x1^2 * x2 + β * x1^3
    xdot[3] = x2 - x3
    return xdot
end

prob = ODEProblem(vanderpolDuffing!, x0, (0.0, tf), [β, ϵ, k])
sol_tmp = solve(prob, Tsit5(), saveat = tk)
sol = reduce(hcat, sol_tmp.u)

## Plotting results
fig = Figure(size=(1200, 600))
f1 = fig[1, 1] = GridLayout()
f2 = fig[2, 1] = GridLayout()
ax1 = Axis3(f1[1,1], xlabel = "x1", ylabel="x2", zlabel="x3")
scatter!(ax1, sol[1, :], sol[2, :], sol[3, :], color=:transparent, strokewidth=1, 
        strokecolor=:red, label="Analytical", marker=:circle)
scatter!(ax1, KSol[1, :], KSol[2, :], KSol[3, :], color=:blue, label="Koopman", marker=:cross)
axislegend(ax1)
ax1.title = "Phase Space"

# Create line plot for time series
ax21 = Axis(f2[1,1:5])
lines!(ax21, tk, sol[1, :], color = :red, label = "x1")
lines!(ax21, tk, KSol[1, :], color = :blue, label = "x1_k", linestyle=:dash)
axislegend(ax21)

ax22 = Axis(f2[1,6:10])
lines!(ax22, tk, sol[2, :], color = :purple, label = "x2")
lines!(ax22, tk, KSol[2, :], color = :green, label = "x2_k", linestyle=:dash)
axislegend(ax22)

ax23 = Axis(f2[1,11:15])
lines!(ax23, tk, sol[3, :], color = :orange, label = "x3")
lines!(ax23, tk, KSol[3, :], color = :brown, label = "x3_k", linestyle=:dash)
axislegend(ax23)

linkxaxes!(ax21, ax22, ax23)
linkyaxes!(ax21, ax22, ax23)

ax21.ylabel = "x"
Label(f2[2,:], "t", justification=:center, width=1.5)

# Display the figure
fig