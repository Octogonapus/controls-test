using DifferentialEquations
using LinearAlgebra
using ControlSystems
using Plots

g = 9.81
nom_v = 12
stall_t = 2.42
stall_c = 133
free_c = 2.7
free_s = 5310

R = nom_v / stall_c
Kₜ = stall_t / stall_c
Kᵥ = free_s / (nom_v - R * free_c)
m = 2.2675
l = 1.2192
J = (1/3) * m * l^2
G = 1 / 2

Gₗ = G
Gᵣ = G
r = 0.08255 / 2.0
rb = 0.59055 / 2.0
J = 6.0
m = 52

function diff_drive(x′, x, u, t)
    θ = x[3]
    vl = x[4]
    vr = x[5]
    v = (vl + vr) / 2

    C1 = -(Gₗ^2 * Kₜ)/(Kᵥ * R * r^2)
    C2 = (Gₗ * Kₜ) / (R * r)
    C3 = -(Gᵣ^2 * Kₜ)/(Kᵥ * R * r^2)
    C4 = (Gᵣ * Kₜ) / (R * r)

    A = [0 0 -v*sin(θ) cos(θ)/2 cos(θ)/2
         0 0 v*cos(θ) sin(θ)/2 sin(θ)/2
         0 0 0 -1/(2*rb) 1/(2*rb)
         0 0 0 (1/m + (rb^2)/J)*C1 (1/m - (rb^2)/J)*C3
         0 0 0 (1/m - (rb^2)/J)*C1 (1/m + (rb^2)/J)*C3]
    B = [0 0
         0 0
         0 0
         (1/m + (rb^2)/J)*C2 (1/m - (rb^2)/J)*C4
         (1/m - (rb^2)/J)*C2 (1/m + (rb^2)/J)*C4]

    lqr_Q = I
    lqr_R = I
    K = dlqr(A, B, lqr_Q, lqr_R)
    @show A
    @show B
    @show K
    @show u(x, K)

    x′ .= A * x + B * u(x, K)
end

function uff(x, K)
    θ = x[3]
    ft = [cos(θ) sin(θ) 0 0 0
          -sin(θ) cos(θ) 0 0 0
          0 0 1 0 0
          0 0 0 1 0
          0 0 0 0 1]
    ref = transpose([0 0 0 1 1])

    @show K * ft
    @show ref - x
    @show K * ft * (ref - x)
    K * ft * (ref - x)
end

x0 = transpose([0 0 0 0 0])
tspan = (0.0, 4.0)
prob = ODEProblem(diff_drive, x0, tspan, uff)
sol = solve(prob)
plot(sol,linewidth=2,xaxis="t",label=["x [m]" "y [m]" "θ [rad]" "vₗ [m/s]" "vᵣ [m/s]"],layout=(5,1),size=(640,480*1.8))
