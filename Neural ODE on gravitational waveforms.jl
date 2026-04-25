### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 85a9406c-c2d6-4ceb-b6ff-fa7c369af160
using Lux, LineSearches, ComponentArrays

# ╔═╡ bbf4a648-db8a-488a-b9fa-268ef0666ffe
using OrdinaryDiffEqLowOrderRK,SciMLSensitivity

# ╔═╡ 2bd3454f-ff05-41df-8940-513269ac97f5
using Optimization, OptimizationOptimJL

# ╔═╡ d958151e-bdb3-411b-9487-124d08d5b893
using Unitful, UnitfulAstro

# ╔═╡ 1be88be7-8568-4175-b4b8-9721e78685ce
using StaticArrays

# ╔═╡ 038b3e95-ab3f-4324-9e34-73f6f42f8ae2
using CairoMakie

# ╔═╡ dd147f9e-38d4-4aee-944b-6e65846b3a12
using CairoMakie: Axis

# ╔═╡ 27eabe02-8fc6-40ec-bb5a-3443aa69c44e
using Printf

# ╔═╡ d8508c3d-1369-424c-a260-4a790777ec27
using Random

# ╔═╡ 879b0551-cd0b-40d8-abc8-8d99ed5bbd88
using PlutoUI

# ╔═╡ 962c289a-c730-11ef-0d6d-b3795820aac8
md"""
# Training a Neural ODE to Model Gravitational Waveforms
"""

# ╔═╡ 786644e5-78ad-4aa3-983c-ce7e338f19f6
md"""
Original tutorial link: <https://lux.csail.mit.edu/stable/tutorials/advanced/1_GravitationalWaveForm>
"""

# ╔═╡ 044849e1-3e71-48c3-bac0-fbc5800a11c1
TableOfContents()

# ╔═╡ ccbbe6f2-b454-4ebc-8722-d361fd0ff6e2
md"""
## Utility functions
"""

# ╔═╡ 892ab4fe-d722-4be2-800e-c8bb0b27ea0c
function one2two(path, m₁, m₂)
    M = m₁ + m₂
    r₁ = m₂ / M .* path
    r₂ = -m₁ / M .* path
    return r₁, r₂
end

# ╔═╡ f4424715-54a1-458b-9574-358d08239b76
@views function soln2orbit(soln, model_params = nothing)
    if size(soln, 1) ∉ [2, 4]
        throw(ArgumentError("size(soln, 1) must be either 2 or 4"))
    end

    if size(soln, 1) == 2
        χ = soln[1,:]
        φ = soln[2,:]

        if length(model_params) != 3
            throw(ArgumentError("length(model_params) must be 3 when size(soln, 1) == 2"))
        end
        p, M, e = model_params
    else
        χ = soln[1,:]
        φ = soln[2,:]
        p = soln[3,:]
        e = soln[4,:]
    end

    r = @. p / (1 + e * cos(χ))
    x = @. r * cos(φ)
    y = @. r * sin(φ)

    orbit = vcat(x', y')
    return orbit
end

# ╔═╡ 01a7e22c-4ace-4653-91ae-dcb5632340ec
"""
    Dₜ(v::AbstractVector, Δt)

Approximate derivative with respect to t.
"""
function Dₜ(v::AbstractVector, Δt)
    weights = SVector(3/2, -2, 1/2)
    a = -weights' * v[begin:begin+2]
    b = (v[begin+2:end] - v[begin:end-2]) / 2
    c = weights' * v[end:-1:end-2]
    return [a; b; c] / Δt
end

# ╔═╡ 1a390f0c-fdff-482a-bb99-e01232796622
"""
    Dₜ²(v::AbstractVector, Δt)

Approximate second derivative with respect to t.
"""
function Dₜ²(v::AbstractVector, Δt)
    weights = SVector(2, -5, 4, -1)
    a = weights' * v[begin:begin+3]
    b = v[begin:end-2] - 2 * v[begin+1:end-1] + v[begin+2:end]
    c = weights' * v[end:-1:end-3]
    return [a; b; c] / Δt^2
end

# ╔═╡ b1616d12-ba27-4c57-861a-2a1fc11f3ce7
function orbit2tensor(orbit, component, mass = 1)
    x = orbit[1, :]
    y = orbit[2, :]

    lxx = x .^ 2
    lyy = y .^ 2
    lxy = x .* y
    trace = lxx .+ lyy

    if component[1] == 1 && component[2] == 1
        tmp = lxx - trace/3
    elseif component[1] == 2 && component[2] == 2
        tmp = lyy - trace/3
    else
        tmp = lxy
    end

    return mass .* tmp
end

# ╔═╡ 1a63dc3c-0a87-4303-a894-142ff80bc802
function h_22_quadrupole_components(Δt, orbit, component, mass = 1)
    M = orbit2tensor(orbit, component, mass) # mass tensor
    return 2 * Dₜ²(M, Δt)
end

# ╔═╡ 5e54004d-c021-4ffc-863e-1ad0f9c057a0
begin
"""
    h_22_quadrupole(Δt, orbit, mass = 1)

One-body quadrupole
"""
function h_22_quadrupole(Δt, orbit, mass = 1)
    h11 = h_22_quadrupole_components(Δt, orbit, (1, 1), mass)
    h22 = h_22_quadrupole_components(Δt, orbit, (2, 2), mass)
    h12 = h_22_quadrupole_components(Δt, orbit, (1, 2), mass)

    return h11, h12, h22
end

"""
    h_22_quadrupole(Δt, orbit1, mass1, orbit2, mass2)

Two-body quadrupole
"""
function h_22_quadrupole(Δt, orbit1, mass1, orbit2, mass2)
    return (
        h_22_quadrupole(Δt, orbit1, mass1)
        .+
        h_22_quadrupole(Δt, orbit2, mass2)
    )
end
end;

# ╔═╡ a63a0e20-d150-428c-9241-8c5d9c5f6b87
begin
"""
    h_22_strain(Δt, orbit)
One-body strain
"""
function h_22_strain(Δt::T, orbit) where {T}
    h11, h12, h22 = h_22_quadrupole(Δt, orbit)

    h₊ = h11 - h22
    hₓ = T(2) * h12

    scale = √(T(π) / 5)
    return scale * h₊, -scale * hₓ
end

"""
    h_22_strain(Δt, orbit1, mass1, orbit2, mass2)

Two-body strain
"""
function h_22_strain(Δt::T, orbit1, mass1, orbit2, mass2) where {T}
    if abs(mass1 + mass2 - 1.0) > 1e-12
        throw(ArgumentError("masses do not sum to unity"))
    end

    h11, h12, h22 = h_22_quadrupole(Δt, orbit1, mass1, orbit2, mass2)

    h₊ = h11 - h22
    hₓ = T(2) * h12

    scale = √(T(π) / 5)
    return scale * h₊, -scale * hₓ
end
end;

# ╔═╡ be1d332a-8e6e-4ff9-8b4c-b93e7b8c9da8
function compute_waveform(Δt::T, soln, mass_ratio, model_params = nothing) where {T}
    mass_ratio > 1 && throw(DomainError(mass_ratio, "mass_ratio must be ≤ 1"))
    mass_ratio < 0 && throw(DomainError(mass_ratio, "mass_ratio must be non-negative"))

    orbit = soln2orbit(soln, model_params)

    if mass_ratio > 0
        m₂ = inv(one(T) + mass_ratio)
        m₁ = mass_ratio * m₂

        orbit₁, orbit₂ = one2two(orbit, m₁, m₂)
        waveform = h_22_strain(Δt, orbit₁, m₁, orbit₂, m₂)
    else
        waveform = h_22_strain(Δt, orbit)
    end

    return waveform
end

# ╔═╡ 48fdc879-a635-4d0e-ace1-8b26e7a0bb4a
md"""
## Simulating the True Model
"""

# ╔═╡ c3e6e4c8-99a1-4750-b6e6-6f8b2f1f3c7a
function RelativisticOrbitModel(u, (; p, M, e), t)
    χ, φ = u

    scale = (
        (p - 2 - 2e * cos(χ)) * (1 + e * cos(χ))^2
        /
        √((p - 2)^2 - 4 * e^2)
    )

    χ̇ = scale * √(p - 6 - 2e * cos(χ)) / (M * p^2)
    φ̇ = scale / (M * p^(3//2))

    return [χ̇, φ̇]
end

# ╔═╡ e1d3e6fa-3137-42be-8202-6f9d1c00b375
const mass_ratio = 0.0; # test particle

# ╔═╡ d3c17256-db15-4009-beca-73cb0407aef5
const u₀ = @SVector Float64[π, 0]; # initial conditions

# ╔═╡ e53e1c46-3a2d-4a8e-b052-d788dc5bd8bb
const datasize = 250;

# ╔═╡ 22667f5e-3372-462c-a145-e6dfacd7a26a
const tspan = (0.0f0, 6.0f4); # timespace for GW waveform

# ╔═╡ 8cb9b29b-4f6a-472d-b24e-61b13aca4a26
const tsteps = range(tspan..., length = datasize);

# ╔═╡ abb72eb0-c0b7-407e-9d67-548a7d0671c2
const Δt_data = step(tsteps)

# ╔═╡ 83ee915f-f2c3-401f-ba76-066317a5bb6e
const Δt = 100.0;

# ╔═╡ ff64c23a-bf1e-47bb-ad85-f21644e8b8f8
const ode_model_params = (p = 100.0, M = 1.0, e = 0.5);

# ╔═╡ bd5eeec1-58b6-48f1-b09f-1d650f8e72e7
md"""
Let's simulate the true model and plot the results using OrdinaryDiffEq.jl
"""

# ╔═╡ 75a864ae-9cec-45d4-8fef-41d3bde1fd63
prob = ODEProblem(RelativisticOrbitModel, u₀, tspan, ode_model_params)

# ╔═╡ 420bac2b-d748-4142-b946-8142af8f6cd7
soln = solve(prob, RK4(), dt = Δt, saveat = tsteps, adaptive = false)

# ╔═╡ 4b2eb1e1-55b0-4dc5-ad24-485249ae11b0
waveform = compute_waveform(Δt_data, Array(soln), mass_ratio, ode_model_params) |> first

# ╔═╡ 213a5866-3d3e-4c16-9b05-d6e6605d9723
let
    fig = Figure()
    ax = Axis(fig[1,1],  xlabel = "Time", ylabel = "Waveform")

    l = lines!(ax, tsteps, waveform, linewidth = 2, alpha = 0.75)
    s = scatter!(ax, tsteps, waveform, marker = :circle, markersize = 10, alpha = 0.25)

    axislegend(ax, [[l, s]], ["Waveform Data"], position = :lt)

    fig
end

# ╔═╡ 603cf46b-9b3b-4212-8fdd-91282b417209
md"""
## Defining a Neural Network Model
"""

# ╔═╡ c822feb6-9fa1-498d-bb7e-6eea21a82ec6
initparams = (init_weight = truncated_normal(std = 1e-4), init_bias = zeros32)

# ╔═╡ 7bd4c04b-8108-417f-814d-ac9d17e41fb5
const nn = Chain(
    Base.Fix1(fast_activation, cos),
    Dense(1 => 32, cos; initparams...),
    Dense(32 => 32, cos; initparams...),
    Dense(32 => 2; initparams...)
)

# ╔═╡ cd960168-dedc-4423-b220-3af753281cf0
ps, st = Lux.setup(MersenneTwister(1000), nn) |> f64

# ╔═╡ 3f9e4d37-a32f-402d-b865-f134af3212cc
const params = ComponentArray(ps |> f64)

# ╔═╡ 46ec1269-92e1-4d77-9a6e-8ed4deec00b5
const nn_model = StatefulLuxLayer{true}(nn, nothing, st)

# ╔═╡ a59e794e-869d-41f4-b18f-cb75ab320df3
function ODE_model(u, nn_params, t)
    χ, φ = u

    # In this example we know that `st` is an empty NamedTuple so
    # we can safelty ignore it. However, in general, we should use `st`
    # to store the state of the neural network.
    y = 1 .+ nn_model([first(u)], nn_params)

    p, M, e = ode_model_params
    scale = (1 + e * cos(χ)^2) / (M * p^(3//2))

    χ̇ = scale * y[1]
    φ̇ = scale * y[2]

    return [χ̇, φ̇]
end

# ╔═╡ fb080d62-fc65-42ce-8055-2b33de730632
prob_nn = ODEProblem(ODE_model, u₀, tspan, ps)

# ╔═╡ 623eb0be-c3e7-400e-b98a-184cf122202a
soln_nn = solve(
    prob_nn, RK4(), u0 = u₀, p = ps,
    saveat = tsteps, dt = Δt, adaptive = false)

# ╔═╡ be6a07c1-d220-4908-9741-37b90e5df708
waveform_nn = compute_waveform(
    Δt_data, Array(soln_nn), mass_ratio, ode_model_params
) |> first

# ╔═╡ 7630650a-14fb-47b3-8362-352d2dd1990d
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time", ylabel = "Waveform")

    lineoptions = (; linewidth = 2, alpha = 0.75)
    markeroptions = (; marker = :circle, markersize = 12, alpha = 0.5, strokewidth = 2)

    l1 = lines!(ax, tsteps, waveform; lineoptions...)
    s1 = scatter!(ax, tsteps, waveform; markeroptions...)

    l2 = lines!(ax, tsteps, waveform_nn; lineoptions...)
    s2 = scatter!(ax, tsteps, waveform_nn; markeroptions...)

    axislegend(ax, [[l1, s1], [l2, s2]], ["Waveform Data", "Waveform Neural Net (Untrained)"], position = :lb)

    fig
end

# ╔═╡ 0ea60643-17fd-4b15-9881-eb04594ad7fb
md"""
## Setting up for Neural Network training
"""

# ╔═╡ 45c88449-f285-4361-a4ff-5f4149aa885e
const mseloss = MSELoss();

# ╔═╡ f8387f39-5051-425b-97a9-31de79b6bf3d
function loss(θ)
    pred = solve(
        prob_nn, RK4(), u0 = u₀, p = θ,
        saveat = tsteps, dt = Δt, adaptive = false)
    pred_waveform = compute_waveform(
        Δt_data, Array(pred), mass_ratio, ode_model_params
    ) |> first
    return mseloss(pred_waveform, waveform)
end

# ╔═╡ ac9ca6c8-af63-4fe8-9072-2ab3bf9c81e2
loss(ps)

# ╔═╡ b0342dc9-f714-436e-a305-8c74ad388c0d
loss(params)

# ╔═╡ 54291b2f-b249-4fa7-8fd4-9d792497ab2a
begin
    const losses = Float64[]

    function callback(θ, loss)
        push!(losses, loss)
        @printf("Training \t Iteration: %5d \t Loss: %.10f\n", θ.iter, l)
        return false
    end
end;

# ╔═╡ c19508f1-b633-4532-b8ad-bff6f618c9ed
md"""
## Training the Neural Network
"""

# ╔═╡ 42474de7-ab48-4c46-81be-062ee0628f84
md"""
Using BFGS optimizers
"""

# ╔═╡ 8dcf26eb-4f25-4e33-b39c-f8180fec4c98
const optf = Optimization.OptimizationFunction((x, p) -> loss(x), AutoZygote())

# ╔═╡ b2204893-c6ed-419f-80b6-7f5aac130304
const optprob = Optimization.OptimizationProblem(optf, params)

# ╔═╡ f4bef33e-7b03-4f8f-8744-adb30c4202d2
res = Optimization.solve(
    optprob, BFGS(initial_stepnorm = 0.01, linesearch = LineSearches.BackTracking());
    callback, maxiters = 1000
)

# ╔═╡ 0a38a17a-fefd-446f-b7ae-0f8bbaaa89e8
md"""
## Visiualizing the Results
"""

# ╔═╡ b1d576dc-fb97-4fe2-b4df-9d15deb1e23f
let
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "Iteration", ylabel = "Loss")

    lines!(ax, losses, linewidth = 4, alpha = 0.75)
    scatter!(ax, losses, marker = :circle, markersize = 12, strokewidth = 2)

    fig
end

# ╔═╡ Cell order:
# ╟─962c289a-c730-11ef-0d6d-b3795820aac8
# ╟─786644e5-78ad-4aa3-983c-ce7e338f19f6
# ╠═85a9406c-c2d6-4ceb-b6ff-fa7c369af160
# ╠═bbf4a648-db8a-488a-b9fa-268ef0666ffe
# ╠═2bd3454f-ff05-41df-8940-513269ac97f5
# ╠═d958151e-bdb3-411b-9487-124d08d5b893
# ╠═1be88be7-8568-4175-b4b8-9721e78685ce
# ╠═038b3e95-ab3f-4324-9e34-73f6f42f8ae2
# ╠═dd147f9e-38d4-4aee-944b-6e65846b3a12
# ╠═27eabe02-8fc6-40ec-bb5a-3443aa69c44e
# ╠═d8508c3d-1369-424c-a260-4a790777ec27
# ╠═879b0551-cd0b-40d8-abc8-8d99ed5bbd88
# ╠═044849e1-3e71-48c3-bac0-fbc5800a11c1
# ╟─ccbbe6f2-b454-4ebc-8722-d361fd0ff6e2
# ╠═892ab4fe-d722-4be2-800e-c8bb0b27ea0c
# ╠═f4424715-54a1-458b-9574-358d08239b76
# ╟─01a7e22c-4ace-4653-91ae-dcb5632340ec
# ╟─1a390f0c-fdff-482a-bb99-e01232796622
# ╠═b1616d12-ba27-4c57-861a-2a1fc11f3ce7
# ╠═1a63dc3c-0a87-4303-a894-142ff80bc802
# ╠═5e54004d-c021-4ffc-863e-1ad0f9c057a0
# ╠═a63a0e20-d150-428c-9241-8c5d9c5f6b87
# ╠═be1d332a-8e6e-4ff9-8b4c-b93e7b8c9da8
# ╟─48fdc879-a635-4d0e-ace1-8b26e7a0bb4a
# ╠═c3e6e4c8-99a1-4750-b6e6-6f8b2f1f3c7a
# ╠═e1d3e6fa-3137-42be-8202-6f9d1c00b375
# ╠═d3c17256-db15-4009-beca-73cb0407aef5
# ╠═e53e1c46-3a2d-4a8e-b052-d788dc5bd8bb
# ╠═22667f5e-3372-462c-a145-e6dfacd7a26a
# ╠═8cb9b29b-4f6a-472d-b24e-61b13aca4a26
# ╠═abb72eb0-c0b7-407e-9d67-548a7d0671c2
# ╠═83ee915f-f2c3-401f-ba76-066317a5bb6e
# ╠═ff64c23a-bf1e-47bb-ad85-f21644e8b8f8
# ╟─bd5eeec1-58b6-48f1-b09f-1d650f8e72e7
# ╠═75a864ae-9cec-45d4-8fef-41d3bde1fd63
# ╠═420bac2b-d748-4142-b946-8142af8f6cd7
# ╠═4b2eb1e1-55b0-4dc5-ad24-485249ae11b0
# ╟─213a5866-3d3e-4c16-9b05-d6e6605d9723
# ╟─603cf46b-9b3b-4212-8fdd-91282b417209
# ╠═7bd4c04b-8108-417f-814d-ac9d17e41fb5
# ╠═c822feb6-9fa1-498d-bb7e-6eea21a82ec6
# ╠═cd960168-dedc-4423-b220-3af753281cf0
# ╠═3f9e4d37-a32f-402d-b865-f134af3212cc
# ╠═46ec1269-92e1-4d77-9a6e-8ed4deec00b5
# ╠═a59e794e-869d-41f4-b18f-cb75ab320df3
# ╠═fb080d62-fc65-42ce-8055-2b33de730632
# ╠═623eb0be-c3e7-400e-b98a-184cf122202a
# ╠═be6a07c1-d220-4908-9741-37b90e5df708
# ╟─7630650a-14fb-47b3-8362-352d2dd1990d
# ╟─0ea60643-17fd-4b15-9881-eb04594ad7fb
# ╠═45c88449-f285-4361-a4ff-5f4149aa885e
# ╠═f8387f39-5051-425b-97a9-31de79b6bf3d
# ╠═ac9ca6c8-af63-4fe8-9072-2ab3bf9c81e2
# ╠═b0342dc9-f714-436e-a305-8c74ad388c0d
# ╠═54291b2f-b249-4fa7-8fd4-9d792497ab2a
# ╟─c19508f1-b633-4532-b8ad-bff6f618c9ed
# ╟─42474de7-ab48-4c46-81be-062ee0628f84
# ╠═8dcf26eb-4f25-4e33-b39c-f8180fec4c98
# ╠═b2204893-c6ed-419f-80b6-7f5aac130304
# ╠═f4bef33e-7b03-4f8f-8744-adb30c4202d2
# ╟─0a38a17a-fefd-446f-b7ae-0f8bbaaa89e8
# ╠═b1d576dc-fb97-4fe2-b4df-9d15deb1e23f
