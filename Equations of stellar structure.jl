### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ ae5820a2-26c3-4bf7-9dd5-d136eaf35038
using CairoMakie

# ╔═╡ 93698b5a-1f3b-4c47-9e67-aff44f234c0f
using PlutoUI: TableOfContents, Dump

# ╔═╡ 5accb1c7-bfad-40ab-af03-5d82970a744e
using StaticArrays

# ╔═╡ 1e6a2dc7-ff39-4463-a6d3-7a87ecc6ae53
using BoundaryValueDiffEq

# ╔═╡ dc9d3ae3-5d6e-4507-92a6-821746b42d21
using ADTypes

# ╔═╡ 44d9b444-3ed6-4df0-866b-7d683be2d278
using ForwardDiff

# ╔═╡ 9a09cca6-4112-4ae8-a195-ab3ea297ef8b
using FiniteDiff

# ╔═╡ 81b35b32-2be1-4458-bf91-97780765a6fd
using RecursiveArrayTools

# ╔═╡ f3e73bf0-7a4b-473e-a5c0-7f99a43b95b8
using RecursiveArrayTools: ArrayPartition

# ╔═╡ 5550c36f-dcef-4d95-9c90-46e91d276d37
using OrdinaryDiffEqLowOrderRK: DP5, RK4

# ╔═╡ 1670914b-092e-4bff-a229-c2a000f5461e
using OrdinaryDiffEq

# ╔═╡ 54f413fa-3a90-4e76-a192-75afa6a5bd6d
using SciMLBase: TwoPointBVPFunction

# ╔═╡ 67dead86-da4e-482e-b87e-f87e51ff5e5b
begin
    using Unitful, UnitfulAstro
    using Unitful: FreeUnits
    using Unitful: m, kg, J, K, W, Pa
    using UnitfulAstro: Msun
end

# ╔═╡ 329fa343-39e2-4a0b-adef-9c590a18f19c
using PhysicalConstants.CODATA2022: σ, k_B, c_0 as c, G

# ╔═╡ da3e18e3-0cda-440b-932e-7c09988e28dc
md"""
# Equations of Stellar Structure
"""

# ╔═╡ 03118617-5cca-4400-9e87-1562d42d7d44
md"""
For a non-rotating star in hydrostatic equilibrium (i.e., no time dependent behavior), the four equations of stellar structure are as follows:

```math
\begin{align*}
\frac{dr}{dM_r} &= \frac{1}{4πr^2 ρ}, &
\frac{dT}{dM_r} &= -\frac{3}{256π^2 σ_\text{SB}}\frac{κ L_r}{T^3 r^4}, \\[1ex]
\frac{dP}{dM_r} &= -\frac{GM_r}{4πr^4}, &
\frac{dL_r}{dM_r} &= ε.
\end{align*}
```
"""

# ╔═╡ 9ed8de43-cd5a-4e40-a8bb-b33038d80a8d
md"""
The key variables are the radial coordinate ``r``; ``M_r``, the mass interior to ``r``; ``P``, the local pressure; ``T``, the local temperature; and ``L_r``, the total energy production interior to ``r``. Note that ``P`` and ``T`` are local, while ``M_r`` and ``L_r`` are cumulative. In addition to the key variables just listed, there are numerous auxiliary variables and coefficients used in the equations:

- The radiation constant, ``a = \frac{4}{c} σ_\text{SB} = \mathrm{7.565×{10}^{-16} \, J\,m^{-3}\,K^{-4}}``, where ``σ_\text{SB}`` is the Stefan--Boltzmann constant

- The local density, _ρ_:

  ```math
  ρ(P,T) = \mathrm{9.91\!×\!{10}^{-28} \, kg} \;\; \frac{P - \frac{4σ_\text{SB}}{3c}T^4}{k_\text{B}T},
  ```

- The mass absorption coefficient/opacity, _κ_:

  ```math
  κ(ρ, T) = κ_0 + κ_1 \left(\frac{ρ}{1\,\mathrm{kg/m^3}}\right) \left(\frac{T}{1\,\mathrm{K}}\right)^{-3.5},
  ```

  where
  - ``κ_0 = \mathrm{0.035 \, {m^2}{/kg}}``
  - ``κ_1 = \mathrm{6.44\!×\!{10}^{18} \, {m^2}/{kg}}``

- The specific luminosity, _ε_:
  ```math
   ε(ρ, T) = \mathrm{0.136 \, W\,kg^{-1}} \left(\tilde{ε}_0 + 3.49\!×\!{10}^{12} \tilde{ε}_1 \right),
  ```

  where, in turn,

  ```math
  \begin{align*}
  \tilde{ε}_0 &= \left(\frac{ρ}{\mathrm{1\,kg/m^3}}\right) \frac{e^{-33.80 / T_6^{1/3}}}{T_6^{2/3}}, &
  \tilde{ε}_1 &= \left(\frac{ρ}{\mathrm{1\,kg/m^3}}\right)^2 \frac{e^{-4403/T_6}}{T_6^3}.
  \end{align*}
  ```
"""

# ╔═╡ d61d868e-7edd-46fc-9994-9c1af60b4cc9
md"""
The above formulae all use SI units, _T_₆ ≡ _T_ / (10⁶ K), and the constants ``c`` and ``k_\text{B}`` have their usual values. These have been calculated assuming a primordial star that is 75% hydrogen and 25% helium by mass, and which is fully ionized throughout. This problem is another boundary value problem. Take a primordial Population-III (i.e. pure H/He) star whose mass is ``100 M_⊙``. At its surface you have ``M_r = M_* = 100 M_⊙``, and ``P = 0``. At its core you have ``L_r = r = 0``.
"""

# ╔═╡ a4b9c073-39de-4097-8b30-b0b14617e843
md"""
Boundary conditions/values (known and unknown):

| Quantity | Inner boundary, ``M_r = 0`` | Outer boundary, ``M_r = M_⋆`` |
|:---------|:----------------------------|:------------------------------|
| ``r``    | ``0``                       | ``R_⋆ = ?``                   |
| ``P``    | ``P_c = ?``                 | ``0``                         |
| ``L_r``  | ``0``                       | ``L_⋆ = ?``                   |
| ``T``    | ``T_c = ?``                 |                               |
``P_c``, ``T_c``, ``R_⋆``, and ``L_⋆`` are unknown.
"""

# ╔═╡ a1288d60-52e2-49e4-96b2-46d0ea3800a1
"""Boundary conditions at Mᵣ = 0"""
function stellarboundaryinner(u_a, p)
    r, T, P, Lᵣ = u_a
    return [r, Lᵣ]
end

# ╔═╡ 0b78405a-ad60-44b8-a2b4-3fc57d866540
"""Boundary conditions at Mᵣ = 0"""
function stellarboundaryinner!(residual, u_a, p)
    r, T, P, Lᵣ = u_a
    residual .= (r, Lᵣ)
    return
end

# ╔═╡ 7ddbaf73-430e-462c-8f04-8299d2ef0c87
"""Boundary conditions at ``M_r = M_⋆``"""
function stellarboundaryouter(u_b, p)
    r, T, P, Lᵣ = u_b
    return P
end

# ╔═╡ 54cf2228-0041-476c-b56d-055b39afc596
"""Boundary conditions at ``M_r = M_⋆``"""
function stellarboundaryouter!(residual, u_b, p)
    r, T, P, Lᵣ = u_b
    residual[1] = P
    return
end

# ╔═╡ 04d0caf4-f307-47e5-9c2d-7b97da6e8f80
md"""
## Part a
"""

# ╔═╡ 5885df68-9573-4651-9a05-5d0a5b18d747
md"""
Use the shooting method to estimate the core temperature and pressure of this star, to two significant figures each. I suggest you use RK4 as your ODE solver, with 100 steps between ``M_r = 0`` and ``M_r = M_*``. Your goal is to get ``P`` close to 0 at the surface without going negative (since pressure is a positive quantity). To guide your initial guess, the Sun's core temperature is 16 × 10⁶ K, and its core pressure is 1.2 × 10¹⁶ Pa.
"""

# ╔═╡ 13b416d0-39ef-40bb-8997-8e1b7226e10a
solver = MultipleShooting(4, RK4())

# ╔═╡ 92fd50e9-dfdf-447c-ac60-fb46d0080dbc
u₀_guess = ArrayPartition(0.1m, 16.0e6K, 1.2e16Pa, 0.1W)

# ╔═╡ 6831dbe0-d11b-490e-b31b-211e4767b126
u₀_guess_nounits = ustrip.(u₀_guess) |> collect

# ╔═╡ 55fa150f-3216-42ff-93a8-a61aee84b9c6
T_c_guesses = logrange(3.0e5, 1.0e10, length = 3)

# ╔═╡ 892f5a7b-83ec-4eb9-a11f-6a86b04ee4e4
P_c_guesses = logrange(1.0e12, 1.0e20, length = 3)

# ╔═╡ 0d32114c-f8a7-4a09-8a8a-23e18f89b259
md"""
## Part b
"""

# ╔═╡ ce125c69-bc40-42d0-9ded-937257c9bc39
md"""
What are the radius and luminosity of the star, in units of solar radius (``R_⊙ = 6.96 × 10^8 \,\mathrm{m}``) and solar luminosity (``L_⊙ = 3.83 × 10^{26} \,\mathrm{W}``)?
"""

# ╔═╡ d5d6964c-0f6c-4d20-b97a-c6b8aab857d7
md"""
## Part c
"""

# ╔═╡ d7db0d1a-92fc-4a43-8285-7438f8d6a8e1
md"""
Plot ``r(M)``, ``T(M)``, and ``L(M)`` from ``0`` to ``M_r`` (It will probably be easier to read if you use three different plots rather than one panel for all three curves.) Use logarithmic y axis and a linear x axis.
"""

# ╔═╡ 6a99bdb9-7f3b-4292-8d9b-1c3c01f05e6e
axisproperties = (yscale = log10, xlabel = "Mᵣ (kg)", xminorgridvisible = true, yminorgridvisible = true, xminorticksvisible = true, yminorticksvisible = true)

# ╔═╡ a105e4e9-f732-4080-b2ea-bfce5c739a08
md"""
## Auxiliary functions
"""

# ╔═╡ fd286a3b-b45c-4891-a1b3-7cf0ca6209fe
"""
    massabsorpcoeff(ρ, T)

Mass absorption coefficient
"""
function massabsorpcoeff(ρ, T)
    κ₀ = 0.035m^2/kg
    κ₁ = 6.44e18m^2/kg
    # scale_factor = ustrip(kg/m^3, ρ) / ustrip(K, T)^3.5
    scale_factor = ustrip(kg / m^3, ρ) / ustrip(K, T)^4
    return κ₀ + κ₁ * scale_factor
end

# ╔═╡ a6728349-7323-4e11-a77b-3b05e61992a2
"""
    specificpower(ρ, T)

Specific power ``ε(ρ, T)`` as a function of local density and temperature.
"""
function specificpower(ρ, T)
    ρ = ustrip(kg / m^3, ρ)
    T₆ = ustrip(K, T) / 1.0e6 |> abs
    ε̃₀ = ρ / T₆^(2 // 3) * exp(-33.8 / ∛T₆)
    ε̃₁ = ρ^2 / T₆^3 * exp(-4403 / T₆)
    return 0.136W/kg * (ε̃₀ + 3.49e12 * ε̃₁)
end

# ╔═╡ 1b9eff24-67fa-4b64-9c3e-17288bf385f1
md"""
# Notebook setup
"""

# ╔═╡ f64a51e1-6b82-46e7-b7c9-4050fcc89f06
md"""
## Imports
"""

# ╔═╡ 7066ef95-8257-4950-a0a7-c15832e1725e
TableOfContents()

# ╔═╡ 45caddb4-8dad-451a-a6dc-fcda9e087dfe
md"""
## Constants reference
"""

# ╔═╡ cd9324f0-87be-4e56-b440-955f1fae11ea
const a = 4σ / c |> (J/m^3/K^4)

# ╔═╡ 58810d44-0b8f-4f74-83b9-9b95605cd2e7
σ

# ╔═╡ b0a76914-5907-49ad-9b7c-eb3a3e4614e7
k_B

# ╔═╡ ba8a7b26-cba3-411d-ba6c-7e8fcef66c58
G

# ╔═╡ b44185c8-9325-4aa9-a448-51ff08d95187
c

# ╔═╡ c68905f0-23ce-498d-a478-e6e41c19308b
md"""
Values specific to the problem:
"""

# ╔═╡ 817a6fc5-b327-48ea-9cda-4a4217931cc9
const M₀ = 9.91e-28kg

# ╔═╡ 9e69b9d6-08e1-4a83-a5ea-102620bdd674
density(P, T) = M₀ * (P - a*T^4/3) / (k_B*T)

# ╔═╡ 3806c7de-e09f-40db-aad5-aa7923b161fb
function odefun(u, p, Mᵣ)
    r, T, P, Lᵣ = u
    # @debug("Working with", Mᵣ, r, T, P, Lᵣ)

    ρ = density(P, T)
    κ = massabsorpcoeff(ρ, T)

    r′ = 1 / (4π * r^2 * ρ) |> m / kg
    T′ = -3 / (64*π^2*a*c) * (κ * Lᵣ) / (T^3 * r^4) |> K / kg
    P′ = - G * Mᵣ / (4π * r^4) |> Pa / kg
    Lᵣ′ = specificpower(ρ, T) |> W / kg

    u′ = ArrayPartition(r′, T′, P′, Lᵣ′)

    return u′
end

# ╔═╡ 7f3668b9-9baa-4c32-9854-420d4e65a5f3
bvpfunction_oop = TwoPointBVPFunction(odefun, (stellarboundaryinner, stellarboundaryouter))

# ╔═╡ 366f153f-e666-45f1-b594-28d0b8381830
function odefun!(u′, u, p, Mᵣ)
    r, T, P, Lᵣ = u
    # @debug("Working with", r, T, P, Lᵣ)

    ρ = density(P, T)
    κ = massabsorpcoeff(ρ, T)

    u′[1] = r′ = 1 / (4π * r^2 * ρ) |> m / kg
    u′[2] = T′ = -3 / (64 * π^2 * a * c) * (κ * Lᵣ) / (T^3 * r^4) |> K / kg
    u′[3] = P′ = - G * Mᵣ / (4π * r^4) |> Pa / kg
    u′[4] = Lᵣ′ = specificpower(ρ, T) |> W / kg

    return # u′
end

# ╔═╡ 6ac5b3da-c540-48ea-95b3-b51719ec36b1
# bvpfunction = BVPFunction(odefun, stellarboundary)
bvpfunction = TwoPointBVPFunction(
    odefun!, (stellarboundaryinner!, stellarboundaryouter!);
    bcresid_prototype = (ArrayPartition(0.0m, 0.0W), [0.0Pa]),
)

# ╔═╡ 4a70bdfa-905e-48d2-baee-0fb8e30b43f3
function odefun_nounits!(u′, u, _, Mᵣ)
    r, T, P, Lᵣ = u
    # @debug("Working with", r, T, P, Lᵣ)
    r, T, P, Lᵣ = r * m, T * K, P * Pa, Lᵣ * W
    Mᵣ = Mᵣ * kg

    ρ = density(P, T)
    κ = massabsorpcoeff(ρ, T)

    u′[1] = r′ = ustrip(m / kg, 1 / (4π * r^2 * ρ))
    u′[2] = T′ = ustrip(K / kg, -3 / (64 * π^2 * a * c) * (κ * Lᵣ) / (T^3 * r^4))
    u′[3] = P′ = ustrip(Pa / kg, - G * Mᵣ / (4π * r^4))
    u′[4] = Lᵣ′ = ustrip(W / kg, specificpower(ρ, T))

    return # u′
end

# ╔═╡ 56048471-a8d0-4d3d-896d-1ab02424c14a
bvpfunction_nounits = TwoPointBVPFunction(
    odefun_nounits!, (stellarboundaryinner!, stellarboundaryouter!);
    bcresid_prototype = ([0.0, 0.0], [0.0]),
)

# ╔═╡ 2267ba11-3424-48a5-b219-07ee3b2bc4d9
const M_star = 100Msun

# ╔═╡ e74e279d-8f51-4c3b-88a6-efffe6a9707e
Mᵣ_domain = (0.0kg, M_star |> kg)
# Mᵣ_domain = (0.0, ustrip(kg, M_star))

# ╔═╡ d8af8e17-d953-41ec-ad1e-ed7c246e23ea
bvproblem = BVProblem(bvpfunction, u₀_guess, Mᵣ_domain)

# ╔═╡ aea355d7-d5e8-44be-82fa-cf3b6f92a9d3
bvproblem.u0

# ╔═╡ fccdb8c0-0f63-49c0-bff9-d51c67b95a38
zero(bvproblem.u0)

# ╔═╡ bd4dcf7c-3189-4b0b-a8ff-0cc453dc78a3
solve(bvproblem, solver; adaptive = false, dt = 1.0e20kg)

# ╔═╡ e656b9ba-d115-49d1-a908-6ff14fc5ed6f
bvproblem_nounits = BVProblem(bvpfunction_nounits, u₀_guess_nounits, ustrip.(kg, Mᵣ_domain))

# ╔═╡ aa521e22-d0e6-4775-860c-186409716d15
solve(bvproblem_nounits, Shooting(Tsit5()), dt = 1.0e31, isoutofdomain = (u, _, Mᵣ) -> any(<(0), u))

# ╔═╡ 1dc16800-a81e-415d-b1c2-14c106ef23b4
bvproblem_oop = BVProblem(bvpfunction_oop, u₀_guess, Mᵣ_domain)

# ╔═╡ e22a514f-aea5-4634-a036-6d6cc41776bc
solve(bvproblem_oop, MIRK4(), dt = 1.0e29kg)

# ╔═╡ 7e96c4d4-c2a2-4887-9161-3f5ce9217e63
cld(Mᵣ_domain[2] - Mᵣ_domain[1], 1.0e29kg)

# ╔═╡ 55d47962-3653-432b-a270-96dc49478c77
cld(Mᵣ_domain[2] - Mᵣ_domain[1], 1.0e23kg)

# ╔═╡ 7e33115d-4ccc-4225-8887-4e42613b1685
(Mᵣ_domain[2] - Mᵣ_domain[1]) / 10

# ╔═╡ 6c25b9f7-e709-45c5-b4ed-86a1be0764b9
tstops = logrange(∛nextfloat(ustrip(kg, Mᵣ_domain[1])), ustrip(kg, Mᵣ_domain[2]), length = 1000)kg # |> collect

# ╔═╡ 7118c7f2-07c3-4e79-98b1-6bdab1ccd2c5
M_saves = range(Mᵣ_domain..., length = 10)

# ╔═╡ 79c2f0f9-5cb4-4a19-ac97-14d174825c50
solve(bvproblem, solver, isoutofdomain = (u, p, Mᵣ) -> any(<(0), u), abstol = 0.1, tstops = M_saves, adaptive = false)

# ╔═╡ d4c3d589-617f-4ee8-a254-fab99978ba08
begin
    solutions = Matrix{Any}(undef, length(T_c_guesses), length(P_c_guesses))
    for (T_idx, T_c) in enumerate(T_c_guesses), (P_idx, P_c) in enumerate(P_c_guesses)
        # create ODE IVP
        local u₀_guess = [1.0e-9, T_c, P_c, 1.0e-10]
        bvp = BVProblem(bvpfunction, u₀_guess, Mᵣ_domain)
        # try solving it
        try
            sol = solve(bvp, solver; tstops = range(Mᵣ_domain..., length = 100), adaptive = false)
            # if it's solved, save it to the solution matrix, the index corresponding to this T_c_guess and p_c_guess
            solutions[T_idx, P_idx] = sol
        catch # if it errors, its not solvable
            solutions[T_idx, P_idx] = nothing
        end
    end
end

# ╔═╡ f85a4f6f-9b10-4520-81fa-a1c6cede6f4b
solutions

# ╔═╡ 30039996-5922-4fb9-9ee1-22f98f2685ba
solutions

# ╔═╡ 3d0a13c4-a1df-401f-84df-37d1499af842
ivp = ODEProblem(odefun!, u₀_guess, Mᵣ_domain)

# ╔═╡ 9d0e470d-51f5-4ec2-a49d-99827a08fab0
sol = solve(ivp, RK4())

# ╔═╡ f9b723c3-e720-4c1e-a1a2-1d15648b6300
lines(sol, axis = (; yscale = log10))

# ╔═╡ 919784e5-74b8-406b-b449-e0c3004cecfe
lines(sol.t, sol[1, :], axis = (; axisproperties..., ylabel = "r (m)"))

# ╔═╡ 088afba7-5757-49fe-88c9-850396d882dc
lines(sol.t, sol[2, :], axis = (; axisproperties..., ylabel = "T (K)"))

# ╔═╡ 2d17f6d3-eb87-4454-95b3-ca1a5dc4f52f
lines(sol.t, sol[3, :], axis = (; axisproperties..., ylabel = "P (Pa)"))

# ╔═╡ 7106980e-6555-402c-ba2b-66330e01f342
lines(sol.t, sol[4, :], axis = (; axisproperties..., ylabel = "Lᵣ (W)"))

# ╔═╡ Cell order:
# ╟─da3e18e3-0cda-440b-932e-7c09988e28dc
# ╟─03118617-5cca-4400-9e87-1562d42d7d44
# ╠═3806c7de-e09f-40db-aad5-aa7923b161fb
# ╠═366f153f-e666-45f1-b594-28d0b8381830
# ╠═4a70bdfa-905e-48d2-baee-0fb8e30b43f3
# ╟─9ed8de43-cd5a-4e40-a8bb-b33038d80a8d
# ╟─d61d868e-7edd-46fc-9994-9c1af60b4cc9
# ╠═e74e279d-8f51-4c3b-88a6-efffe6a9707e
# ╟─a4b9c073-39de-4097-8b30-b0b14617e843
# ╠═a1288d60-52e2-49e4-96b2-46d0ea3800a1
# ╠═0b78405a-ad60-44b8-a2b4-3fc57d866540
# ╠═7ddbaf73-430e-462c-8f04-8299d2ef0c87
# ╠═54cf2228-0041-476c-b56d-055b39afc596
# ╟─04d0caf4-f307-47e5-9c2d-7b97da6e8f80
# ╟─5885df68-9573-4651-9a05-5d0a5b18d747
# ╠═6ac5b3da-c540-48ea-95b3-b51719ec36b1
# ╠═56048471-a8d0-4d3d-896d-1ab02424c14a
# ╠═7f3668b9-9baa-4c32-9854-420d4e65a5f3
# ╠═d8af8e17-d953-41ec-ad1e-ed7c246e23ea
# ╠═e656b9ba-d115-49d1-a908-6ff14fc5ed6f
# ╠═aa521e22-d0e6-4775-860c-186409716d15
# ╠═1dc16800-a81e-415d-b1c2-14c106ef23b4
# ╠═7e96c4d4-c2a2-4887-9161-3f5ce9217e63
# ╠═e22a514f-aea5-4634-a036-6d6cc41776bc
# ╠═13b416d0-39ef-40bb-8997-8e1b7226e10a
# ╠═aea355d7-d5e8-44be-82fa-cf3b6f92a9d3
# ╠═fccdb8c0-0f63-49c0-bff9-d51c67b95a38
# ╠═bd4dcf7c-3189-4b0b-a8ff-0cc453dc78a3
# ╠═55d47962-3653-432b-a270-96dc49478c77
# ╠═7e33115d-4ccc-4225-8887-4e42613b1685
# ╠═6c25b9f7-e709-45c5-b4ed-86a1be0764b9
# ╠═7118c7f2-07c3-4e79-98b1-6bdab1ccd2c5
# ╠═92fd50e9-dfdf-447c-ac60-fb46d0080dbc
# ╠═6831dbe0-d11b-490e-b31b-211e4767b126
# ╠═55fa150f-3216-42ff-93a8-a61aee84b9c6
# ╠═892f5a7b-83ec-4eb9-a11f-6a86b04ee4e4
# ╠═d4c3d589-617f-4ee8-a254-fab99978ba08
# ╠═f85a4f6f-9b10-4520-81fa-a1c6cede6f4b
# ╠═3d0a13c4-a1df-401f-84df-37d1499af842
# ╠═9d0e470d-51f5-4ec2-a49d-99827a08fab0
# ╠═30039996-5922-4fb9-9ee1-22f98f2685ba
# ╠═79c2f0f9-5cb4-4a19-ac97-14d174825c50
# ╟─0d32114c-f8a7-4a09-8a8a-23e18f89b259
# ╟─ce125c69-bc40-42d0-9ded-937257c9bc39
# ╟─d5d6964c-0f6c-4d20-b97a-c6b8aab857d7
# ╟─d7db0d1a-92fc-4a43-8285-7438f8d6a8e1
# ╠═ae5820a2-26c3-4bf7-9dd5-d136eaf35038
# ╠═6a99bdb9-7f3b-4292-8d9b-1c3c01f05e6e
# ╠═f9b723c3-e720-4c1e-a1a2-1d15648b6300
# ╠═919784e5-74b8-406b-b449-e0c3004cecfe
# ╠═088afba7-5757-49fe-88c9-850396d882dc
# ╠═2d17f6d3-eb87-4454-95b3-ca1a5dc4f52f
# ╠═7106980e-6555-402c-ba2b-66330e01f342
# ╟─a105e4e9-f732-4080-b2ea-bfce5c739a08
# ╠═9e69b9d6-08e1-4a83-a5ea-102620bdd674
# ╠═fd286a3b-b45c-4891-a1b3-7cf0ca6209fe
# ╠═a6728349-7323-4e11-a77b-3b05e61992a2
# ╟─1b9eff24-67fa-4b64-9c3e-17288bf385f1
# ╟─f64a51e1-6b82-46e7-b7c9-4050fcc89f06
# ╠═93698b5a-1f3b-4c47-9e67-aff44f234c0f
# ╠═7066ef95-8257-4950-a0a7-c15832e1725e
# ╠═5accb1c7-bfad-40ab-af03-5d82970a744e
# ╠═1e6a2dc7-ff39-4463-a6d3-7a87ecc6ae53
# ╠═dc9d3ae3-5d6e-4507-92a6-821746b42d21
# ╠═44d9b444-3ed6-4df0-866b-7d683be2d278
# ╠═9a09cca6-4112-4ae8-a195-ab3ea297ef8b
# ╠═81b35b32-2be1-4458-bf91-97780765a6fd
# ╠═f3e73bf0-7a4b-473e-a5c0-7f99a43b95b8
# ╠═5550c36f-dcef-4d95-9c90-46e91d276d37
# ╠═1670914b-092e-4bff-a229-c2a000f5461e
# ╠═54f413fa-3a90-4e76-a192-75afa6a5bd6d
# ╠═67dead86-da4e-482e-b87e-f87e51ff5e5b
# ╠═329fa343-39e2-4a0b-adef-9c590a18f19c
# ╟─45caddb4-8dad-451a-a6dc-fcda9e087dfe
# ╠═cd9324f0-87be-4e56-b440-955f1fae11ea
# ╠═58810d44-0b8f-4f74-83b9-9b95605cd2e7
# ╠═b0a76914-5907-49ad-9b7c-eb3a3e4614e7
# ╠═ba8a7b26-cba3-411d-ba6c-7e8fcef66c58
# ╠═b44185c8-9325-4aa9-a448-51ff08d95187
# ╟─c68905f0-23ce-498d-a478-e6e41c19308b
# ╠═817a6fc5-b327-48ea-9cda-4a4217931cc9
# ╠═2267ba11-3424-48a5-b219-07ee3b2bc4d9
