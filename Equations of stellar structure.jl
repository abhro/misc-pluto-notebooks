### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ ae5820a2-26c3-4bf7-9dd5-d136eaf35038
using CairoMakie

# ╔═╡ 93698b5a-1f3b-4c47-9e67-aff44f234c0f
using PlutoUI: TableOfContents

# ╔═╡ 1e6a2dc7-ff39-4463-a6d3-7a87ecc6ae53
using BoundaryValueDiffEq

# ╔═╡ 44d9b444-3ed6-4df0-866b-7d683be2d278
using ForwardDiff

# ╔═╡ 81b35b32-2be1-4458-bf91-97780765a6fd
using RecursiveArrayTools

# ╔═╡ 5550c36f-dcef-4d95-9c90-46e91d276d37
using OrdinaryDiffEq: DP5

# ╔═╡ 54f413fa-3a90-4e76-a192-75afa6a5bd6d
using SciMLBase: TwoPointBVPFunction

# ╔═╡ 67dead86-da4e-482e-b87e-f87e51ff5e5b
begin
    using Unitful, UnitfulAstro
    using Unitful: m, kg, J, K, W, Pa
    using UnitfulAstro: Msun
end

# ╔═╡ 329fa343-39e2-4a0b-adef-9c590a18f19c
using PhysicalConstants.CODATA2018: σ, k_B, c_0 as c, G

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
\frac{dT}{dM_r} &= -\frac{3}{64π^2 ac}\frac{κ L_r}{T^3 r^4}, \\[1ex]
\frac{dP}{dM_r} &= -\frac{GM_r}{4πr^4}, &
\frac{dL_r}{dM_r} &= ε.
\end{align*}
```
"""

# ╔═╡ 9ed8de43-cd5a-4e40-a8bb-b33038d80a8d
md"""
The key variables are the radial coordinate ``r``; ``M_r``, the mass interior to ``r``; ``P``, the local pressure; ``T``, the local temperature; and ``L_r``, the total energy production interior to ``r``. Note that ``P`` and ``T`` are local, while ``M_r`` and ``L_r`` are cumulative. In addition to the key variables just listed, there are numerous auxiliary variables and coefficients used in the equations:

- The radiation constant, ``a = \dfrac{4σ_\text{SB}}{c} = \mathrm{7.565×{10}^{-16} \, J\,m^{-3}\,K^{-4}}``, where ``σ_\text{SB}`` is the Stefan--Boltzmann constant

- The local density, _ρ_:

  ```math
  ρ(P,T) = \mathrm{9.91\!×\!{10}^{-28} \, kg} \;\; \frac{P - aT^4/3}{k_\text{B}T},
  ```

- The mass absorption coefficient, _κ_:

  ```math
  κ(ρ, T) = \mathrm{0.035 \, \frac{m^2}{kg}}
    + \mathrm{6.44\!×\!{10}^{18} \, \frac{m^2}{kg}} \left(\frac{ρ}{1\,\mathrm{kg/m^3}}\right) \left(\frac{T}{1\,\mathrm{K}}\right)^{-3.5},
  ```

- The specific luminosity,
  ```math
   ε(ρ, T) = \mathrm{0.136 \, W\,kg^{-1}} \left(\tilde{ε}_0 + 3.49\!×\!{10}^{12} \tilde{ε}_1 \right),
  ```

  where, in turn,

  ```math
  \begin{align*}
  \tilde{ε}_0 &= \left(\frac{ρ}{\mathrm{1\,kg/m^3}}\right) T_6^{-2/3} e^{-33.80 T_6^{-1/3}}, &
  \tilde{ε}_1 &= \left(\frac{ρ}{\mathrm{1\,kg/m^3}}\right)^2 T_6^{-3} e^{-4403 T_6^{-1}}.
  \end{align*}
  ```
"""

# ╔═╡ fd286a3b-b45c-4891-a1b3-7cf0ca6209fe
"""
    massabsorpcoeff(ρ, T)
Mass absorption coefficient
"""
function massabsorpcoeff(ρ, T)
    scale_factor = ustrip(kg/m^3, ρ)/ustrip(K, T)^3.5
    return (0.035 + 6.44e18 * scale_factor) * m^2/kg
end;

# ╔═╡ a6728349-7323-4e11-a77b-3b05e61992a2
"""``ε(ρ, T)``"""
function specificpower(ρ, T)
    ρ = ustrip(kg/m^3, ρ)
    T₆ = T/1e6K |> NoUnits
    ε̃₀ = ρ / T₆^(2//3) * exp(-33.80/∛T₆)
    ε̃₁ = ρ^2 / T₆^3 * exp(-4403/T₆)
    return 0.136W/kg * (ε̃₀ + 3.49e12 * ε̃₁)
end;

# ╔═╡ d61d868e-7edd-46fc-9994-9c1af60b4cc9
md"""
The above formulae all use SI units, _T_₆ ≡ _T_ / (10⁶ K), and the constants ``c`` and ``k_\text{B}`` have their usual values. These have been calculated assuming a primordial star that is 75% hydrogen and 25% helium by mass, and which is fully ionized throughout. This problem is another boundary value problem. Take a primordial Population-III (i.e. pure H/He) star whose mass is ``100 M_⊙``. At its surface you have ``M_r = M_* = 100 M_⊙``, and ``P = 0``. At its core you have ``L_r = r = 0``.
"""

# ╔═╡ a4b9c073-39de-4097-8b30-b0b14617e843
md"""
Boundary conditions/values (known and unknown):

| Quantity | Inner boundary, ``M_r = 0`` | Outer boundary, ``M_r = M_⋆`` |
|:---------|:----------------------------|:------------------------------|
| ``r``    | ``0``                       | ``R_⋆``                       |
| ``P``    | ``P_c = ?``                 | ``0``                         |
| ``L_r``  | ``0``                       | ``L_⋆``                       |
| ``T``    | ``T_c = ?``                 |                               |
``P_c``, ``T_c``, ``R_⋆``, and ``L_⋆`` are unknown.
"""

# ╔═╡ 0b78405a-ad60-44b8-a2b4-3fc57d866540
function bca(u_a, p)
    r, T, P, Lᵣ = u_a
    # r, T, P, Lᵣ = r*m, T*K, P*Pa, Lᵣ*W

    return [r, Lᵣ]
    #return [ustrip(m, r), ustrip(W, Lᵣ)]
end

# ╔═╡ 54cf2228-0041-476c-b56d-055b39afc596
function bcb(u_b, p)
    r, T, P, Lᵣ = u_b
    # r, T, P, Lᵣ = r*m, T*K, P*Pa, Lᵣ*W

    return [P]
    #return [ustrip(Pa, P)]
end

# ╔═╡ 04d0caf4-f307-47e5-9c2d-7b97da6e8f80
md"""
## Part a
"""

# ╔═╡ 5885df68-9573-4651-9a05-5d0a5b18d747
md"""
Use the shooting method to estimate the core temperature and pressure of this star, to two significant figures each. I suggest you use RK4 as your ODE solver, with 100 steps between ``M_r = 0`` and ``M_r = M_*``. Your goal is to get ``P`` close to 0 at the surface without going negative (since pressure is a positive quantity). To guide your initial guess, the Sun's core temperature is 16 × 10⁶ K, and its core pressure is 1.2 × 10¹⁶ Pa.
"""

# ╔═╡ 92fd50e9-dfdf-447c-ac60-fb46d0080dbc
# u₀_guess = ArrayPartition([0.1m], [16e6K], [1.2e16Pa], [0.1W])
u₀_guess = [1e-9, 16e6, 1.2e16, 1e-10]

# ╔═╡ 13b416d0-39ef-40bb-8997-8e1b7226e10a
solver = Shooting(DP5(), jac_alg = AutoForwardDiff())

# ╔═╡ 9d3e56a2-1eb3-46ad-baf9-fd4bb1f3e8de
# odefun()

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
Plot ``r(M)``, ``T(M)``, and ``L(M)`` from ``0`` to ``M_r`` (It will probably be easier to read if you use three different plots rather than one panel for all three curves.) Use logarithmic y axes and a linear x axis.
"""

# ╔═╡ 45caddb4-8dad-451a-a6dc-fcda9e087dfe
md"""
## Constants reference
"""

# ╔═╡ 7066ef95-8257-4950-a0a7-c15832e1725e
TableOfContents()

# ╔═╡ cd9324f0-87be-4e56-b440-955f1fae11ea
const a = 4σ/c |> (J*m^-3*K^-4)

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

# ╔═╡ 366f153f-e666-45f1-b594-28d0b8381830
function odefun(u, p, Mᵣ)
    r, T, P, Lᵣ = u
    r, T, P, Lᵣ = r*m, T*K, P*Pa, Lᵣ*W
    Mᵣ = Mᵣ * kg

    ρ = density(P, T)
    κ = massabsorpcoeff(ρ, T)

    r′ = 1 / (4π*r^2*ρ) |> m/kg
    T′ = -3/(64*π^2*a*c) * (κ*Lᵣ) / (T^3*r^4) |> K/kg
    P′ = - G*Mᵣ / (4π*r^4) |> Pa/kg
    Lᵣ′ = specificpower(ρ, T) |> W/kg

    # u′ = ArrayPartition([r′], [T′], [P′], [Lᵣ′])
    u′ = [ustrip(m/kg, r′), ustrip(K/kg, T′), ustrip(Pa/kg, P′), ustrip(W/kg, Lᵣ′)]

    return u′
end

# ╔═╡ 6ac5b3da-c540-48ea-95b3-b51719ec36b1
# bvpfunction = BVPFunction(odefun, stellarboundary)
bvpfunction = TwoPointBVPFunction(odefun, (bca, bcb))

# ╔═╡ 2267ba11-3424-48a5-b219-07ee3b2bc4d9
const M_star = 100Msun

# ╔═╡ e74e279d-8f51-4c3b-88a6-efffe6a9707e
# Mᵣ_domain = (0kg, M_star |> kg)
Mᵣ_domain = (0, ustrip(kg, M_star))

# ╔═╡ d8af8e17-d953-41ec-ad1e-ed7c246e23ea
bvproblem = BVProblem(bvpfunction, u₀_guess, Mᵣ_domain)
# bvproblem = TwoPointBVProblem(odefun, (bca, bcb), u₀_guess, Mᵣ_domain)

# ╔═╡ 79c2f0f9-5cb4-4a19-ac97-14d174825c50
solution = solve(bvproblem, solver)

# ╔═╡ e4fb8237-d81a-44c2-85f8-b44a0a11bd3c
function stellarboundary(u, _, Mᵣ)
    u₀ = u(0)
    return [u₀[1], u(ustrip(kg, M_star))[3], u₀[4]]
    return ArrayPartition(
        [u₀[1]],             # at the center, r = 0
        [u(ustrip(kg, M_star))[3]],      # at the surface, P = 0
        [u₀[4]],             # at the center, Lᵣ = 0
    )
end

# ╔═╡ 6e5aa9ff-6879-4d04-9ce5-f47afa13e6e0
function stellarboundary!(residual, u, p, Mᵣ)
    u₀ = u(0)
    residual[1] = u₀[1]         # at the center, r = 0
    residual[2] = u(ustrip(kg, M_star))[3]  # at the surface, P = 0
    residual[3] = u₀[4]         # at the center, Lᵣ = 0
    return residual
end

# ╔═╡ Cell order:
# ╟─da3e18e3-0cda-440b-932e-7c09988e28dc
# ╟─03118617-5cca-4400-9e87-1562d42d7d44
# ╠═366f153f-e666-45f1-b594-28d0b8381830
# ╟─9ed8de43-cd5a-4e40-a8bb-b33038d80a8d
# ╠═9e69b9d6-08e1-4a83-a5ea-102620bdd674
# ╠═fd286a3b-b45c-4891-a1b3-7cf0ca6209fe
# ╠═a6728349-7323-4e11-a77b-3b05e61992a2
# ╟─d61d868e-7edd-46fc-9994-9c1af60b4cc9
# ╠═e74e279d-8f51-4c3b-88a6-efffe6a9707e
# ╟─a4b9c073-39de-4097-8b30-b0b14617e843
# ╠═0b78405a-ad60-44b8-a2b4-3fc57d866540
# ╠═54cf2228-0041-476c-b56d-055b39afc596
# ╟─04d0caf4-f307-47e5-9c2d-7b97da6e8f80
# ╟─5885df68-9573-4651-9a05-5d0a5b18d747
# ╠═e4fb8237-d81a-44c2-85f8-b44a0a11bd3c
# ╠═6e5aa9ff-6879-4d04-9ce5-f47afa13e6e0
# ╠═92fd50e9-dfdf-447c-ac60-fb46d0080dbc
# ╠═6ac5b3da-c540-48ea-95b3-b51719ec36b1
# ╠═d8af8e17-d953-41ec-ad1e-ed7c246e23ea
# ╠═13b416d0-39ef-40bb-8997-8e1b7226e10a
# ╠═9d3e56a2-1eb3-46ad-baf9-fd4bb1f3e8de
# ╠═79c2f0f9-5cb4-4a19-ac97-14d174825c50
# ╟─0d32114c-f8a7-4a09-8a8a-23e18f89b259
# ╟─ce125c69-bc40-42d0-9ded-937257c9bc39
# ╟─d5d6964c-0f6c-4d20-b97a-c6b8aab857d7
# ╟─d7db0d1a-92fc-4a43-8285-7438f8d6a8e1
# ╠═ae5820a2-26c3-4bf7-9dd5-d136eaf35038
# ╟─45caddb4-8dad-451a-a6dc-fcda9e087dfe
# ╠═93698b5a-1f3b-4c47-9e67-aff44f234c0f
# ╠═7066ef95-8257-4950-a0a7-c15832e1725e
# ╠═1e6a2dc7-ff39-4463-a6d3-7a87ecc6ae53
# ╠═44d9b444-3ed6-4df0-866b-7d683be2d278
# ╠═81b35b32-2be1-4458-bf91-97780765a6fd
# ╠═5550c36f-dcef-4d95-9c90-46e91d276d37
# ╠═54f413fa-3a90-4e76-a192-75afa6a5bd6d
# ╠═67dead86-da4e-482e-b87e-f87e51ff5e5b
# ╠═329fa343-39e2-4a0b-adef-9c590a18f19c
# ╠═cd9324f0-87be-4e56-b440-955f1fae11ea
# ╠═58810d44-0b8f-4f74-83b9-9b95605cd2e7
# ╠═b0a76914-5907-49ad-9b7c-eb3a3e4614e7
# ╠═ba8a7b26-cba3-411d-ba6c-7e8fcef66c58
# ╠═b44185c8-9325-4aa9-a448-51ff08d95187
# ╟─c68905f0-23ce-498d-a478-e6e41c19308b
# ╠═817a6fc5-b327-48ea-9cda-4a4217931cc9
# ╠═2267ba11-3424-48a5-b219-07ee3b2bc4d9
