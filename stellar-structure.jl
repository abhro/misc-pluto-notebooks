### A Pluto.jl notebook ###
# v0.20.9

using Markdown
using InteractiveUtils

# ╔═╡ 93698b5a-1f3b-4c47-9e67-aff44f234c0f
using PlutoUI

# ╔═╡ 67dead86-da4e-482e-b87e-f87e51ff5e5b
begin
	using Unitful, UnitfulAstro
	using Unitful: m, kg, J, K, W
	using UnitfulAstro: Msun
end

# ╔═╡ 329fa343-39e2-4a0b-adef-9c590a18f19c
using PhysicalConstants.CODATA2018: σ, k_B, c_0 as c, G

# ╔═╡ ae5820a2-26c3-4bf7-9dd5-d136eaf35038
using CairoMakie

# ╔═╡ 7066ef95-8257-4950-a0a7-c15832e1725e
TableOfContents()

# ╔═╡ da3e18e3-0cda-440b-932e-7c09988e28dc
md"""
# Equations of Stellar Structure
"""

# ╔═╡ 03118617-5cca-4400-9e87-1562d42d7d44
md"""
For a non-rotating star in hydrostatic equilibrium (i.e. no time dependent behavior), the four equations of stellar structure are as follows:

```math
\begin{align*}
\frac{dr}{dM_r} &= \frac{1}{4πr^2 ρ} &
\frac{dT}{dM_r} &= -\frac{3 κ L_r}{64π^2 ac T^3 r^4} \\[1em]
\frac{dP}{dM_r} &= -\frac{GM_r}{4πr^4} &
\frac{dL_r}{dM_r} &= ε
\end{align*}
```

The key variables are the radial coordinate ``r``; ``M_r``, the mass interior to ``r``; ``P``, the local pressure; ``T``, the local temperature; and ``L_r``, the total energy production interior to ``r``. Note that ``P`` and ``T`` are local, while ``M_r`` and ``L_r`` are cumulative. In addition to the key variables just listed, there are numerous auxiliary variables and coefficients used in the equations:

```math
\begin{align*}
ρ &= 9.91 \times {10}^{-28} \frac{P - aT^4/3}{k_\text{B}T} ~~ \mathrm{[kg\,m^{-3}]} \\
κ &= 0.035 + 6.44 \times {10}^{18} \frac{ρ}{T^{3.5}} ~~ \mathrm{[m^2\,kg^{-1}]} \\
a &= 7.565 \times {10}^{-16} ~~ \mathrm{[J\,m^{-3}\,K^{-4}]} \\
ε &= 0.136 ρ T_6^{-2/3} e^{-33.80 T_6^{-1/3}} + 4.69 \times {10}^{11} ρ^2 T_6^{-3} e^{-4403 T_6^{-1}} ~~\mathrm{[W\,kg^{-1}]}
\end{align*}
```

The above formulae all use SI units, T₆ ≡ T / (10⁶ K), and the constants ``c`` and ``k_\text{B}`` have their usual values. These have been calculated assuming a primordial star that is 75% hydrogen and 25% helium by mass, and which is fully ionized throughout. This problem is another boundary value problem. Take a primordial Population-III (i.e. pure H/He) star whose mass is ``100 M_⊙``. At its surface you have ``M_r = M_* = 100 M_⊙``, and ``P = 0``. At its core you have ``L_r = r = 0``.
"""

# ╔═╡ cd9324f0-87be-4e56-b440-955f1fae11ea
const a = 4σ/c |> (J*m^-3*K^-4)

# ╔═╡ 9e69b9d6-08e1-4a83-a5ea-102620bdd674
ρfunc(P,T) = 9.91e-28kg/m^3 * NoUnits((P - a*T^4/3) / (k_B*T))

# ╔═╡ fd286a3b-b45c-4891-a1b3-7cf0ca6209fe
κfunc(ρ,T) = (0.035 + 6.44e18 * ρ/ustrip(K, T^3.5)) * m^2*kg^-1

# ╔═╡ 366f153f-e666-45f1-b594-28d0b8381830
function odefun(u, p, Mᵣ)
	r, T, P, Lᵣ = u

	ρ = ρfunc(P, T)
	κ = κfunc(ρ, T)

	dr = 1 / (4π*r^2*ρ)
	dT = -3/(64*π^2*a*c) * (κ*Lᵣ) / (T^3*r^4)
	dP = - G*Mᵣ / (4π*r^4)
	dLᵣ = ε(ρ, T/1e6K)

	du = [dr, dT, dP, dLᵣ]

	return du
end

# ╔═╡ a6728349-7323-4e11-a77b-3b05e61992a2
εfunc(ρ,T₆) = (
	(0.136*ρ/T₆^(2//3) * exp(-33.80/∛T₆)
	+ 4.69e11*ρ^2/T₆^3 * exp(-4403/T₆))W/kg)

# ╔═╡ e74e279d-8f51-4c3b-88a6-efffe6a9707e
Mᵣ_domain = (0kg, 100Msun |> kg)

# ╔═╡ 04d0caf4-f307-47e5-9c2d-7b97da6e8f80
md"""
## Part a
"""

# ╔═╡ 5885df68-9573-4651-9a05-5d0a5b18d747
md"""
Use the shooting method to estimate the core temperature and pressure of this star, to two significant figures each. I suggest you use RK4 as your ODE solver, with 100 steps between ``M_r = 0`` and ``M_r = M_*``. Your goal is to get ``P`` close to 0 at the surface without going negative (since pressure is a positive quantity). To guide your initial guess, the Sun's core temperature is 16 × 10⁶ K, and its core pressure is 1.2 × 10¹⁶ Pa.
"""

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

# ╔═╡ 58810d44-0b8f-4f74-83b9-9b95605cd2e7
σ

# ╔═╡ b0a76914-5907-49ad-9b7c-eb3a3e4614e7
k_B

# ╔═╡ ba8a7b26-cba3-411d-ba6c-7e8fcef66c58
G

# ╔═╡ Cell order:
# ╠═93698b5a-1f3b-4c47-9e67-aff44f234c0f
# ╠═7066ef95-8257-4950-a0a7-c15832e1725e
# ╟─da3e18e3-0cda-440b-932e-7c09988e28dc
# ╟─03118617-5cca-4400-9e87-1562d42d7d44
# ╠═366f153f-e666-45f1-b594-28d0b8381830
# ╠═67dead86-da4e-482e-b87e-f87e51ff5e5b
# ╠═329fa343-39e2-4a0b-adef-9c590a18f19c
# ╠═cd9324f0-87be-4e56-b440-955f1fae11ea
# ╠═9e69b9d6-08e1-4a83-a5ea-102620bdd674
# ╠═fd286a3b-b45c-4891-a1b3-7cf0ca6209fe
# ╠═a6728349-7323-4e11-a77b-3b05e61992a2
# ╠═e74e279d-8f51-4c3b-88a6-efffe6a9707e
# ╟─04d0caf4-f307-47e5-9c2d-7b97da6e8f80
# ╟─5885df68-9573-4651-9a05-5d0a5b18d747
# ╟─0d32114c-f8a7-4a09-8a8a-23e18f89b259
# ╟─ce125c69-bc40-42d0-9ded-937257c9bc39
# ╟─d5d6964c-0f6c-4d20-b97a-c6b8aab857d7
# ╟─d7db0d1a-92fc-4a43-8285-7438f8d6a8e1
# ╠═ae5820a2-26c3-4bf7-9dd5-d136eaf35038
# ╟─45caddb4-8dad-451a-a6dc-fcda9e087dfe
# ╠═58810d44-0b8f-4f74-83b9-9b95605cd2e7
# ╠═b0a76914-5907-49ad-9b7c-eb3a3e4614e7
# ╠═ba8a7b26-cba3-411d-ba6c-7e8fcef66c58
