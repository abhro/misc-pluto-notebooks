### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ d44f2cfb-1db0-4f96-9baa-d0a17958fc7a
using DifferentialEquations

# ╔═╡ 37d900fd-17be-4a8f-a7bc-94bd3ed65f43
using LinearAlgebra: norm

# ╔═╡ 5e94d6c6-88e7-444e-9f45-cf3ef46b3c16
using DiffEqNoiseProcess

# ╔═╡ 77815dba-ede8-11ee-393f-f77d44dc5d35
md"""
Differential equation for each of the phase space components
"""

# ╔═╡ 2ee2907f-daed-4e91-88ad-f9ccab853980
md"""
```math
\frac{∂f}{∂t}
=
\boldsymbol{∇} ⋅ \boldsymbol{κ}_⟂ \cdot \boldsymbol{∇} f
- \left(vμ \hat{\mathbf{b}} + \mathbf{V} + \mathbf{V}_d\right) \cdot \boldsymbol{∇} f
+ \frac{∂}{∂μ} D_{μμ} \frac{∂f}{∂μ}
- \frac{dμ}{dt} \frac{∂f}{∂μ}
- \frac{dp}{dt} \frac{∂f}{∂p}.
```
"""

# ╔═╡ 122c5bff-e8d1-44a2-befc-db783073d82c
md"""
Can be recast as 5 stochastic differential equations
```math
\begin{align}
d\mathbf{x}(t) &= \sqrt{2 \boldsymbol{κ}_⟂} ⋅ d\mathbf{w}(t) + \left(\boldsymbol{∇} ⋅ \boldsymbol{κ}_⟂ + vμ\hat{\mathbf{b}} + \mathbf{V} + \mathbf{V}_d\right) dt \\
dμ(t) &= \sqrt{2 \max(D_{μμ},0)} \, dw(t) + \left(\frac{∂D_{μμ}}{∂μ} + \frac{dμ}{dt}\right) dt \\
dp(t) &= \frac{dp}{dt} \, dt
\end{align}
```
"""

# ╔═╡ 4db02831-a76a-49be-9d19-3f449597f5f7
md"""
``w(t)`` is a normally distributed Wiener noise process.
"""

# ╔═╡ 3732a5b2-5810-48d5-953d-7ae4c906ff47
md"""
Perpendicular diffusion coefficient:

```math
κ_⟂ = \frac{v}{2} D_\text{FLRW}
```

Assume ``κ_⟂`` is a diagonal tensor.

```math
\boldsymbol{κ}_⟂ = \begin{bmatrix} κ_⟂ & ⋅ & ⋅ \\ ⋅ & κ_⟂ & ⋅ \\ ⋅ & ⋅ & κ_⟂ \end{bmatrix}
```

Identity: ``κ_⟂ = \frac{1}{3} \operatorname{tr}(\boldsymbol{κ}_⟂)``

Divergence of the tensor is the gradient of the scalar:

```math
\boldsymbol{∇} κ_⟂ = \boldsymbol{∇} ⋅ \boldsymbol{κ}_⟂
= \tfrac{1}{2} \boldsymbol{∇} (v D_\text{FLRW})
= \tfrac{1}{2} \left(\boldsymbol{∇} v ⋅ D_\text{FLRW} + v ⋅ \boldsymbol{∇} D_\text{FLRW} \right)
```

Question: can the square root of the tensor be taken as the square root of the scalar?
```math
\sqrt{2 \boldsymbol{κ}_⟂} ⋅ d\mathbf{w}(t) ~\overset{?}{=}~ \sqrt{2 κ_⟂} \, d\mathbf{w}(t)
```

Answer: Yes? If it's a multiple of the identity tensor?

```math
\sqrt{λ \mathbf{1}} = \sqrt{λ} \sqrt{\mathbf{1}} = \sqrt{λ} \mathbf{1}
```
"""

# ╔═╡ 177c7d6e-a0ed-42c5-bc5d-31504bbfee97
md"""
Model the equations as a 5d system:

```math
\begin{align*}
\mathbf{u}(t) &= \begin{bmatrix} \mathbf{x}(t) \\ μ(t) \\ p(t) \end{bmatrix} \\[2ex]
d\mathbf{u} &= f\,dt + g\,dw \\[1.5ex]
d\mathbf{u} &= \begin{bmatrix}
    \boldsymbol{∇} ⋅ \boldsymbol{κ}_⟂ + vμ\hat{\mathbf{b}} + \mathbf{V} + \mathbf{V}_d \\[0.5ex]
    \displaystyle \frac{∂D_{μμ}}{∂μ} + \frac{dμ}{dt} \\[0.5ex]
    \displaystyle \frac{dp}{dt} \end{bmatrix} dt + \begin{bmatrix} \sqrt{2κ_⟂} \\[3.3ex]
    \sqrt{2\max(D_{μμ},0)} \\[3.3ex]
    0
\end{bmatrix} dw
\end{align*}
```
"""

# ╔═╡ 7f35e344-f19f-4c1a-bf95-63bd735395b4
function f!(du, u, p, t)
    (r, θ, φ, μ, p) = u
    b = magfield(r, θ, φ)
    du[1:3] .= grad(κ_⊥) + norm(v) * μ * b/norm(b) + 0 + V_d
    du[4] = ∂D_μμ∂t + dμdt
    du[5] = dpdt
    return du
end

# ╔═╡ 3bf7487a-7c84-443b-9cf9-9e99ea4d4a90
function g!(du, u, p, t)
    du[1:3] .= √(2κ_⊥)
    du[4] = √(2 * max(D_μμ, 0))
    du[5] = 0
    return du
end

# ╔═╡ bbf5d07b-9d78-42bb-8562-a7d342107f4e
# TEMPORARY
tspan = (0, 1)

# ╔═╡ 9af9b0ee-b588-41ec-aac6-2f68fa5c35b5
# TEMPORARY
u_0 = zeros(5)

# ╔═╡ 01bd5f8a-4c79-4fe3-9575-e8a1ef578190
prob = SDEProblem(f!, g!, u_0, tspan)

# ╔═╡ Cell order:
# ╟─77815dba-ede8-11ee-393f-f77d44dc5d35
# ╟─2ee2907f-daed-4e91-88ad-f9ccab853980
# ╟─122c5bff-e8d1-44a2-befc-db783073d82c
# ╟─4db02831-a76a-49be-9d19-3f449597f5f7
# ╟─3732a5b2-5810-48d5-953d-7ae4c906ff47
# ╠═d44f2cfb-1db0-4f96-9baa-d0a17958fc7a
# ╠═37d900fd-17be-4a8f-a7bc-94bd3ed65f43
# ╠═5e94d6c6-88e7-444e-9f45-cf3ef46b3c16
# ╟─177c7d6e-a0ed-42c5-bc5d-31504bbfee97
# ╠═7f35e344-f19f-4c1a-bf95-63bd735395b4
# ╠═3bf7487a-7c84-443b-9cf9-9e99ea4d4a90
# ╠═bbf5d07b-9d78-42bb-8562-a7d342107f4e
# ╠═9af9b0ee-b588-41ec-aac6-2f68fa5c35b5
# ╠═01bd5f8a-4c79-4fe3-9575-e8a1ef578190
