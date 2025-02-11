### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ c8314c51-ff6d-459a-bdd6-0ff217a2ef47
using PlutoUI

# ╔═╡ b38d823e-0ab3-4e08-afe9-5c9b46e83ec5
using GLMakie

# ╔═╡ 41cb3e44-f842-4e58-84df-420c0ce6416a
using Unitful

# ╔═╡ dc2be4b2-ed93-4bcf-b6be-3aceddcbc39d
using StaticArrays

# ╔═╡ 87d06c27-9d29-489b-9fba-233c6d5d9203
using PlutoUI: Slider

# ╔═╡ d60444c0-ddce-11ef-2e4e-234174458fd9
md"""
# Motion of single charged particle
"""

# ╔═╡ 9b918a1b-40f4-4c4b-859c-b69dce887a7e
TableOfContents()

# ╔═╡ 10c298cf-7d23-46bd-8708-ce92ad2595eb
md"""
General force law:
```math
m \frac{d\mathbf{v}}{dt} = q\mathbf{E} + q\mathbf{v}×\mathbf{B} + \mathbf{F}_\text{ext}
```
"""

# ╔═╡ 7429f193-2695-4dcb-be2e-5a9f49716621
md"""
List of velocity notations:
- ``\mathbf{u}``: Gyromotion
- ``\mathbf{v}_g``: Guiding center motion
- ``\mathbf{v}_E``: **E**×**B** drift
- ``\mathbf{v}_D``: Guiding center drift
"""

# ╔═╡ f909adaa-cd08-47c6-8adb-ef2a3796ffe3
md"""
Plot projection: $(@bind do_proj CheckBox(default = true))
"""

# ╔═╡ 16a381a7-e3ac-42be-b94d-05a1fba351e5
md"""
Plot guiding center: $(@bind do_guiding_center CheckBox(default = true))
"""

# ╔═╡ ee41070b-b68b-4f5f-98b6-c9843f9587eb
md"""
## Constant and uniform **B** field
"""

# ╔═╡ 289316fd-59ee-4c9e-9bad-3dd2a93fcf21
md"""
```math
\begin{align}
\mathbf{E} &= \mathbf{0} \\
\mathbf{B} &= B \hat{\mathbf{z}} \\
\mathbf{F}_\text{ext} &= \mathbf{0}
\end{align}
```
"""

# ╔═╡ 8e214cc4-c04f-4739-a573-b539014c55aa
md"""
Solution:
"""

# ╔═╡ 6743bdc7-8f4f-494f-bdf8-33002a2b68b6
md"""
```math
\begin{align}
\mathbf{x}(t) &= \mathbf{x}_0 + \boldsymbol{ρ} + \boldsymbol{ξ}(t)
= \mathbf{R}_C + \boldsymbol{ξ}(t) \\[1.5ex]
&= \begin{pmatrix}
    x_0 \\
    y_0 \\
    z_0
\end{pmatrix}
+ \begin{pmatrix}
    \hphantom{-} \frac{v_{0⟂}}{ω_c} \cos(γ_0) \\
    - \frac{v_{0⟂}}{ω_c} \sin(γ_0) \\
    0
\end{pmatrix}
+
\begin{pmatrix}
    - \frac{v_{0⟂}}{ω_c} \cos(ω_c t + γ_0) \\
    \hphantom{-}\frac{v_{0⟂}}{ω_c} \sin(ω_c t + γ_0) \\
    v_{0∥}t
\end{pmatrix} \\[2ex]
&= \begin{pmatrix}
    x_0 + \frac{v_{0⟂}}{ω_c} \cos(γ_0) - \frac{v_{0⟂}}{ω_c} \cos(ω_c t + γ_0) \\
    y_0 - \frac{v_{0⟂}}{ω_c} \sin(γ_0) + \frac{v_{0⟂}}{ω_c} \sin(ω_c t + γ_0) \\
    z_0 + v_{0∥}t
\end{pmatrix}
\end{align}
```
"""

# ╔═╡ c4ffcbdc-02dc-46bf-8f7b-ca480ce2a034
md"""
```math
\begin{align*}
\mathbf{v}(t)
&= \mathbf{u}(t) + \mathbf{v}_\parallel(t) \\
&=
\begin{pmatrix}
    v_{0⟂} \sin(ω_c t + γ_0) \\
    v_{0⟂} \cos(ω_c t + γ_0) \\
    0
\end{pmatrix}
+
\begin{pmatrix} 0 \\ 0 \\ v_{0∥} \end{pmatrix} \\
&= \begin{pmatrix}
    v_{0⟂} \sin(ω_c t + γ_0) \\
    v_{0⟂} \cos(ω_c t + γ_0) \\
    v_{0∥}
\end{pmatrix}
\end{align*}
```
"""

# ╔═╡ 866c51f0-ee41-4147-866b-ad1cc386fe1f
md"""
- ``ω_c = \dfrac{qB}{m}``
- ``v_{0⟂} = \sqrt{v_{0,x}^2 + v_{0,y}^2}``
- ``γ_0 = \arctan\left(\dfrac{v_{0,y}}{v_{0,x}}\right)``
"""

# ╔═╡ a8810170-7c78-4e63-b8aa-bc7c3d6b3f04
md"""
Choose system parameters:
"""

# ╔═╡ b2c4a082-0b67-4dee-90c2-849e73edf3cc
md"""
Particle charge = $(@bind q Select([Unitful.q => "+e", -Unitful.q => "-e"]))
"""

# ╔═╡ 771624c1-ade2-4584-a841-3db0132d1be6
md"""
Particle mass = $(@bind m Select([Unitful.mp => "proton mass", Unitful.me => "electron mass"]))
"""

# ╔═╡ 9a791de7-144a-4537-b258-08613a1920a3
md"""
Magnetic field strength:
"""

# ╔═╡ 36146cf5-c26a-4759-839a-13c2adc29d96
B = 3.0u"nT";

# ╔═╡ 73d642db-b214-45c3-a16f-558c18a30011
ω_c = q * B / m |> u"rad/s"

# ╔═╡ f10f18d0-8f85-44be-8b5a-474fc2a90d2a
md"""
Choose initial conditions:
"""

# ╔═╡ 1a8691be-4fa9-4a21-afcc-10d0d6baf716
md"""
Initial position:
"""

# ╔═╡ 62dce489-601c-4d04-a999-068face9b775
x₀ = SVector(0.3, 0.3, 0.3)u"m";

# ╔═╡ c6243506-4a2a-4f37-ad15-a33e4d50d7ea
v_range = range(0.0, 5, step = 0.1)u"m/s"

# ╔═╡ b06f0f86-eaa1-4dbb-9c21-1baa251d4c10
md"""
Initial velocity:
"""

# ╔═╡ 55b22b9e-0727-4ed0-80c5-095175ad7be0
md"""
Perpendicular component: v\_0⟂ = $(
    @bind v_perp Slider(v_range, show_value = true, default = 1u"m/s")
)
"""

# ╔═╡ 98891a54-cee8-4b20-9c62-34aa51557e17
md"""
Parallel component: v\_0∥ = $(
    @bind v_parallel Slider(v_range, show_value = true, default = 1u"m/s")
)
"""

# ╔═╡ ec75d67b-b248-471c-912e-f68096e4133c
md"""
Velocity phase angle: γ₀ = $(@bind γ₀ Slider((0:359)u"°", show_value = true))
"""

# ╔═╡ e3382bd2-d27f-4434-8417-4d91cbb7b9be
v₀ = SVector(v_perp*cos(γ₀), v_perp*sin(γ₀), v_parallel)

# ╔═╡ 96b90227-a92a-45d7-a37b-fd07f3ea0b6c
md"""
Choose animation parameters:
"""

# ╔═╡ 88ada895-2d26-46b4-abfa-4e0dd7187339
md"""
Time span:
"""

# ╔═╡ 178b6fa9-4066-4acb-8cbc-73bae83930a3
tspan = range(0.0u"s", 300.0u"s", length = 2000);

# ╔═╡ 84d3fd17-d184-4051-a932-0b1784fc2cba
md"""
TODO: Choose bounds from each 3 sets of points (particle, gyrocenter, projection)
"""

# ╔═╡ 2509dad4-0058-42f5-8a24-1dbec0b6772b
function bounds(points)
    d = length(eltype(points))
    T = eltype(eltype(points))
    bounds = Vector{Tuple{T,T}}(undef, d)
    for i in 1:d
        bounds[i] = extrema(x[i] for x in points)
    end
    bounds
end

# ╔═╡ 756c7bb1-fb78-4b21-83b2-ccb6b71d1ba3
function track_motion(fig, ax, points, t = eachindex(points))
    traj = Observable(points[begin:begin+1])
    lead_point = Observable(points[begin+1])

    p = lines!(ax, traj)
    scatter!(ax, lead_point)

    r = Record(fig) do io
        for (t_i, x_i) in zip(t, points)
            push!(traj[], x_i) # add new point to trajectory
            lead_point[] = x_i # update leading point
            traj[] = traj[]
            recordframe!(io)
        end
    end
    return (p, r)
end

# ╔═╡ e36633d5-648c-41fc-a553-7b373e6be6f2
function track_motion(fig, ax, p::NTuple{2}, t = eachindex(p[1]))
    points1, points2 = p
    traj1 = Observable(points1[begin:begin+1])
    lead_point1 = Observable(points1[begin+1])
    traj2 = Observable(points2[begin:begin+1])
    lead_point2 = Observable(points2[begin+1])

    p1 = lines!(ax, traj1)
    scatter!(ax, lead_point1)
    p2 = lines!(ax, traj2)
    scatter!(ax, lead_point2)

    r = Record(fig) do io
        for (t_i, x1, x2) in zip(t, points1, points2)
            push!(traj1[], x1) # add new point to trajectory
            push!(traj2[], x2) # add new point to trajectory
            lead_point1[] = x1 # update leading point
            lead_point2[] = x2 # update leading point
            traj1[] = traj1[]
            traj2[] = traj2[]
            recordframe!(io)
        end
    end
    return (p1, p2, r)
end

# ╔═╡ 7fc0454a-25a3-4a7d-9b67-30b3679c4fb0
function track_motion(fig, ax, p::NTuple{3}, t = eachindex(p[1]))
    points1, points2, points3 = p
    traj1 = Observable(points1[begin:begin+1])
    lead_point1 = Observable(points1[begin+1])
    traj2 = Observable(points2[begin:begin+1])
    lead_point2 = Observable(points2[begin+1])
    traj3 = Observable(points3[begin:begin+1])
    lead_point3 = Observable(points3[begin+1])

    p1 = lines!(ax, traj1)
    scatter!(ax, lead_point1)
    p2 = lines!(ax, traj2)
    scatter!(ax, lead_point2)
    p3 = lines!(ax, traj3)
    scatter!(ax, lead_point3)

    r = Record(fig) do io
        for (t_i, x1, x2, x3) in zip(t, points1, points2, points3)
            push!(traj1[], x1) # add new point to trajectory
            push!(traj2[], x2) # add new point to trajectory
            push!(traj3[], x3) # add new point to trajectory
            lead_point1[] = x1 # update leading point
            lead_point2[] = x2 # update leading point
            lead_point3[] = x3 # update leading point
            traj1[] = traj1[]
            traj2[] = traj2[]
            traj3[] = traj3[]
            recordframe!(io)
        end
    end
    return (p1, p2, p3, r)
end

# ╔═╡ 6d736b0b-9e26-4800-8925-f287fc5d292b
function gyromotion(time, params)
    v_perp = params.v_perp
    v_parallel = params.v_parallel
    ω_c = params.ω_c
    γ₀ = params.γ₀
    x₀ = params.x₀
    return SVector(
        -v_perp/ω_c * cos(ω_c * time + γ₀),
        -v_perp/ω_c * sin(ω_c * time + γ₀),
        zero(params.x₀.z))
end

# ╔═╡ fdc85a17-c829-4f0e-8f4e-9fb5909aaea3
x = gyromotion.(tspan, Ref((; ω_c, γ₀, v_perp, v_parallel, x₀, v₀)))

# ╔═╡ e120ae58-3593-45e2-b884-f940bdca716d
bounds(x)

# ╔═╡ c4914c73-bbf5-4dab-861a-db045f0fc5dd
let
    f = Figure()
    ax = Axis3(f[1,1])

    #time = Observable(first(tspan))

    #circ = @lift(gyromotion($time, (; ω_c, γ₀), x₀, v₀))
    #scatter!(ax, Point3(circ))

    gc = SVector.(0u"m", 0u"m", v_perp .* tspan)
    ξ = gyromotion.(tspan, Ref((; ω_c, γ₀, v_perp, v_parallel, x₀, v₀)))
    x = Vector{Point3f}(undef, length(tspan))
    for i in eachindex(x)
        x[i] = Point3f(ustrip.(u"m", gc[i] + ξ[i]))
    end
    xall = x
    gc = let
        gc_stripped = Vector{Point3f}(undef, length(gc))
        for i in eachindex(gc)
            gc_stripped[i] = Point3f(ustrip.(u"m", gc[i]))
        end
        gc_stripped
    end
    #xall = x.(tspan) + ustrip.(u"m", gc)

    projection = Vector{Point3f}(undef, length(x))
    for i in eachindex(projection)
        projection[i] = Point3f(x[i][1], x[i][2], 0)
    end

    #traj = Observable(xall[1:2])
    #lead_point = Observable(xall[2])

    b = bounds(xall)
    limits!(ax, b...)

    #lines!(ax, traj)
    #scatter!(ax, lead_point)

    r = if do_guiding_center
        if do_proj
            p1, p2, p3, r = track_motion(f, ax, (xall, projection, gc))
            axislegend(ax, [p1, p2, p3], ["Particle", "Projection", "Gyrocenter"])
            r
        else
            _, _, r = track_motion(f, ax, (xall, gc))
            r
        end
    else
        if do_proj
            _, _, r = track_motion(f, ax, (xall, projection))
            r
        else
            _, r = track_motion(f, ax, xall)
            r
        end
    end
    r
end

# ╔═╡ 71069a27-0467-4352-9aff-20db6a27a663
md"""
## Constant and uniform **B** field and **F**
"""

# ╔═╡ 6c2aebe1-9893-42b2-850e-2a1bed21180d
md"""
```math
\begin{align}
\mathbf{E} &= \mathbf{0} \\
\mathbf{B} &= B \hat{\mathbf{z}} \\
\mathbf{F}_\text{ext} &= \text{const.}
\end{align}
```
"""

# ╔═╡ 57ff82cd-5b15-4121-aa7c-fd37e5e712a5
md"""
Decompose ``\mathbf{F}`` into:
```math
\begin{align*}
\mathbf{F}_⟂ &= F_x \hat{\mathbf{x}} + F_y \hat{\mathbf{y}} \\
\mathbf{F}_∥ &= F_z \hat{\mathbf{z}}
\end{align*}
```
"""

# ╔═╡ c8f06a99-2adb-406c-aa32-c3347616e480
md"""
Solution:
"""

# ╔═╡ cd2517b7-8507-45f3-8d52-ce8e33818bae
md"""
```math
\mathbf{x}(t) = \begin{pmatrix}
    \\
    \\
    z_0 + v_{\parallel,0}t + \dfrac{F_\parallel}{2m} t^2
\end{pmatrix}
```
"""

# ╔═╡ bde4ef61-4c56-440f-9082-8decd9b0137a
md"""
```math
\begin{align}
\mathbf{v}(t)
&= \mathbf{u} + \mathbf{v}_\text{D} + \left(\frac{F_∥}{m} t + v_{∥0}\right) \hat{\mathbf{z}} \\
&= \begin{pmatrix}
\hphantom{-} \frac{F_y}{qB} \\
-\frac{F_x}{qB} \\
0
\end{pmatrix}
+
\begin{pmatrix}
0 \\
0 \\
\frac{F_∥}{m} t + v_{∥0}
\end{pmatrix}
\end{align}
```
"""

# ╔═╡ 8f9070d9-9412-4ccd-b84a-d8a20b6f7059
md"""
## Constant and uniform **E** and **B** fields
"""

# ╔═╡ 629f429f-8b33-4f25-af62-823c89da67a5
md"""
```math
\begin{align}
\mathbf{E} &= \text{const.} \\
\mathbf{B} &= B \hat{\mathbf{z}} \\
\mathbf{F}_\text{ext} &= 0
\end{align}
```
"""

# ╔═╡ cb013904-1a68-4655-848b-4d520e66fc48
md"""
Solution:
"""

# ╔═╡ b345fb68-3b7f-4ec2-aa85-77f3f5b94fa1
md"""
```math
\mathbf{x}(t) = \begin{pmatrix}
    \\
    \\
    z_0 + v_{\parallel,0}t + \dfrac{E_\parallel}{2m} t^2
\end{pmatrix}
```
"""

# ╔═╡ a56195b7-3170-4c59-99d1-1ea8bae7bf69
md"""
```math
\begin{align}
\mathbf{v}(t)
&= \mathbf{u} + \mathbf{v}_E + \left(\frac{qE_∥}{m} t + v_{∥0}\right) \hat{\mathbf{z}} \\
&= \begin{pmatrix}
\\
\\

\end{pmatrix}
\end{align}
```
"""

# ╔═╡ bd6f39b4-a999-44d2-a8b0-646c1e0a26a3
md"""
## Constant and uniform **F**, **E**, and **B**
"""

# ╔═╡ d9f14d74-164f-42ab-aac3-4d52a1b8917b
function guidingcenter(time, params)
    # TODO generalize for non-zero E, Fₑₓₜ
    v_perp = hypot(v₀.x, v₀.y)
    ω_c = params.ω_c
    γ₀ = params.γ₀
    return x₀ + v_perp/ω_c * SVector(cos(γ₀), -sin(γ₀), 0)
end


# ╔═╡ Cell order:
# ╟─d60444c0-ddce-11ef-2e4e-234174458fd9
# ╠═c8314c51-ff6d-459a-bdd6-0ff217a2ef47
# ╠═9b918a1b-40f4-4c4b-859c-b69dce887a7e
# ╠═b38d823e-0ab3-4e08-afe9-5c9b46e83ec5
# ╠═41cb3e44-f842-4e58-84df-420c0ce6416a
# ╠═dc2be4b2-ed93-4bcf-b6be-3aceddcbc39d
# ╠═87d06c27-9d29-489b-9fba-233c6d5d9203
# ╟─10c298cf-7d23-46bd-8708-ce92ad2595eb
# ╟─7429f193-2695-4dcb-be2e-5a9f49716621
# ╟─f909adaa-cd08-47c6-8adb-ef2a3796ffe3
# ╟─16a381a7-e3ac-42be-b94d-05a1fba351e5
# ╟─ee41070b-b68b-4f5f-98b6-c9843f9587eb
# ╟─289316fd-59ee-4c9e-9bad-3dd2a93fcf21
# ╟─8e214cc4-c04f-4739-a573-b539014c55aa
# ╟─6743bdc7-8f4f-494f-bdf8-33002a2b68b6
# ╟─c4ffcbdc-02dc-46bf-8f7b-ca480ce2a034
# ╟─866c51f0-ee41-4147-866b-ad1cc386fe1f
# ╟─a8810170-7c78-4e63-b8aa-bc7c3d6b3f04
# ╟─b2c4a082-0b67-4dee-90c2-849e73edf3cc
# ╟─771624c1-ade2-4584-a841-3db0132d1be6
# ╟─9a791de7-144a-4537-b258-08613a1920a3
# ╠═36146cf5-c26a-4759-839a-13c2adc29d96
# ╟─73d642db-b214-45c3-a16f-558c18a30011
# ╟─f10f18d0-8f85-44be-8b5a-474fc2a90d2a
# ╟─1a8691be-4fa9-4a21-afcc-10d0d6baf716
# ╠═62dce489-601c-4d04-a999-068face9b775
# ╠═c6243506-4a2a-4f37-ad15-a33e4d50d7ea
# ╟─b06f0f86-eaa1-4dbb-9c21-1baa251d4c10
# ╟─55b22b9e-0727-4ed0-80c5-095175ad7be0
# ╟─98891a54-cee8-4b20-9c62-34aa51557e17
# ╟─ec75d67b-b248-471c-912e-f68096e4133c
# ╠═e3382bd2-d27f-4434-8417-4d91cbb7b9be
# ╟─96b90227-a92a-45d7-a37b-fd07f3ea0b6c
# ╟─88ada895-2d26-46b4-abfa-4e0dd7187339
# ╠═178b6fa9-4066-4acb-8cbc-73bae83930a3
# ╠═fdc85a17-c829-4f0e-8f4e-9fb5909aaea3
# ╟─84d3fd17-d184-4051-a932-0b1784fc2cba
# ╠═2509dad4-0058-42f5-8a24-1dbec0b6772b
# ╠═e120ae58-3593-45e2-b884-f940bdca716d
# ╠═c4914c73-bbf5-4dab-861a-db045f0fc5dd
# ╠═756c7bb1-fb78-4b21-83b2-ccb6b71d1ba3
# ╠═e36633d5-648c-41fc-a553-7b373e6be6f2
# ╠═7fc0454a-25a3-4a7d-9b67-30b3679c4fb0
# ╠═6d736b0b-9e26-4800-8925-f287fc5d292b
# ╟─71069a27-0467-4352-9aff-20db6a27a663
# ╟─6c2aebe1-9893-42b2-850e-2a1bed21180d
# ╟─57ff82cd-5b15-4121-aa7c-fd37e5e712a5
# ╟─c8f06a99-2adb-406c-aa32-c3347616e480
# ╠═cd2517b7-8507-45f3-8d52-ce8e33818bae
# ╠═bde4ef61-4c56-440f-9082-8decd9b0137a
# ╟─8f9070d9-9412-4ccd-b84a-d8a20b6f7059
# ╟─629f429f-8b33-4f25-af62-823c89da67a5
# ╟─cb013904-1a68-4655-848b-4d520e66fc48
# ╟─b345fb68-3b7f-4ec2-aa85-77f3f5b94fa1
# ╠═a56195b7-3170-4c59-99d1-1ea8bae7bf69
# ╟─bd6f39b4-a999-44d2-a8b0-646c1e0a26a3
# ╠═d9f14d74-164f-42ab-aac3-4d52a1b8917b
