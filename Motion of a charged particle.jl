### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
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

# ╔═╡ 91302cb8-87fa-43d4-a41a-db6c409e3378
md"""
List of position variables:
- ``\mathbf{x}(t)``: position of particle
- ``\mathbf{R}_C(t)``: position of gyrocenter
- ``\boldsymbol{ξ}(t)``: position of particle relative to gyrocenter, i.e., ``\boldsymbol{ξ} = \mathbf{x} - \mathbf{R}_C``
- ``\mathbf{x}_0``: initial position of the particle
- ``\mathbf{X}_0``: initial position of the gyrocenter (``\mathbf{X}_0 = \mathbf{x}_0 - \frac{v_{0⟂}}{ω_c} (\cos(γ_0)\hat{\mathbf{x}} - \sin(γ_0)\hat{\mathbf{y}})``)
- ``\boldsymbol{ρ}(t)``: relative displacement of the gyrocenter from initial position? (``\boldsymbol{ρ}(t) = \mathbf{X}_0 - \mathbf{R}_C(t)``)
"""

# ╔═╡ 7429f193-2695-4dcb-be2e-5a9f49716621
md"""
List of velocity variables:
- ``\mathbf{v}``: Velocity of particle
- ``\mathbf{u}``: Gyromotion (velocity of particle relative to guiding center/gyrocenter)
- ``\mathbf{v}_g``: Guiding center velocity
- ``\mathbf{v}_E``: **E**×**B** drift of guiding center
- ``\mathbf{v}_\text{D}``: Guiding center drift
- ``γ_0 = \arctan\left(\dfrac{v_{0,y}}{v_{0,x}}\right)``: Initial phase of the velocity vector in the xy-plane.
"""

# ╔═╡ e9c4f014-11c1-4c6d-b0e0-24ade00f394c
md"""
Generally true that ``\mathbf{v} = \mathbf{v}_g + \mathbf{u}``.
"""

# ╔═╡ db44bb1c-7263-4fb8-82cc-effefab75153
md"""
## Control panel
"""

# ╔═╡ 2d8bf6cc-0751-4b4e-bd3c-57ef04448d20
guiding_center_plot_binder = @bind do_guiding_center CheckBox(default = true);

# ╔═╡ baf9b28c-b16e-421a-8813-572bfb2cab35
projection_plot_binder = @bind do_proj CheckBox(default = true);

# ╔═╡ acdea016-7f32-4f37-a5dc-8e7998a6c3b9
mass_binder = @bind q Select([u"q" => "+e", -u"q" => "-e"]);

# ╔═╡ 44eaf06c-1681-42f9-90e5-958bc822d760
charge_binder = @bind m Select([u"mp" => "proton mass", u"me" => "electron mass"]);

# ╔═╡ f909adaa-cd08-47c6-8adb-ef2a3796ffe3
md"""
Plot projection: $projection_plot_binder
"""

# ╔═╡ 16a381a7-e3ac-42be-b94d-05a1fba351e5
md"""
Plot guiding center: $guiding_center_plot_binder
"""

# ╔═╡ 88ada895-2d26-46b4-abfa-4e0dd7187339
md"""
Time span:
"""

# ╔═╡ 178b6fa9-4066-4acb-8cbc-73bae83930a3
tspan = range(0.0u"s", 300.0u"s", length = 2000);

# ╔═╡ ee41070b-b68b-4f5f-98b6-c9843f9587eb
md"""
## Constant and uniform **B** field
"""

# ╔═╡ 289316fd-59ee-4c9e-9bad-3dd2a93fcf21
md"""
```math
\begin{align}
\mathbf{E} &= \mathbf{0}, &
\mathbf{B} &= B \hat{\mathbf{z}}, &
\mathbf{F}_\text{ext} &= \mathbf{0}.
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
= \mathbf{R}_C(t) + \boldsymbol{ξ}(t) \\[1.5ex]
&= \begin{pmatrix}
    x_0 \\
    y_0 \\
    z_0
\end{pmatrix}
+
\begin{pmatrix}
    \hphantom{-} \frac{v_{0⟂}}{ω_c} \cos(γ_0) - \frac{v_{0⟂}}{ω_c} \cos(ω_c t + γ_0) \\
    - \frac{v_{0⟂}}{ω_c} \sin(γ_0) + \frac{v_{0⟂}}{ω_c} \sin(ω_c t + γ_0) \\
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
= \mathbf{u}(t) + \mathbf{v}_\parallel(t)
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
where
- ``ω_c = {qB}/{m}`` is the gyrofrequency
- ``v_{0⟂} = \sqrt{v_{0,x}^2 + v_{0,y}^2}``
"""

# ╔═╡ a8810170-7c78-4e63-b8aa-bc7c3d6b3f04
md"""
Choose system parameters:
"""

# ╔═╡ b2c4a082-0b67-4dee-90c2-849e73edf3cc
md"""
Particle charge = $mass_binder
"""

# ╔═╡ 771624c1-ade2-4584-a841-3db0132d1be6
md"""
Particle mass = $charge_binder
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

# ╔═╡ 4687d979-afe8-4116-b4e7-614daa2d46be
md"""
Plot projection: $projection_plot_binder
"""

# ╔═╡ 5659e5a1-ba89-4ab7-86a2-bbc59ef48903
md"""
Plot guiding center: $guiding_center_plot_binder
"""

# ╔═╡ 71069a27-0467-4352-9aff-20db6a27a663
md"""
## Constant and uniform **B** field and **F**
"""

# ╔═╡ 6c2aebe1-9893-42b2-850e-2a1bed21180d
md"""
```math
\begin{align}
\mathbf{E} &= \mathbf{0}, &
\mathbf{B} &= B \hat{\mathbf{z}}, &
\mathbf{F}_\text{ext} &= \mathbf{F}_⟂ + \mathbf{F}_∥,
\end{align}
```
where ``\mathbf{F}_⟂ = F_x \hat{\mathbf{x}} + F_y \hat{\mathbf{y}}`` and ``\mathbf{F}_∥ = F_z \hat{\mathbf{z}}``.
"""

# ╔═╡ c8f06a99-2adb-406c-aa32-c3347616e480
md"""
Solution:
"""

# ╔═╡ cd2517b7-8507-45f3-8d52-ce8e33818bae
md"""
```math
\mathbf{x}(t) = \begin{pmatrix}
    x_0 + \frac{v_{0⟂}}{ω_c} \cos(γ_0) - \frac{v_{0⟂}}{ω_c} \cos(ω_c t + γ_0) + \frac{F_y}{qB} t \\
    y_0 - \frac{v_{0⟂}}{ω_c} \sin(γ_0) + \frac{v_{0⟂}}{ω_c} \sin(ω_c t + γ_0) - \frac{F_x}{qB} t \\
    z_0 + v_{\parallel,0}t + \frac{F_\parallel}{2m} t^2
\end{pmatrix}
```
"""

# ╔═╡ bde4ef61-4c56-440f-9082-8decd9b0137a
md"""
```math
\begin{align}
\mathbf{v}(t)
&= \mathbf{u} + \underbrace{\mathbf{v}_\text{D} + \left(\frac{F_∥}{m} t + v_{∥0}\right) \hat{\mathbf{z}}}_{\mathbf{v}_g} \\
&= \begin{pmatrix}
    v_{0⟂} \sin(ω_c t + γ_0) \\
    v_{0⟂} \cos(ω_c t + γ_0) \\
    0 \end{pmatrix}
+ \begin{pmatrix}
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

# ╔═╡ 1b93bdea-c76f-4bc2-9b89-08217ca8249d
md"""
where ``\mathbf{v}_\text{D} = \frac{1}{qB} (F_y \hat{\mathbf{x}} - F_x \hat{\mathbf{y}})``.
"""

# ╔═╡ f9c390fd-f0b9-471f-a4fc-10f8125f9388
md"""
Choose force vector:
"""

# ╔═╡ db7a3387-57e5-43bf-b077-83ca08c585a5
F = SVector(0.5, 0.3, 3.3)u"N"

# ╔═╡ 538efe14-e0c1-427a-8440-1870777c0618
r1 = v_perp/ω_c * SVector(cos(γ₀), -sin(γ₀), 0);

# ╔═╡ 427c7557-f201-4279-8c67-d080f180b6d3
r2 = Ref(SVector(F.y, -F.x, 0u"N")/(q*B)) .* tspan;

# ╔═╡ 2ec96af8-4e67-4756-8236-71b10c6f146e
r3 = Ref(SVector(0u"m/s", 0u"m/s", v_parallel)) .* tspan;

# ╔═╡ 82ba8cd1-c7d3-4b29-8ec8-cfd5999d3b84
r4 = Ref(SVector(0u"m/s^2", 0u"m/s^2", F.z/2m)) .* tspan.^2;

# ╔═╡ 8f9070d9-9412-4ccd-b84a-d8a20b6f7059
md"""
## Constant and uniform **E** and **B** fields
"""

# ╔═╡ 629f429f-8b33-4f25-af62-823c89da67a5
md"""
```math
\begin{align}
\mathbf{E} &= \mathbf{E}_⟂ + E_∥ \hat{\mathbf{z}}, &
\mathbf{B} &= B \hat{\mathbf{z}}, &
\mathbf{F}_\text{ext} &= 0.
\end{align}
```
where ``\mathbf{E}_⟂ = E_x \hat{\mathbf{x}} + E_y \hat{\mathbf{y}}`` and ``\mathbf{E}_∥ = E_z \hat{\mathbf{z}}``.
"""

# ╔═╡ cb013904-1a68-4655-848b-4d520e66fc48
md"""
Solution:
"""

# ╔═╡ b345fb68-3b7f-4ec2-aa85-77f3f5b94fa1
md"""
```math
\mathbf{x}(t) = \begin{pmatrix}
    x_0 + \frac{v_{0⟂}}{ω_c} \cos(γ_0) - \frac{v_{0⟂}}{ω_c} \cos(ω_c t + γ_0) + \frac{E_y}{B}t \\
    y_0 - \frac{v_{0⟂}}{ω_c} \sin(γ_0) + \frac{v_{0⟂}}{ω_c} \sin(ω_c t + γ_0) - \frac{E_x}{B} t \\
    z_0 + v_{∥,0}t + \frac{qE_∥}{2m} t^2
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
    v_{0⟂} \sin(ω_c t + γ_0) \\
    v_{0⟂} \cos(ω_c t + γ_0) \\
    0
\end{pmatrix}
+ \begin{pmatrix} \hphantom{-} \frac{E_y}{B} \\ - \frac{E_x}{B} \\ 0 \end{pmatrix}
+ \begin{pmatrix} 0 \\ 0 \\ \frac{qE_∥}{m} t + v_{∥0} \end{pmatrix}
\end{align}
```
where ``\mathbf{v}_E = \frac{1}{B} \left(E_y \hat{\mathbf{x}} - E_x \hat{\mathbf{y}}\right)``.
"""

# ╔═╡ bd6f39b4-a999-44d2-a8b0-646c1e0a26a3
md"""
## Constant and uniform **F**, **E**, and **B**
"""

# ╔═╡ 22f446e1-e229-41f2-bb32-85fa8b799e28
md"""
```math
\begin{align}
\mathbf{E} &= \mathbf{E}_∥ + E_⟂ \hat{\mathbf{z}}, &
\mathbf{B} &= B \hat{\mathbf{z}}, &
\mathbf{F}_\text{ext} &=  \mathbf{F}_⟂ + \mathbf{F}_∥,
\end{align}
```
where ``\mathbf{F}_⟂ = F_x \hat{\mathbf{x}} + F_y \hat{\mathbf{y}}`` and ``\mathbf{F}_∥ = F_z \hat{\mathbf{z}}``.
"""

# ╔═╡ 60329b62-cb14-4471-8177-4b0537a99bfb
md"""
Solution:
"""

# ╔═╡ a07549ec-36f0-4cb1-afe6-362684b2eedb
md"""
```math
\mathbf{x}(t) = \begin{pmatrix}
    \\
    \\
\end{pmatrix}
```
"""

# ╔═╡ be32cbaa-a55f-4440-b8a5-c16c2d43cb7d
md"""
```math
\begin{align}
\mathbf{v}(t)
&= \mathbf{u} + \mathbf{v}_\text{D} + \mathbf{v}_E + \left(\frac{qE_∥ + F_∥}{m} t + v_{∥0}\right) \hat{\mathbf{z}} \\[1ex]
&= \begin{pmatrix}
    v_{0⟂} \sin(ω_c t + γ_0) \\
    v_{0⟂} \cos(ω_c t + γ_0) \\
    0
\end{pmatrix}
+ \begin{pmatrix} \hphantom{-} \frac{F_y}{qB} \\ - \frac{F_x}{qB} \\ 0 \end{pmatrix}
+ \begin{pmatrix} \hphantom{-} \frac{E_y}{B} \\ - \frac{E_x}{B} \\ 0 \end{pmatrix}
+ \begin{pmatrix} 0 \\ 0 \\ \frac{qE_∥ + F_∥}{m} t + v_{∥0} \end{pmatrix}
\end{align}
```
"""

# ╔═╡ d9f14d74-164f-42ab-aac3-4d52a1b8917b
function guidingcenter(time, params)
    # TODO generalize for non-zero E, Fₑₓₜ
    v_perp = hypot(v₀.x, v₀.y)
    (; ω_c, γ₀) = params
    return x₀ + v_perp/ω_c * SVector(cos(γ₀), -sin(γ₀), 0)
end


# ╔═╡ da308570-8d28-4d9f-8986-73706ab56983
md"""
## Helper functions
"""

# ╔═╡ 135f1d0d-8403-47c5-9aae-5efa2b995baf
"""
    track_motion(fig, ax, trajectories, t = eachindex(trajectories[begin]))

### Arguments
- `fig`: Makie `Figure`
- `ax`: Makie `Axis`/`Axis3`
- `trajectories`: list of lists, i.e., similar to `Vector{Vector{Point}}`. Each item of trajectories is `Vector{Point}`, which defines a trajectory.
- `t`: parameter with which to define parameteric curve
"""
function track_motion(fig, ax, trajectories, t = eachindex(trajectories[begin]))
    n = length(trajectories)
    live_traj = Observable[]
    lead_points = Observable[]
    for traj in trajectories
        push!(live_traj, Observable(traj[begin:begin+1]))
        push!(lead_points, Observable(traj[begin+1]))

        # draw the start of this new trajectory
        lines!(ax, live_traj[end])
        scatter!(ax, lead_points[end])
    end

    r = Record(fig) do io
        for (i, t_i) in enumerate(t)
            for k_traj in 1:n
                x_i = trajectories[k_traj][i]
                # add new point to trajectory
                push!(live_traj[k_traj][], x_i)
                # update leading point
                lead_points[k_traj][] = x_i
                live_traj[k_traj][] = live_traj[k_traj][]
            end
            recordframe!(io)
        end
    end
    return r
end

# ╔═╡ 6d736b0b-9e26-4800-8925-f287fc5d292b
function gyromotion(time, params)
    (; v_perp, v_parallel, ω_c, γ₀, x₀) = params
    return SVector(
        -v_perp/ω_c * cos(ω_c * time + γ₀),
        -v_perp/ω_c * sin(ω_c * time + γ₀),
        zero(params.x₀.z))
end

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

# ╔═╡ 9c4ad1df-9076-4193-97d9-4e058776443c
Unitful.uconvert(u, p::Point) = uconvert.(u, p)

# ╔═╡ 27707139-f329-4908-9342-d3acc560b7d7
"""
    projectxy(p::Point3)

Project a point in 3d to the xy-plane. (Zero out the last component).
"""
projectxy(p::Point3{T}) where T = Point3(p[1], p[2], zero(T));

# ╔═╡ c4914c73-bbf5-4dab-861a-db045f0fc5dd
let
    gc = SVector.(0u"m", 0u"m", v_perp .* tspan)
    ξ = gyromotion.(tspan, Ref((; ω_c, γ₀, v_perp, v_parallel, x₀, v₀)))
    x = ustrip.(u"m", Point3.(gc + ξ))
    gc = let
        gc_stripped = Vector{Point3f}(undef, length(gc))
        for i in eachindex(gc)
            gc_stripped[i] = ustrip(u"m", Point3(gc[i]))
        end
        gc_stripped
    end

    projection = projectxy.(x)

    fig = Figure()
    ax = Axis3(fig[1,1])

    # leg = Legend(fig[1,2], ax, ["Particle", "Projection", "Gyrocenter"])

    b = bounds(vcat(x, gc, projection))
    limits!(ax, b...)

    trajectories = (x,)
    if do_proj
        trajectories = (trajectories..., projection)
    end
    if do_guiding_center
        trajectories = (trajectories..., gc)
    end

    r = track_motion(fig, ax, trajectories)
end

# ╔═╡ 2d038925-3025-415d-9198-750101aeb703
let
    gc = (Ref(r1) .+ r2 .+ r3 .+ r4)
    # gc = SVector.(0u"m", 0u"m", v_perp .* tspan)
    ξ = gyromotion.(tspan, Ref((; ω_c, γ₀, v_perp, v_parallel, x₀, v₀)))
    x = ustrip.(u"m", Point3.(gc + ξ))
    gc = let
        gc_stripped = Vector{Point3f}(undef, length(gc))
        for i in eachindex(gc)
            gc_stripped[i] = ustrip(u"m", Point3(gc[i]))
        end
        gc_stripped
    end

    projection = projectxy.(x)

    fig = Figure()
    ax = Axis3(fig[1,1])

    # leg = Legend(fig[1,2], ax, ["Particle", "Projection", "Gyrocenter"])

    b = bounds(vcat(x, gc, projection))
    limits!(ax, b...)

    trajectories = (x,)
    if do_proj
        trajectories = (trajectories..., projection)
    end
    if do_guiding_center
        trajectories = (trajectories..., gc)
    end

    r = track_motion(fig, ax, trajectories)
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
# ╟─91302cb8-87fa-43d4-a41a-db6c409e3378
# ╟─7429f193-2695-4dcb-be2e-5a9f49716621
# ╟─e9c4f014-11c1-4c6d-b0e0-24ade00f394c
# ╟─db44bb1c-7263-4fb8-82cc-effefab75153
# ╠═2d8bf6cc-0751-4b4e-bd3c-57ef04448d20
# ╠═baf9b28c-b16e-421a-8813-572bfb2cab35
# ╠═acdea016-7f32-4f37-a5dc-8e7998a6c3b9
# ╠═44eaf06c-1681-42f9-90e5-958bc822d760
# ╟─f909adaa-cd08-47c6-8adb-ef2a3796ffe3
# ╟─16a381a7-e3ac-42be-b94d-05a1fba351e5
# ╟─88ada895-2d26-46b4-abfa-4e0dd7187339
# ╠═178b6fa9-4066-4acb-8cbc-73bae83930a3
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
# ╟─c6243506-4a2a-4f37-ad15-a33e4d50d7ea
# ╟─b06f0f86-eaa1-4dbb-9c21-1baa251d4c10
# ╟─55b22b9e-0727-4ed0-80c5-095175ad7be0
# ╟─98891a54-cee8-4b20-9c62-34aa51557e17
# ╟─ec75d67b-b248-471c-912e-f68096e4133c
# ╠═e3382bd2-d27f-4434-8417-4d91cbb7b9be
# ╟─4687d979-afe8-4116-b4e7-614daa2d46be
# ╟─5659e5a1-ba89-4ab7-86a2-bbc59ef48903
# ╠═c4914c73-bbf5-4dab-861a-db045f0fc5dd
# ╟─71069a27-0467-4352-9aff-20db6a27a663
# ╟─6c2aebe1-9893-42b2-850e-2a1bed21180d
# ╟─c8f06a99-2adb-406c-aa32-c3347616e480
# ╟─cd2517b7-8507-45f3-8d52-ce8e33818bae
# ╟─bde4ef61-4c56-440f-9082-8decd9b0137a
# ╟─1b93bdea-c76f-4bc2-9b89-08217ca8249d
# ╟─f9c390fd-f0b9-471f-a4fc-10f8125f9388
# ╠═db7a3387-57e5-43bf-b077-83ca08c585a5
# ╠═538efe14-e0c1-427a-8440-1870777c0618
# ╠═427c7557-f201-4279-8c67-d080f180b6d3
# ╠═2ec96af8-4e67-4756-8236-71b10c6f146e
# ╠═82ba8cd1-c7d3-4b29-8ec8-cfd5999d3b84
# ╠═2d038925-3025-415d-9198-750101aeb703
# ╟─8f9070d9-9412-4ccd-b84a-d8a20b6f7059
# ╟─629f429f-8b33-4f25-af62-823c89da67a5
# ╟─cb013904-1a68-4655-848b-4d520e66fc48
# ╟─b345fb68-3b7f-4ec2-aa85-77f3f5b94fa1
# ╟─a56195b7-3170-4c59-99d1-1ea8bae7bf69
# ╟─bd6f39b4-a999-44d2-a8b0-646c1e0a26a3
# ╟─22f446e1-e229-41f2-bb32-85fa8b799e28
# ╟─60329b62-cb14-4471-8177-4b0537a99bfb
# ╟─a07549ec-36f0-4cb1-afe6-362684b2eedb
# ╟─be32cbaa-a55f-4440-b8a5-c16c2d43cb7d
# ╠═d9f14d74-164f-42ab-aac3-4d52a1b8917b
# ╟─da308570-8d28-4d9f-8986-73706ab56983
# ╠═135f1d0d-8403-47c5-9aae-5efa2b995baf
# ╠═6d736b0b-9e26-4800-8925-f287fc5d292b
# ╠═2509dad4-0058-42f5-8a24-1dbec0b6772b
# ╠═9c4ad1df-9076-4193-97d9-4e058776443c
# ╠═27707139-f329-4908-9342-d3acc560b7d7
