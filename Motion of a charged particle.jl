### A Pluto.jl notebook ###
# v0.20.24

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
md"# Motion of single charged particle"

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
md"Generally true that ``\mathbf{v} = \mathbf{v}_g + \mathbf{u}``."

# ╔═╡ db44bb1c-7263-4fb8-82cc-effefab75153
md"## Control panel"

# ╔═╡ 2d8bf6cc-0751-4b4e-bd3c-57ef04448d20
guiding_center_plot_checkbox = @bind do_guiding_center CheckBox(default = true);

# ╔═╡ baf9b28c-b16e-421a-8813-572bfb2cab35
projection_plot_checkbox = @bind do_proj CheckBox(default = true);

# ╔═╡ acdea016-7f32-4f37-a5dc-8e7998a6c3b9
charge_selector = @bind q Select([u"q" => "+e", -u"q" => "-e"]);

# ╔═╡ 44eaf06c-1681-42f9-90e5-958bc822d760
mass_selector = @bind m Select([u"mp" => "proton mass", u"me" => "electron mass"]);

# ╔═╡ 854ba5e9-5eae-42a7-8378-c4c5d01e306f
md"Controllers for initial velocity"

# ╔═╡ 59d363ae-24e9-4b5d-b83b-b9df4f4203ad
γ₀_slider = @bind γ₀ Slider((0:359)u"°", show_value = true);

# ╔═╡ f909adaa-cd08-47c6-8adb-ef2a3796ffe3
md"Plot projection: $projection_plot_checkbox"

# ╔═╡ 16a381a7-e3ac-42be-b94d-05a1fba351e5
md"Plot guiding center: $guiding_center_plot_checkbox"

# ╔═╡ 88ada895-2d26-46b4-abfa-4e0dd7187339
md"Time span:"

# ╔═╡ 6068010f-df5c-41a8-94d3-288eb71220d5
tspan = range(0.0u"s", 30.0u"s", length = 601);

# ╔═╡ ee41070b-b68b-4f5f-98b6-c9843f9587eb
md"## Constant and uniform **B** field"

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
md"Solution:"

# ╔═╡ 6743bdc7-8f4f-494f-bdf8-33002a2b68b6
md"""
```math
\begin{align}
\mathbf{x}(t) &= \mathbf{x}_0 + \boldsymbol{ρ} + \boldsymbol{ξ}(t)
= \mathbf{R}_C(t) + \boldsymbol{ξ}(t) \\[1.5ex]
&=
\begin{pmatrix} x_0 \\ y_0 \\ z_0 \end{pmatrix}
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
md"Choose system parameters:"

# ╔═╡ b2c4a082-0b67-4dee-90c2-849e73edf3cc
md"Particle charge = $charge_selector"

# ╔═╡ 771624c1-ade2-4584-a841-3db0132d1be6
md"Particle mass = $mass_selector"

# ╔═╡ 9a791de7-144a-4537-b258-08613a1920a3
md"Magnetic field strength:"

# ╔═╡ 36146cf5-c26a-4759-839a-13c2adc29d96
B = 30.0u"nT";

# ╔═╡ 73d642db-b214-45c3-a16f-558c18a30011
ω_c = q * B / m |> u"rad/s"

# ╔═╡ f10f18d0-8f85-44be-8b5a-474fc2a90d2a
md"Choose initial conditions:"

# ╔═╡ 1a8691be-4fa9-4a21-afcc-10d0d6baf716
md"Initial position:"

# ╔═╡ 62dce489-601c-4d04-a999-068face9b775
x₀ = SVector(0.3, 0.3, 0.3)u"m";

# ╔═╡ c6243506-4a2a-4f37-ad15-a33e4d50d7ea
v_range = range(0.0, 500, step = 0.1)u"m/s"

# ╔═╡ 45b2d155-1a54-4169-8c0f-15e4690a2cdf
v₀_perp_slider = @bind v_perp Slider(v_range, show_value = true, default = 1u"m/s");

# ╔═╡ f0d527de-413a-4a8a-89d6-3ff63f927996
v₀_parallel_slider = @bind v_parallel Slider(v_range, show_value = true, default = 1u"m/s");

# ╔═╡ 32347eb2-3316-4e9c-bd06-434b13e38e2a
md"""
Larmor radius (gyroradius): $(hypot(v_perp, v_parallel)/ω_c |> u"m")

Gyroperiod: $(2π*u"rad"/ω_c)
"""

# ╔═╡ e2690fc7-39f8-436f-973b-ef19a28c8fa0
gcbase = Point3.(0u"m", 0u"m", v_perp .* tspan);

# ╔═╡ b06f0f86-eaa1-4dbb-9c21-1baa251d4c10
md"""
Initial velocity:
- Parallel component: ``v_{0∥}`` = $(v₀_parallel_slider)
- Perpendicular component: ``v_{0⟂}`` = $(v₀_perp_slider)
- Velocity phase angle: ``γ_0`` = $(γ₀_slider)
"""

# ╔═╡ e3382bd2-d27f-4434-8417-4d91cbb7b9be
v₀ = SVector(v_perp*cos(γ₀), v_perp*sin(γ₀), v_parallel)

# ╔═╡ 4687d979-afe8-4116-b4e7-614daa2d46be
md"Plot projection: $projection_plot_checkbox"

# ╔═╡ 5659e5a1-ba89-4ab7-86a2-bbc59ef48903
md"Plot guiding center: $guiding_center_plot_checkbox"

# ╔═╡ 71069a27-0467-4352-9aff-20db6a27a663
md"## Constant and uniform **B** field and **F**"

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
md"Solution:"

# ╔═╡ cd2517b7-8507-45f3-8d52-ce8e33818bae
md"""
```math
\mathbf{x}(t) = \begin{pmatrix}
    x_0 + \frac{v_{0⟂}}{ω_c} \cos(γ_0) - \frac{v_{0⟂}}{ω_c} \cos(ω_c t + γ_0) + \frac{F_y}{qB} t \\
    y_0 - \frac{v_{0⟂}}{ω_c} \sin(γ_0) + \frac{v_{0⟂}}{ω_c} \sin(ω_c t + γ_0) - \frac{F_x}{qB} t \\
    z_0 + v_{0∥}t + \frac{F_\parallel}{2m} t^2
\end{pmatrix}
```
"""

# ╔═╡ bde4ef61-4c56-440f-9082-8decd9b0137a
md"""
```math
\begin{align}
\mathbf{v}(t)
&= \mathbf{u} + \underbrace{\mathbf{v}_\text{D} + \left(\frac{F_∥}{m} t + v_{0∥}\right) \hat{\mathbf{z}}}_{\mathbf{v}_g} \\
&= \begin{pmatrix}
    v_{0⟂} \sin(ω_c t + γ_0) \\
    v_{0⟂} \cos(ω_c t + γ_0) \\
    0 \end{pmatrix}
+
\begin{pmatrix} \hphantom{-} \frac{F_y}{qB} \\ -\frac{F_x}{qB} \\ 0 \end{pmatrix}
+
\begin{pmatrix} 0 \\ 0 \\ \frac{F_∥}{m} t + v_{0∥} \end{pmatrix}
\end{align}
```
"""

# ╔═╡ 1b93bdea-c76f-4bc2-9b89-08217ca8249d
md"""
where ``\mathbf{v}_\text{D} = \frac{1}{qB} (F_y \hat{\mathbf{x}} - F_x \hat{\mathbf{y}})``.
"""

# ╔═╡ f9c390fd-f0b9-471f-a4fc-10f8125f9388
md"Choose force vector:"

# ╔═╡ db7a3387-57e5-43bf-b077-83ca08c585a5
F = SVector(0.5e-23, 8.4e-24, 3.3e-26)u"μN"

# ╔═╡ 3fe0112b-ca82-447b-b086-6370b9530319
md"""
Decompose the guiding center ``\mathbf{r}_g(t)`` into its additive parts:
```math
\mathbf{r}_g(t) = \mathbf{r}_1 + \mathbf{r}_2(t) + \mathbf{r}_3(t) + \mathbf{r}_4(t)
```
where
```math
\begin{align*}
&\mathbf{r}_1 = \frac{v_⟂}{ω_c} \begin{pmatrix} \hphantom{-} \cos(γ_0) \\ - \sin(γ_0) \\ 0 \end{pmatrix}, &
&\mathbf{r}_2(t) = \begin{pmatrix} 0 \\ 0 \\ v_∥ t \end{pmatrix} = v_∥ t \hat{\mathbf{z}}, \\
&\mathbf{r}_3(t) = \frac{1}{qB} \begin{pmatrix}  \hphantom{-} F_y t \\ - F_x t \\ 0 \end{pmatrix}, &
&\mathbf{r}_4(t) = \begin{pmatrix} 0 \\ 0 \\ \frac{F_∥}{2m} t^2 \end{pmatrix} = \frac{F_∥}{2m} t^2 \hat{\mathbf{z}}.
\end{align*}
```
"""

# ╔═╡ 538efe14-e0c1-427a-8440-1870777c0618
r1 = v_perp/ω_c * Point3(cos(γ₀), -sin(γ₀), 0);

# ╔═╡ 2ec96af8-4e67-4756-8236-71b10c6f146e
r2 = Ref(Point3(0u"m/s", 0u"m/s", v_parallel)) .* tspan;

# ╔═╡ 427c7557-f201-4279-8c67-d080f180b6d3
r3F = Ref(Point3(F.y, -F.x, 0u"N")/(q*B)) .* tspan;

# ╔═╡ 82ba8cd1-c7d3-4b29-8ec8-cfd5999d3b84
r4F = Ref(Point3(0u"m/s^2", 0u"m/s^2", F.z/2m)) .* tspan.^2;

# ╔═╡ f1f7aade-2a49-4a92-be9c-889d74e39f1d
gcF = Ref(r1) .+ r2 .+ r3F .+ r4F;

# ╔═╡ 68833be2-9688-4a16-ab58-51606ed8115b
md"""
Initial velocity:
- Parallel component: ``v_{0∥}`` = $(v₀_parallel_slider)
- Perpendicular component: ``v_{0⟂}`` = $(v₀_perp_slider)
- Velocity phase angle: ``γ_0`` = $(γ₀_slider)
"""

# ╔═╡ 79e7ffab-7643-4ae8-88ab-c99ae649060a
md"Plot projection: $projection_plot_checkbox"

# ╔═╡ cb5fdeae-d036-48da-9e5b-7d56a8d6bdb0
md"Plot guiding center: $guiding_center_plot_checkbox"

# ╔═╡ 8f9070d9-9412-4ccd-b84a-d8a20b6f7059
md"## Constant and uniform **E** and **B** fields"

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
md"Solution:"

# ╔═╡ b345fb68-3b7f-4ec2-aa85-77f3f5b94fa1
md"""
```math
\mathbf{x}(t) = \begin{pmatrix}
    x_0 + \frac{v_{0⟂}}{ω_c} \cos(γ_0) - \frac{v_{0⟂}}{ω_c} \cos(ω_c t + γ_0) + \frac{E_y}{B}t \\
    y_0 - \frac{v_{0⟂}}{ω_c} \sin(γ_0) + \frac{v_{0⟂}}{ω_c} \sin(ω_c t + γ_0) - \frac{E_x}{B} t \\
    z_0 + v_{0∥}t + \frac{qE_∥}{2m} t^2
\end{pmatrix}
```
"""

# ╔═╡ a56195b7-3170-4c59-99d1-1ea8bae7bf69
md"""
```math
\begin{align}
\mathbf{v}(t)
&= \mathbf{u} + \mathbf{v}_E + \left(\frac{qE_∥}{m} t + v_{0∥}\right) \hat{\mathbf{z}} \\
&= \begin{pmatrix} v_{0⟂} \sin(ω_c t + γ_0) \\ v_{0⟂} \cos(ω_c t + γ_0) \\ 0 \end{pmatrix}
+ \begin{pmatrix} \hphantom{-} \frac{E_y}{B} \\ - \frac{E_x}{B} \\ 0 \end{pmatrix}
+ \begin{pmatrix} 0 \\ 0 \\ \frac{qE_∥}{m} t + v_{0∥} \end{pmatrix}
\end{align}
```
where ``\mathbf{v}_E = \frac{1}{B} \left(E_y \hat{\mathbf{x}} - E_x \hat{\mathbf{y}}\right)``.
"""

# ╔═╡ 90a32382-7d24-400d-aa06-545b29a53177
E_perp = Point2f(0.0, 0.0003)u"mV/m";

# ╔═╡ 75baec21-d3fa-487a-b8da-e41d88e35c8d
E_parallel = 0.0u"V/m";

# ╔═╡ ad936609-ba30-48a0-9bf7-491b8d963d1c
E = SVector(E_perp[1], E_perp[2], E_parallel) .|> u"V/m"

# ╔═╡ 11f230a0-2384-47e6-b88c-a95cad2d3c85
md"""
Decompose the guiding center ``\mathbf{r}_g(t)`` into its additive parts:
```math
\mathbf{r}_g(t) = \mathbf{r}_1 + \mathbf{r}_2(t) + \mathbf{r}_3(t) + \mathbf{r}_4(t)
```
where
```math
\begin{align*}
&\mathbf{r}_1 = \frac{v_{0⟂}}{ω_c} \begin{pmatrix} \hphantom{-} \cos(γ_0) \\ - \sin(γ_0) \\ 0 \end{pmatrix}, &
&\mathbf{r}_2(t) = \begin{pmatrix} 0 \\ 0 \\ v_∥ t \end{pmatrix} = v_∥ t \hat{\mathbf{z}}, \\
&\mathbf{r}_3(t) = \frac{1}{B} \begin{pmatrix}  \hphantom{-} E_y t \\ - E_x t \\ 0 \end{pmatrix}, &
&\mathbf{r}_4(t) = \begin{pmatrix} 0 \\ 0 \\ \frac{qE_∥}{2m} t^2 \end{pmatrix} = \frac{qE_∥}{2m} t^2 \hat{\mathbf{z}}.
\end{align*}
```
"""

# ╔═╡ a1fc3ec1-9b51-4fee-80e8-40b1d309df59
r3E = Ref(SVector(E.y, -E.x, 0u"V/m")/B) .* tspan;

# ╔═╡ 40dd493c-21fb-4a29-a7b3-89a884345ce5
r4E = Ref(SVector(0u"m/s^2", 0u"m/s^2", q*E.z/2m)) .* tspan.^2;

# ╔═╡ 3d48e09d-95bf-45a8-a3e2-b6595f9e6441
md"""
Initial velocity:
- Parallel component: ``v_{0∥}`` = $(v₀_parallel_slider)
- Perpendicular component: ``v_{0⟂}`` = $(v₀_perp_slider)
- Velocity phase angle: ``γ_0`` = $(γ₀_slider)
"""

# ╔═╡ d99f34e3-8c27-4362-b745-64b0abf64f21
md"Plot projection: $projection_plot_checkbox"

# ╔═╡ ea8f6887-368a-4d13-89f8-d0f60e29af1d
md"Plot guiding center: $guiding_center_plot_checkbox"

# ╔═╡ f6343498-ad2e-4dc8-a015-afe24ad4e129
gcE = Point3.(Ref(r1) .+ r2 .+ r3E .+ r4E);

# ╔═╡ bd6f39b4-a999-44d2-a8b0-646c1e0a26a3
md"## Constant and uniform **F**, **E**, and **B**"

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
    x_0 + \frac{v_{0⟂}}{ω_c} \cos(γ_0) - \frac{v_{0⟂}}{ω_c} \cos(ω_c t + γ_0) + \frac{E_y}{B}t + \frac{F_y}{qB}t \\
    y_0 - \frac{v_{0⟂}}{ω_c} \sin(γ_0) + \frac{v_{0⟂}}{ω_c} \sin(ω_c t + γ_0) - \frac{E_x}{B} t - \frac{F_y}{qB}t \\
    z_0 + v_{0∥}t + \frac{qE_∥ + F_∥}{2m} t^2
\end{pmatrix}
```
"""

# ╔═╡ be32cbaa-a55f-4440-b8a5-c16c2d43cb7d
md"""
```math
\begin{align}
\mathbf{v}(t)
&= \mathbf{u} + \mathbf{v}_\text{D} + \mathbf{v}_E + \left(\frac{qE_∥ + F_∥}{m} t + v_{0∥}\right) \hat{\mathbf{z}} \\[1ex]
&= \begin{pmatrix}
    v_{0⟂} \sin(ω_c t + γ_0) \\
    v_{0⟂} \cos(ω_c t + γ_0) \\
    0
\end{pmatrix}
+ \begin{pmatrix} \hphantom{-} \frac{F_y}{qB} \\ - \frac{F_x}{qB} \\ 0 \end{pmatrix}
+ \begin{pmatrix} \hphantom{-} \frac{E_y}{B} \\ - \frac{E_x}{B} \\ 0 \end{pmatrix}
+ \begin{pmatrix} 0 \\ 0 \\ \frac{qE_∥ + F_∥}{m} t + v_{0∥} \end{pmatrix}
\end{align}
```
"""

# ╔═╡ 1ab6fd2b-6c1a-4a6d-954c-4b181735de71
md"""
Initial velocity:
- Parallel component: ``v_{0∥}`` = $(v₀_parallel_slider)
- Perpendicular component: ``v_{0⟂}`` = $(v₀_perp_slider)
- Velocity phase angle: ``γ_0`` = $(γ₀_slider)
"""

# ╔═╡ 01781a48-99f5-4023-abf5-df6cadf35597
gcEF = Point3.(Ref(r1) .+ r2 .+ r3E .+ r4E .+ r3F .+ r4F);

# ╔═╡ 135b50ad-9a9c-45d3-af9b-e0fae2d4e71e
md"Plot projection: $projection_plot_checkbox"

# ╔═╡ 04e8bf8e-e5d3-40b9-91f2-0ce58fb16905
md"Plot guiding center: $guiding_center_plot_checkbox"

# ╔═╡ da308570-8d28-4d9f-8986-73706ab56983
md"## Helper functions"

# ╔═╡ d9f14d74-164f-42ab-aac3-4d52a1b8917b
function guidingcenter(time, params)
    # TODO generalize for non-zero E, Fₑₓₜ
    v_perp = hypot(v₀.x, v₀.y)
    (; ω_c, γ₀) = params
    return x₀ + v_perp/ω_c * SVector(cos(γ₀), -sin(γ₀), 0)
end

# ╔═╡ 44ca7c41-ecdc-4f32-99b4-834b89e60ad6
"""
    track_motion(fig, ax, trajectory)

Track motion of a single trajectory (one particle only).

### Arguments
- `fig`: Makie `Figure`
- `ax`: Makie `Axis`/`Axis3`
- `trajectory`: list of points
"""
function track_motion(fig, ax, trajectory::AbstractVector{T}) where T <: Union{Point, SVector}
    live_traj = Observable([trajectory[begin]])
    lead_point = Observable(trajectory[begin])
    sizehint!(live_traj[], length(trajectory))
    # draw the start of this new trajectory
    lines!(ax, live_traj)
    scatter!(ax, lead_point)

    r = Record(fig) do io
        # over each timestep/frame
        for r in trajectory
            # update leading point
            lead_point[] = r
            # add new point to trajectory
            live_traj[] = push!(live_traj[], r)
            recordframe!(io)
        end
    end
    return r
end

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
        # over each timestep/frame
        for (i, t_i) in enumerate(t)
            for k_traj in 1:n
                x_i = trajectories[k_traj][i]
                # update leading point
                lead_points[k_traj][] = x_i
                # add new point to trajectory
                push!(live_traj[k_traj][], x_i)
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
    return Point3(
        -v_perp/ω_c * cos(ω_c * time + γ₀),
        -v_perp/ω_c * sin(ω_c * time + γ₀),
        zero(x₀.z))
end

# ╔═╡ fb68ec9d-9112-47f2-afdd-5cd772e61e49
ξ = gyromotion.(tspan, Ref((; ω_c, γ₀, v_perp, v_parallel, x₀, v₀)));

# ╔═╡ 6587b76e-8b93-4d0a-96f8-5bbcd2faafeb
xbase = Point3.(gcbase + ξ);

# ╔═╡ ce558592-d7c8-4723-aa9a-cb7956e5aed3
xF = Point3.(gcF + ξ);

# ╔═╡ 35d14f11-4363-4392-a4e2-8bb9ff6db174
xE = Point3.(gcE + ξ);

# ╔═╡ 531171f9-5b84-4ff8-891a-f66f3fe4d58e
xEF = Point3.(gcEF + ξ);

# ╔═╡ 2509dad4-0058-42f5-8a24-1dbec0b6772b
function bounds(points, eps = 0.1)
    d = length(eltype(points))
    T = eltype(eltype(points))
    bounds = Vector{Tuple{T,T}}(undef, d)
    for i in 1:d
        min, max = extrema(x[i] for x in points)
        min, max = min-T(eps), max+T(eps)   # add some padding
        bounds[i] = (min, max)
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
projectxy(p::Point3{T}) where T = Point3(p[1], p[2], zero(T))

# ╔═╡ 953d901e-63c8-4987-8655-a304521e34bb
projectionbase = projectxy.(xbase);

# ╔═╡ c4914c73-bbf5-4dab-861a-db045f0fc5dd
let
    gc = ustrip.(u"m", gcbase)
    x = ustrip.(u"m", xbase)
    projection = ustrip.(u"m", projectionbase)

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

# ╔═╡ 42023bec-6832-4d27-8bca-17c637fd2f99
projectionF = projectxy.(xF);

# ╔═╡ 2d038925-3025-415d-9198-750101aeb703
let
    gc = ustrip.(u"m", gcF)
    x = ustrip.(u"m", xF)
    projection = ustrip.(u"m", projectionF)

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

# ╔═╡ 21b479b1-27bb-4cdb-9da6-84f4265db7ab
projectionE = projectxy.(xE);

# ╔═╡ 99274bb9-7598-455a-b450-c02782184e29
let
    x = ustrip.(u"m", xE)
    gc = ustrip.(u"m", gcE)
    projection = ustrip.(u"m", projectionE)

    fig = Figure()
    ax = Axis3(fig[1,1])

    b = bounds(vcat(x, gc, projection))
    limits!(ax, b...)

    trajectories = (x,)
    labels = ["Particle"]
    if do_proj
        trajectories = (trajectories..., projection)
        push!(labels, "Projection")
    end
    if do_guiding_center
        trajectories = (trajectories..., gc)
        push!(labels, "Gyrocenter")
    end
    # leg = Legend(fig[1,2], ax, labels)
    r = track_motion(fig, ax, trajectories)
end

# ╔═╡ f9793228-8d2a-4a00-867e-a72988d90ec8
projectionEF = projectxy.(xEF);

# ╔═╡ 2208063b-ce6f-43d8-97b1-5274203c87d3
let
    x = ustrip.(u"m", xEF)
    gc = ustrip.(u"m", gcEF)
    projection = ustrip.(u"m", projectionEF)

    fig = Figure()
    ax = Axis3(fig[1,1])

    b = bounds(vcat(x, gc, projection))
    limits!(ax, b...)

    trajectories = (x,)
    labels = ["Particle"]
    if do_proj
        trajectories = (trajectories..., projection)
        push!(labels, "Projection")
    end
    if do_guiding_center
        trajectories = (trajectories..., gc)
        push!(labels, "Gyrocenter")
    end
    # leg = Legend(fig[1,2], ax, labels)
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
# ╟─854ba5e9-5eae-42a7-8378-c4c5d01e306f
# ╠═45b2d155-1a54-4169-8c0f-15e4690a2cdf
# ╠═f0d527de-413a-4a8a-89d6-3ff63f927996
# ╠═59d363ae-24e9-4b5d-b83b-b9df4f4203ad
# ╟─f909adaa-cd08-47c6-8adb-ef2a3796ffe3
# ╟─16a381a7-e3ac-42be-b94d-05a1fba351e5
# ╟─88ada895-2d26-46b4-abfa-4e0dd7187339
# ╠═6068010f-df5c-41a8-94d3-288eb71220d5
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
# ╠═73d642db-b214-45c3-a16f-558c18a30011
# ╟─32347eb2-3316-4e9c-bd06-434b13e38e2a
# ╟─f10f18d0-8f85-44be-8b5a-474fc2a90d2a
# ╟─1a8691be-4fa9-4a21-afcc-10d0d6baf716
# ╠═62dce489-601c-4d04-a999-068face9b775
# ╟─c6243506-4a2a-4f37-ad15-a33e4d50d7ea
# ╠═e2690fc7-39f8-436f-973b-ef19a28c8fa0
# ╠═fb68ec9d-9112-47f2-afdd-5cd772e61e49
# ╠═6587b76e-8b93-4d0a-96f8-5bbcd2faafeb
# ╠═953d901e-63c8-4987-8655-a304521e34bb
# ╟─b06f0f86-eaa1-4dbb-9c21-1baa251d4c10
# ╠═e3382bd2-d27f-4434-8417-4d91cbb7b9be
# ╟─4687d979-afe8-4116-b4e7-614daa2d46be
# ╟─5659e5a1-ba89-4ab7-86a2-bbc59ef48903
# ╟─c4914c73-bbf5-4dab-861a-db045f0fc5dd
# ╟─71069a27-0467-4352-9aff-20db6a27a663
# ╟─6c2aebe1-9893-42b2-850e-2a1bed21180d
# ╟─c8f06a99-2adb-406c-aa32-c3347616e480
# ╟─cd2517b7-8507-45f3-8d52-ce8e33818bae
# ╟─bde4ef61-4c56-440f-9082-8decd9b0137a
# ╟─1b93bdea-c76f-4bc2-9b89-08217ca8249d
# ╟─f9c390fd-f0b9-471f-a4fc-10f8125f9388
# ╠═db7a3387-57e5-43bf-b077-83ca08c585a5
# ╟─3fe0112b-ca82-447b-b086-6370b9530319
# ╠═538efe14-e0c1-427a-8440-1870777c0618
# ╠═2ec96af8-4e67-4756-8236-71b10c6f146e
# ╠═427c7557-f201-4279-8c67-d080f180b6d3
# ╠═82ba8cd1-c7d3-4b29-8ec8-cfd5999d3b84
# ╠═f1f7aade-2a49-4a92-be9c-889d74e39f1d
# ╠═ce558592-d7c8-4723-aa9a-cb7956e5aed3
# ╠═42023bec-6832-4d27-8bca-17c637fd2f99
# ╟─68833be2-9688-4a16-ab58-51606ed8115b
# ╟─79e7ffab-7643-4ae8-88ab-c99ae649060a
# ╟─cb5fdeae-d036-48da-9e5b-7d56a8d6bdb0
# ╟─2d038925-3025-415d-9198-750101aeb703
# ╟─8f9070d9-9412-4ccd-b84a-d8a20b6f7059
# ╟─629f429f-8b33-4f25-af62-823c89da67a5
# ╟─cb013904-1a68-4655-848b-4d520e66fc48
# ╟─b345fb68-3b7f-4ec2-aa85-77f3f5b94fa1
# ╟─a56195b7-3170-4c59-99d1-1ea8bae7bf69
# ╠═90a32382-7d24-400d-aa06-545b29a53177
# ╠═75baec21-d3fa-487a-b8da-e41d88e35c8d
# ╠═ad936609-ba30-48a0-9bf7-491b8d963d1c
# ╟─11f230a0-2384-47e6-b88c-a95cad2d3c85
# ╠═a1fc3ec1-9b51-4fee-80e8-40b1d309df59
# ╠═40dd493c-21fb-4a29-a7b3-89a884345ce5
# ╠═f6343498-ad2e-4dc8-a015-afe24ad4e129
# ╠═35d14f11-4363-4392-a4e2-8bb9ff6db174
# ╠═21b479b1-27bb-4cdb-9da6-84f4265db7ab
# ╟─3d48e09d-95bf-45a8-a3e2-b6595f9e6441
# ╟─d99f34e3-8c27-4362-b745-64b0abf64f21
# ╟─ea8f6887-368a-4d13-89f8-d0f60e29af1d
# ╟─99274bb9-7598-455a-b450-c02782184e29
# ╟─bd6f39b4-a999-44d2-a8b0-646c1e0a26a3
# ╟─22f446e1-e229-41f2-bb32-85fa8b799e28
# ╟─60329b62-cb14-4471-8177-4b0537a99bfb
# ╟─a07549ec-36f0-4cb1-afe6-362684b2eedb
# ╟─be32cbaa-a55f-4440-b8a5-c16c2d43cb7d
# ╠═01781a48-99f5-4023-abf5-df6cadf35597
# ╠═531171f9-5b84-4ff8-891a-f66f3fe4d58e
# ╠═f9793228-8d2a-4a00-867e-a72988d90ec8
# ╟─1ab6fd2b-6c1a-4a6d-954c-4b181735de71
# ╟─135b50ad-9a9c-45d3-af9b-e0fae2d4e71e
# ╟─04e8bf8e-e5d3-40b9-91f2-0ce58fb16905
# ╟─2208063b-ce6f-43d8-97b1-5274203c87d3
# ╟─da308570-8d28-4d9f-8986-73706ab56983
# ╠═d9f14d74-164f-42ab-aac3-4d52a1b8917b
# ╟─44ca7c41-ecdc-4f32-99b4-834b89e60ad6
# ╟─135f1d0d-8403-47c5-9aae-5efa2b995baf
# ╟─6d736b0b-9e26-4800-8925-f287fc5d292b
# ╟─2509dad4-0058-42f5-8a24-1dbec0b6772b
# ╠═9c4ad1df-9076-4193-97d9-4e058776443c
# ╟─27707139-f329-4908-9342-d3acc560b7d7
