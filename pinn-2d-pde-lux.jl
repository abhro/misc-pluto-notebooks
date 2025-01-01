### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 6a3c45a4-c856-11ef-21e8-2f3a5eecf4da
using PlutoUI

# ╔═╡ 98c588bf-4496-4764-966d-211778c90d4f
using ADTypes, Lux, Optimisers, Zygote

# ╔═╡ b220fb81-9813-47f0-8d9f-de40f59980b6
using Random, Statistics, MLUtils, OnlineStats

# ╔═╡ bcca1d71-9f5b-4969-b66f-984a2c8a0214
using Printf

# ╔═╡ 4ae6eda5-0efc-4a38-a126-3b4d71066fb7
using CairoMakie

# ╔═╡ bff355d3-e47f-45cf-9737-b72058651335
using CUDA, LuxCUDA; CUDA.allowscalar(false)

# ╔═╡ 7d404ce3-2ffe-47d4-9620-04dfd4c84e77
md"""
# Training a PINN on 2D PDE
"""

# ╔═╡ b9baa435-91d5-4a93-b4ac-4b2c4086aaa1
TableOfContents()

# ╔═╡ 06447c74-512b-4607-9ab8-0f96ed7031f8
md"""
Original tutorial from Lux.jl at <https://lux.csail.mit.edu/stable/tutorials/intermediate/4_PINN2DPDE>
"""

# ╔═╡ 12cb9044-82cd-487d-bac8-dee10f62e5f8
md"""
## Package Imports
"""

# ╔═╡ cd15c58f-429f-41e5-a04a-293680fb1104
const gdev = gpu_device();

# ╔═╡ cb734e39-fdca-45fb-abb9-a093e2d3a33e
const cdev = cpu_device();

# ╔═╡ 0c7e1a3b-1483-4094-a615-60c1747e40b1
md"""
## Problem Definition
"""

# ╔═╡ 8c9e02d5-2da0-49bc-a668-8fa7d8152c6d
md"""
The networks should take in 3 values and output one scalar.
"""

# ╔═╡ 485ff842-f862-44be-aa8a-f0119b02b958
md"""
## Define the Neural Netwroks
"""

# ╔═╡ 4de2cde0-8af0-44fc-a2fd-dc8aacfdb607
create_mlp(activation, hidden_dims) = Chain(
	Dense(3 => hidden_dims, activation),
	Dense(hidden_dims => hidden_dims, activation),
	Dense(hidden_dims => hidden_dims, activation),
	Dense(hidden_dims => 1))

# ╔═╡ 8c3e3a18-e99f-4cc5-a50f-102c86be019d
begin
	struct PINN{U, V, W} <: Lux.AbstractLuxContainerLayer{(:u, :v, :w)}
		u::U
		v::V
		w::W
	end
	PINN(; hidden_dims::Int = 32) = PINN(
		create_mlp(tanh, hidden_dims),
		create_mlp(tanh, hidden_dims),
		create_mlp(tanh, hidden_dims))
end

# ╔═╡ e9c6ec31-22c0-4de6-a4be-ecbf6fc4d2bf
md"""
## Define the Loss Functions
"""

# ╔═╡ c7ce3ed4-6e94-4eb3-9282-be820a0cf8e9
"""
	PILoss(u, v, w, xyt::AbstractArray)

Physics informed loss function. `u`, `v`, `w` must be a `StatefulLuxLayer`.
"""
@views function PILoss(
	u::StatefulLuxLayer, v::StatefulLuxLayer, w::StatefulLuxLayer,
	xyt::AbstractArray)

	∇u = Zygote.gradient(sum ∘ u, xyt) |> only
	∂u_∂x, ∂u_∂y, ∂u_∂t = ∇u[1:1, :], ∇u[2:2, :], ∇u[3,:, :]
	∂v_∂x = (Zygote.gradient(sum ∘ v, xyt) |> only)[1:1, :]
	v_xyt = v(xyt)
	∂w_∂y = (Zygote.gradient(sum ∘ w, xyt) |> only)[2:2, :]
	w_xyt = w(xyt)

	return (
		mean(abs2, ∂u_∂t - ∂v_∂x - ∂w_∂y) +
		mean(abs2, v_xyt - ∂y_∂x) +
		mean(abs2, w_xyt - ∂u_∂y)
	)
end

# ╔═╡ c5627c81-aadf-4411-92b3-5ae1a5d72555
md"""
Additionally, we need to compute the loss with respect to the boundary conditions.
"""

# ╔═╡ 64192143-4d35-4ec6-bd49-fccb17145bc1
mse_loss_function(
	u::StatefulLuxLayer, target::AbstractArray, xyt::AbstractArray
) = MSELoss()(u(xyt), target);

# ╔═╡ 0889cc63-78a2-4dbb-81cd-17401c936ec2
function loss_function(model, ps, st, (xyt, target, xyt_bc, target_bc))
	u_net = StatefulLuxLayer{true}(model.u, ps.u, st.u)
	v_net = StatefulLuxLayer{true}(model.v, ps.v, st.v)
	w_net = StatefulLuxLayer{true}(model.v, ps.w, st.w)
	physics_loss = PILoss(u_net, v_net, w_net, xyt)
	data_loss = mse_loss_function(u_net, target_data, xyt)
	bc_loss = mse_loss_function(u_net, target_bc, xyt_bc)
	loss = physics_Loss + data_loss + bc_loss

	return (
		loss,
		(; u = u_net.st, v = v_net.st, w = w_net.st),
		(; physics_loss, data_loss, bc_loss),
	)
end

# ╔═╡ edcf06cb-247e-4446-94bd-28c9cfbf11f1
md"""
## Generate the Data
"""

# ╔═╡ 5b01b39e-c503-474c-b688-19466b30f9a9


# ╔═╡ 903fff49-d4a2-48d8-b689-2561e6152a97
md"""
## Training
"""

# ╔═╡ 34efddda-66e0-4541-899c-beb7504a8362
md"""
## Visualizing the Results
"""

# ╔═╡ 6d461a63-9f12-4282-b3d2-fdb352be2ecc
md"""
## Appendix
"""

# ╔═╡ Cell order:
# ╟─7d404ce3-2ffe-47d4-9620-04dfd4c84e77
# ╟─6a3c45a4-c856-11ef-21e8-2f3a5eecf4da
# ╟─b9baa435-91d5-4a93-b4ac-4b2c4086aaa1
# ╟─06447c74-512b-4607-9ab8-0f96ed7031f8
# ╟─12cb9044-82cd-487d-bac8-dee10f62e5f8
# ╠═98c588bf-4496-4764-966d-211778c90d4f
# ╠═b220fb81-9813-47f0-8d9f-de40f59980b6
# ╠═bcca1d71-9f5b-4969-b66f-984a2c8a0214
# ╠═4ae6eda5-0efc-4a38-a126-3b4d71066fb7
# ╠═bff355d3-e47f-45cf-9737-b72058651335
# ╠═cd15c58f-429f-41e5-a04a-293680fb1104
# ╠═cb734e39-fdca-45fb-abb9-a093e2d3a33e
# ╟─0c7e1a3b-1483-4094-a615-60c1747e40b1
# ╟─8c9e02d5-2da0-49bc-a668-8fa7d8152c6d
# ╟─485ff842-f862-44be-aa8a-f0119b02b958
# ╠═8c3e3a18-e99f-4cc5-a50f-102c86be019d
# ╠═4de2cde0-8af0-44fc-a2fd-dc8aacfdb607
# ╟─e9c6ec31-22c0-4de6-a4be-ecbf6fc4d2bf
# ╠═c7ce3ed4-6e94-4eb3-9282-be820a0cf8e9
# ╟─c5627c81-aadf-4411-92b3-5ae1a5d72555
# ╠═64192143-4d35-4ec6-bd49-fccb17145bc1
# ╠═0889cc63-78a2-4dbb-81cd-17401c936ec2
# ╟─edcf06cb-247e-4446-94bd-28c9cfbf11f1
# ╠═5b01b39e-c503-474c-b688-19466b30f9a9
# ╟─903fff49-d4a2-48d8-b689-2561e6152a97
# ╟─34efddda-66e0-4541-899c-beb7504a8362
# ╟─6d461a63-9f12-4282-b3d2-fdb352be2ecc
