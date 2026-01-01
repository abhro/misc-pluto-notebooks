### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 22b190dd-fba5-4e33-aad8-e4bba5d08073
using Pkg

# ╔═╡ 6a3c45a4-c856-11ef-21e8-2f3a5eecf4da
using PlutoUI: TableOfContents, with_terminal

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

# ╔═╡ d2d252ab-d68a-4878-908a-ea13d9adf324
using InteractiveUtils

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
    ∂u_∂x, ∂u_∂y, ∂u_∂t = ∇u[1:1, :], ∇u[2:2, :], ∇u[3:3, :]

    ∂v_∂x = (Zygote.gradient(sum ∘ v, xyt) |> only)[1:1, :]
    v_xyt = v(xyt)

    ∂w_∂y = (Zygote.gradient(sum ∘ w, xyt) |> only)[2:2, :]
    w_xyt = w(xyt)

    return (
        mean(abs2, ∂u_∂t - ∂v_∂x - ∂w_∂y) +
        mean(abs2, v_xyt - ∂u_∂x) +
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
    data_loss = mse_loss_function(u_net, target, xyt)
    bc_loss = mse_loss_function(u_net, target_bc, xyt_bc)
    loss = physics_loss + data_loss + bc_loss

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

# ╔═╡ be4aade8-1f7f-46df-b006-d301c0835ac1
analytical_solution(x, y, t) = exp(x + y) * cos(x + y + 4t);

# ╔═╡ 222bd364-af5f-47ee-8120-bea7513ae40a
analytical_solution(xyt) = analytical_solution.(xyt[1, :], xyt[2, :], xyt[3, :]);

# ╔═╡ d0c8d9dc-faf0-41ee-bc5e-e1696d118442
grid_len = 16;

# ╔═╡ cfbd31b8-8882-43c6-9c66-6cacd5d5668e
grid = range(0.0f0, 2.0f0, length = grid_len);

# ╔═╡ 4a73be53-0e19-44dd-ab43-5fced8f188e1
xyt = stack([[elem...] for elem in vec(collect(Iterators.product(grid, grid, grid)))])

# ╔═╡ 7df36b77-65ea-4356-95a7-dd6b025780c4
target_data = reshape(analytical_solution(xyt), 1, :)

# ╔═╡ cf94c793-4fdd-4c31-8a69-dcc0786f98d7
bc_len = 512;

# ╔═╡ 24fcec22-807d-4790-a4a2-9506dc48a044
x = y = t = range(0.0f0, 2.0f0; length = bc_len);

# ╔═╡ ed721e8b-eb50-4d36-b293-9ab7fc1cf31c
xyt_bc = hcat(
    stack((x, y, zeros(Float32, bc_len)), dims = 1),

    stack((zeros(Float32, bc_len), y, t), dims = 1),
    stack((2ones(Float32, bc_len), y, t), dims = 1),

    stack((x, zeros(Float32, bc_len), t), dims = 1),
    stack((x, 2ones(Float32, bc_len), t), dims = 1)
)

# ╔═╡ 46cf6020-33c8-4354-bd4e-164c8743c5e0
target_bc = reshape(analytical_solution(xyt_bc), 1, :)

# ╔═╡ a94c00e3-1ba9-4baf-a31f-fc1528aa788d
min_target_bc, max_target_bc = extrema(target_bc)

# ╔═╡ 1bbcbbc7-41bb-49f5-90f8-2f0552afec10
min_data, max_data = extrema(target_data)

# ╔═╡ c65994ca-0e8b-4c72-85df-c62fb7306461
begin
    """Apply min-max normalization"""
    minmaxnormalize(a, min, max) = (a .- min) / (max - min)
    minmaxnormalize(a) = (a .- minimum(a)) / (maximum(a) - minimum(a))
end

# ╔═╡ 613cdd28-f832-4763-8cef-173c074b51d6
min_pde_val, max_pde_val = min(min_data, min_target_bc), max(max_data, max_target_bc)

# ╔═╡ 576374ec-224b-4f47-b498-00d699e9e5f8
xyt_scaled = minmaxnormalize(xyt)

# ╔═╡ cb0d5b82-eb2f-4698-99b8-6e3ee0cee6d9
xyt_bc_scaled = minmaxnormalize(xyt_bc)

# ╔═╡ 99b730a3-479e-4129-a642-3a71dc8700c2
target_bc_scaled = minmaxnormalize(target_bc, min_pde_val, max_pde_val)

# ╔═╡ 5ff62ed3-ab40-428f-9ba6-8be4ab1bf3da
target_data_scaled = minmaxnormalize(target_data, min_pde_val, max_pde_val)

# ╔═╡ 903fff49-d4a2-48d8-b689-2561e6152a97
md"""
## Training
"""

# ╔═╡ ddd1cbc7-7f1e-4a93-9f4f-07f5a4288a95
function train_model(
    xyt, target_data, xyt_bc, target_bc;
    seed::Int = 0, maxiters::Int = 50_000, hidden_dims::Int = 32)

    rng = Random.default_rng()
    Random.seed!(rng, seed)

    pinn = PINN(; hidden_dims)
    ps, st = Lux.setup(rng, pinn) |> gdev

    bc_dataloader = DataLoader((xyt_bc, target_bc), batchsize = 32, shuffle = true) |> gdev
    pde_dataloader = DataLoader((xyt, target_data), batchsize = 32, shuffle = true) |> gdev

    train_state = Training.TrainState(pinn, ps, st, Adam(5f-2))
    # adaptive learning rate
    lr = i -> i < 5000 ? 5f-2 : (i < 10_000 ? 5f-3 : 5f-4)

    total_loss_tracker, physics_loss_tracker, data_loss_tracker, bc_loss_tracker = ntuple(
        _ -> Lag(Float32, 32), 4)

    iter = 1
    dataiterator = zip(
        Iterators.cycle(pde_dataloader), Iterators.cycle(bc_dataloader))

    for ((xyt_batch, u_batch), (xyt_bc_batch, u_bc_batch)) in dataiterator

        Optimisers.adjust!(train_state, lr(iter))

        _, loss, stats, train_state = Training.single_train_step!(
            AutoZygote(), loss_function,
            (xyt_batch, u_batch, xyt_bc_batch, u_bc_batch),
            train_state)

        fit!(total_loss_tracker, loss)
        fit!(physics_loss_tracker, stats.physics_loss)
        fit!(data_loss_tracker, stats.data_loss)
        fit!(bc_loss_tracker, stats.bc_loss)

        mean_loss = mean(OnlineStats.value(total_loss_tracker))
        mean_physics_loss = mean(OnlineStats.value(physics_loss_tracker))
        mean_data_loss = mean(OnlineStats.value(data_loss_tracker))
        mean_bc_loss = mean(OnlineStats.value(bc_loss_tracker))

        isnan(loss) && throw(ArgumentError("NaN Loss Detected"))

        if iter == 1 || iter % 1000 == 0 || iter == maxiters
            @info(@sprintf(
                "Iteration: [%5d / %5d] \t Loss: %.9f (%.9f) \t Physics Loss: %.9f (%.9f) \t Data Loss: %.9F (%.9f) \t BC Loss: %.9f (%.9f))",
                iter, maxiters, loss, mean_loss,
                stats.physics_loss, mean_physics_loss,
                stats.data_loss, mean_data_loss,
                stats.bc_loss, mean_bc_loss
            ))
        end

        iter ≥ maxiters && break
        iter += 1
    end

    return StatefulLuxLayer{true}(
        pinn, cdev(train_state.parameters), cdev(train_state.states))
end

# ╔═╡ 4fe21a00-c083-4fd6-a695-c52240c3a060
trained_model = train_model(xyt_scaled, target_data_scaled, xyt_bc_scaled, target_bc_scaled)

# ╔═╡ b67d7b58-c79a-45bc-a92e-81be68febf2d
trained_u = Lux.testmode(StatefulLuxLayer{true}(
    trained_model.model.u, trained_model.ps.u, trained_model.st.u))

# ╔═╡ 34efddda-66e0-4541-899c-beb7504a8362
md"""
## Visualizing the Results
"""

# ╔═╡ dea86dbf-0ffd-416a-88dd-b9300b792424
ts = 0.0f0:0.05f0:2.0f0;

# ╔═╡ 1b3b4087-3283-40fc-9eee-62b719e54766
xs = ys = 0.0f0:0.02f0:2.0f0;

# ╔═╡ 24c21e66-f8dd-4327-ad65-e039673ab51b
viz_grid = stack([[elem...] for elem in vec(collect(Iterators.product(xs, ys, ts)))])

# ╔═╡ 01694d01-3ad6-4076-878a-25651fc1d7f6
u_real = reshape(analytical_solution(viz_grid), length(xs), length(ys), length(ts))

# ╔═╡ ccf7acaf-64ed-46b9-9ed1-70227d5c5559
viz_grid_scaled = minmaxnormalize(viz_grid)

# ╔═╡ c86af698-884d-46d2-8f29-6d0ce7e7d56e
u_pred = (
    reshape(trained_u(viz_grid_scaled), length(xs), length(ys), length(ts))
) * (max_pde_val - min_pde_val) .+ min_pde_val

# ╔═╡ dc8bbc94-9454-4336-a0ee-12085884c5f0
let fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y")
    errs = [abs.(u_pred[:,:,i] - u_real[:,:,i]) for i in eachindex(ts)]
    Colorbar(fig[1,2], limits = extrema(stack(errs)))

    CairoMakie.record(fig, "pinn_nested_ad.gif", eachindex(ts), framerate = 5) do i
        ax.title = @sprintf("Abs. Predictor Error | Time: %.2f", ts[i])
        err = errs[i]
        contour!(ax, xs, ys, err; levels = 10, linewidth = 2)
        heatmap!(ax, xs, ys, err)
        return fig
    end

    fig
end

# ╔═╡ 8f3b2acc-3440-4abb-a94e-ae14b6d0b05f
md"""
![](pinn_nested_ad.gif)
"""

# ╔═╡ 6d461a63-9f12-4282-b3d2-fdb352be2ecc
md"""
## Appendix
"""

# ╔═╡ d9b0b263-8633-44d9-b849-6878e32a187a
InteractiveUtils.versioninfo |> with_terminal

# ╔═╡ e35fa2b5-3e77-4976-b398-d22a0f62b4c3
with_terminal() do
    if @isdefined(MLDataDevices)
        if @isdefined(CUDA) && MLDataDevices.functional(CUDADevice)
            CUDA.versioninfo()
        end

        if @isdefined(AMDGPU) && MLDataDevices.functional(AMDGPUDevice)
            AMDGPU.versioninfo()
        end
    end
end

# ╔═╡ 5691fa54-6c11-43e9-89eb-420964841da0
Pkg.status |> with_terminal

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
# ╠═be4aade8-1f7f-46df-b006-d301c0835ac1
# ╠═222bd364-af5f-47ee-8120-bea7513ae40a
# ╠═d0c8d9dc-faf0-41ee-bc5e-e1696d118442
# ╠═cfbd31b8-8882-43c6-9c66-6cacd5d5668e
# ╠═4a73be53-0e19-44dd-ab43-5fced8f188e1
# ╠═7df36b77-65ea-4356-95a7-dd6b025780c4
# ╠═cf94c793-4fdd-4c31-8a69-dcc0786f98d7
# ╠═24fcec22-807d-4790-a4a2-9506dc48a044
# ╠═ed721e8b-eb50-4d36-b293-9ab7fc1cf31c
# ╠═46cf6020-33c8-4354-bd4e-164c8743c5e0
# ╠═a94c00e3-1ba9-4baf-a31f-fc1528aa788d
# ╠═1bbcbbc7-41bb-49f5-90f8-2f0552afec10
# ╠═c65994ca-0e8b-4c72-85df-c62fb7306461
# ╠═613cdd28-f832-4763-8cef-173c074b51d6
# ╠═576374ec-224b-4f47-b498-00d699e9e5f8
# ╠═cb0d5b82-eb2f-4698-99b8-6e3ee0cee6d9
# ╠═99b730a3-479e-4129-a642-3a71dc8700c2
# ╠═5ff62ed3-ab40-428f-9ba6-8be4ab1bf3da
# ╟─903fff49-d4a2-48d8-b689-2561e6152a97
# ╠═ddd1cbc7-7f1e-4a93-9f4f-07f5a4288a95
# ╠═4fe21a00-c083-4fd6-a695-c52240c3a060
# ╠═b67d7b58-c79a-45bc-a92e-81be68febf2d
# ╟─34efddda-66e0-4541-899c-beb7504a8362
# ╠═dea86dbf-0ffd-416a-88dd-b9300b792424
# ╠═1b3b4087-3283-40fc-9eee-62b719e54766
# ╠═24c21e66-f8dd-4327-ad65-e039673ab51b
# ╠═01694d01-3ad6-4076-878a-25651fc1d7f6
# ╠═ccf7acaf-64ed-46b9-9ed1-70227d5c5559
# ╠═c86af698-884d-46d2-8f29-6d0ce7e7d56e
# ╠═dc8bbc94-9454-4336-a0ee-12085884c5f0
# ╠═8f3b2acc-3440-4abb-a94e-ae14b6d0b05f
# ╟─6d461a63-9f12-4282-b3d2-fdb352be2ecc
# ╠═22b190dd-fba5-4e33-aad8-e4bba5d08073
# ╠═d2d252ab-d68a-4878-908a-ea13d9adf324
# ╠═d9b0b263-8633-44d9-b849-6878e32a187a
# ╠═e35fa2b5-3e77-4976-b398-d22a0f62b4c3
# ╠═5691fa54-6c11-43e9-89eb-420964841da0
