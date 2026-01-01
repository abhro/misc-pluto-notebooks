### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 6b4b4e10-fce3-4888-8326-3e3eebb42ff0
using Lux, LuxCUDA

# ╔═╡ 75e7943b-6b1e-4f81-a5c2-00ccf126bda7
using ADTypes, MLUtils, Optimisers, Random

# ╔═╡ d25a81f7-2718-4bfd-8a53-74957445ca30
using Printf

# ╔═╡ 63820f11-90d4-4124-858f-410c2d0a9121
using JLD2

# ╔═╡ 3a3224dc-e90c-498c-bee0-ba26e4289667
using Zygote

# ╔═╡ bcaa13fe-f331-44a6-9dea-7963756a5e7c
using PlutoUI

# ╔═╡ eaad4dc0-67f5-11f0-0d39-bf6ac3456cff
md"""
# Training a Simple LSTM
"""

# ╔═╡ 3a8dfd7b-df83-401f-b8a9-11531892e8dd
md"""
LSTM: Long Short-Term Memory. An extension of Recurrent Neural Networks (RNNs).
"""

# ╔═╡ 02bb568c-48cc-457c-9ee6-1920a13a9e18
md"""
Original tutorial from <https://lux.csail.mit.edu/dev/tutorials/beginner/3_SimpleRNN>.
"""

# ╔═╡ 405f1ff7-e662-4958-ace4-b2d21b424c95
TableOfContents()

# ╔═╡ b90cb8ab-2dba-4e1c-845b-6cc3ccb96155
md"""
## Dataset
"""

# ╔═╡ 09a1ebaf-4566-4c54-abb4-e5023e5dffc6
function create_dataset(; dataset_size = 1000, sequence_length = 50)
    # create the spirals
    data = [MLUtils.Datasets.make_spiral(sequence_length) for _ in 1:dataset_size]
    # Get the labels
    labels = vcat(fill(0.0f0, dataset_size÷2), fill(1.0f0, dataset_size÷2))

    clockwise_spirals = [
        reshape(d[1][:, 1:sequence_length], :, sequence_length, 1)
        for d in data[1:(dataset_size÷2)]
    ]
    anticlockwise_spirals = [
        reshape(d[1][:, (sequence_length+1):end], :, sequence_length, 1)
        for d in data[(dataset_size÷2 + 1):end]
    ]

    x_data = Float32.(cat(clockwise_spirals..., anticlockwise_spirals...; dims=3))

    return x_data, labels
end

# ╔═╡ 0acda80a-e3b1-4a49-936d-f861df0a6469
function get_dataloaders(features, labels)
    # Split the dataset
    (Xtrain, ytrain), (Xval, yval) = splitobs((features, labels); at = 0.8, shuffle=true)

    # Create DataLoaders

    # Use DataLoader to automatically minibatch and shuffle the data
    trainloader = DataLoader(
        collect.((Xtrain, ytrain)); batchsize = 128, shuffle = true, partial = false)
    # Don't shuffle the validation data
    valloader = DataLoader(
        collect.((Xval, yval)); batchsize = 128, shuffle = false,  partial = false)

    return (trainloader, valloader)
end

# ╔═╡ ca07e5a0-b325-4be8-ba60-7b1e411dd550
md"""
## Creating a Classifier
"""

# ╔═╡ fb2405a9-9ea5-42d3-936f-ad0388e46932
@doc AbstractLuxContainerLayer

# ╔═╡ 8c3d9196-44f6-4421-ad41-4b20fdf152f1
begin
    struct SpiralClassifier{L, C} <: AbstractLuxContainerLayer{(:lstm_cell, :classifier)}
        lstm_cell::L
        classifier::C
    end

    function SpiralClassifier(in_dims, hidden_dims, out_dims)
        return SpiralClassifier(
            LSTMCell(in_dims => hidden_dims),
            Dense(hidden_dims => out_dims, sigmoid),
        )
    end

    function (sc::SpiralClassifier)(x::AbstractArray{T,3}, ps::NamedTuple, st::NamedTuple) where T

        # First we will have to run the sequence through the LSTM Cell.
        # The first call to LSTM Cell will create the initial hidden state.
        # See that the parameters and states are automatically populated into a field called `lstm_cell`.
        # We use `eachslice()` t et the elements in the sequence without copying, and
        # `Iterators.peel` to split out the first element for LSTM initialization.
        X_init, X_rest = Iterators.peel(LuxOps.eachslice(x, Val(2)))
        (y, carry), st_lstm = sc.lstm_cell(X_init, ps.lstm_cell, st.lstm_cell)

        # Now that we have the hidden state and memory in `carry` we will pass the input and `carry` jointly
        for x in X_rest
            (y, cary), st_lstm = sc.lstm_cell((x, carry), ps.lstm_cell, st_lstm)
        end

        # After rnning through the sequence we will pass the output through the classifier
        y, st_classifier = sc.classifier(y, ps.classifier, st.classifier)
        # Finally remember to create the updated state
        st = merge(st, (; classifier = st_classifier, lstm_cell = st_lstm))

        return vec(y), st
    end
end

# ╔═╡ 200c0b5c-ef80-4a74-9069-8abc930abf4c
md"""
## Using the `@compact` API
"""

# ╔═╡ 25f19c69-991f-4286-8aad-8e85a47b2a5e
function SpiralClassifierCompact(in_dims, hidden_dims, out_dims)
    lstm_cell = LSTMCell(in_dims => hidden_dims)
    classifier = Dense(hidden_dims => out_dims, sigmoid)
    return @compact(; lstm_cell, classifier) do x::AbstractArray{T,3} where T
        x_init, x_rest = Iterators.peel(LuxOps.eachslice(x, Val(2)))
        y, carry = lstm_cell(x_init)
        for x in x_rest
            y, carry = lstm_cell((x, carry))
        end
        @return vec(classifier(y))
    end
end

# ╔═╡ 36c5e673-4ba6-41e1-84ec-699d568191f3
SpiralClassifier(1, 2, 3)

# ╔═╡ cce42bf0-6242-4252-abbd-1863e0892fa0
SpiralClassifierCompact(1, 2, 3)

# ╔═╡ 36c5f77d-7717-44a8-8399-5946c03fd92a
md"""
## Defining Accuracy and Loss
"""

# ╔═╡ bcfad818-c154-4939-9b74-311c615c87bb
const lossfn = BinaryCrossEntropyLoss()

# ╔═╡ 4abe4041-9ccb-4674-ad24-edb5d4867f33
function compute_loss(model, ps, st, (X, y))
    ŷ, st_post = model(x, ps, st)
    loss = lossfn(ŷ, y)
    return loss, st_post, (; y_pred = ŷ)
end

# ╔═╡ 784dcb44-5ffc-4207-b0c1-10564f9264c4
matches(ŷ, y) = sum((ŷ .> 0.5f0) .== y)

# ╔═╡ 3d3af6ba-5513-414e-a83d-a59cffd0eddf
accuracy(ŷ, y) = matches(ŷ, y) / length(ŷ)

# ╔═╡ 28d9a055-0847-4584-ad4a-e75282ab7d58
cdev = cpu_device();

# ╔═╡ 8ae75e1c-36d3-4d3f-a8f5-9ba36d54914a
gdev = gpu_device();

# ╔═╡ aaaf486c-f2fb-47ce-8e40-ef2e0761e967
model = SpiralClassifier(2, 8, 1)

# ╔═╡ d7c3c72e-b1cb-4f32-bafc-5d64da678886
function train_one_epoch!(ad, lossfn, dataloader, train_state)
    total_loss = 0.0f0
    total_samples = 0
    for (X, y) in dataloader
        (_, loss, _, train_state) = Training.single_train_step!(ad, lossfn, (X, y), train_state)
        total_loss += loss * length(y)
        total_samples += length(y)
    end

    return train_state, total_loss, total_samples
end

# ╔═╡ 2fa93fd9-b0b4-4649-8128-1e953b74d3de
function validate_one_epoch(model, dataloader, train_state)
    total_acc = 0.0f0
    total_loss = 0.0f0
    total_samples = 0

    st_frozen = Lux.testmode(train_state.states)
    for (X, y) in dataloader
        ŷ, st_ = model(X, train_state.parameters, st_frozen)
        ŷ, y = cdev(ŷ), cdev(y)
        total_acc += accuracy(ŷ, y) * length(y)
        total_loss += lossfn(ŷ, y) * length(y)
        total_samples += length(y)
    end

    return train_state, total_loss, total_samples, total_acc
end

# ╔═╡ 7bf4a845-566b-4269-9a69-bd5df474301a
function train(model, dataloaders; nepochs=25)
    ps_init, st_init = Lux.setup(Xoshiro(1), model) |> gdev

    train_loader, val_loader = dataloaders |> gdev

    train_state = Training.TrainState(model, ps_init, st_init, Adam(0.01f0))

    ad = AutoZygote()

    for epoch in 1:nepochs
        # Train the model
        (train_state, total_loss, total_samples) = train_one_epoch!(
            ad, lossfn, train_loader, train_state)
        @printf("Epoch [%3d]: Loss %4.5f\n", epoch, total_loss/total_samples)

        # Validate the model
        train_state, total_loss, total_samples, total_acc = validate_one_epoch(
                model, val_loader, train_state)
        @printf("Validation:\tLoss %4.5f\tAccuracy %4.5f\n",
                total_loss/total_samples, total_acc/total_samples)
    end

    return (train_state.parameters, train_state.states) |> cdev
end

# ╔═╡ f663e6ee-da89-42d6-913d-59ee6ea98a9a
ps_trained, st_trained = train(model, get_dataloaders(create_dataset()...))

# ╔═╡ b9e6d8b8-2d2e-4040-9c64-bcb56e2f28a2
ps_trained2, st_trained2 = train(SpiralClassifierCompact(2, 8, 1), get_dataloaders(create_dataset()...))

# ╔═╡ bcb423a7-2f43-4dbf-a21e-a3592c8685bd
md"""
## Saving the model
"""

# ╔═╡ Cell order:
# ╟─eaad4dc0-67f5-11f0-0d39-bf6ac3456cff
# ╟─3a8dfd7b-df83-401f-b8a9-11531892e8dd
# ╟─02bb568c-48cc-457c-9ee6-1920a13a9e18
# ╠═6b4b4e10-fce3-4888-8326-3e3eebb42ff0
# ╠═75e7943b-6b1e-4f81-a5c2-00ccf126bda7
# ╠═d25a81f7-2718-4bfd-8a53-74957445ca30
# ╠═63820f11-90d4-4124-858f-410c2d0a9121
# ╠═3a3224dc-e90c-498c-bee0-ba26e4289667
# ╠═bcaa13fe-f331-44a6-9dea-7963756a5e7c
# ╠═405f1ff7-e662-4958-ace4-b2d21b424c95
# ╟─b90cb8ab-2dba-4e1c-845b-6cc3ccb96155
# ╠═09a1ebaf-4566-4c54-abb4-e5023e5dffc6
# ╠═0acda80a-e3b1-4a49-936d-f861df0a6469
# ╟─ca07e5a0-b325-4be8-ba60-7b1e411dd550
# ╠═fb2405a9-9ea5-42d3-936f-ad0388e46932
# ╠═8c3d9196-44f6-4421-ad41-4b20fdf152f1
# ╟─200c0b5c-ef80-4a74-9069-8abc930abf4c
# ╠═25f19c69-991f-4286-8aad-8e85a47b2a5e
# ╠═36c5e673-4ba6-41e1-84ec-699d568191f3
# ╠═cce42bf0-6242-4252-abbd-1863e0892fa0
# ╟─36c5f77d-7717-44a8-8399-5946c03fd92a
# ╠═bcfad818-c154-4939-9b74-311c615c87bb
# ╠═4abe4041-9ccb-4674-ad24-edb5d4867f33
# ╠═784dcb44-5ffc-4207-b0c1-10564f9264c4
# ╠═3d3af6ba-5513-414e-a83d-a59cffd0eddf
# ╠═28d9a055-0847-4584-ad4a-e75282ab7d58
# ╠═8ae75e1c-36d3-4d3f-a8f5-9ba36d54914a
# ╠═aaaf486c-f2fb-47ce-8e40-ef2e0761e967
# ╠═7bf4a845-566b-4269-9a69-bd5df474301a
# ╠═d7c3c72e-b1cb-4f32-bafc-5d64da678886
# ╠═2fa93fd9-b0b4-4649-8128-1e953b74d3de
# ╠═f663e6ee-da89-42d6-913d-59ee6ea98a9a
# ╠═b9e6d8b8-2d2e-4040-9c64-bcb56e2f28a2
# ╟─bcb423a7-2f43-4dbf-a21e-a3592c8685bd
