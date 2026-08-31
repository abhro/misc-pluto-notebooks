### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ f8722ee8-4a4e-41ee-8acb-b1c373eed434
using CairoMakie

# ╔═╡ 6555358a-0e7e-11f0-0a63-e716c46957ec
md"""
# Tracing magnetic field lines
"""

# ╔═╡ bf129681-b331-43f7-969a-a8b281e692e6
md"""
From <https://book.magneticearth.org/geomag-obs-models/03a_magnetic-field-line-tracing>.

Original authors: Ashley Smith, David Kerridge, Grace Cox, Alexandre Fournier, Chris Finlay, Alexander Grayver.
"""

# ╔═╡ 4bea5744-631e-4cf0-b3cf-0a0dde25bdca
md"""
## Axisymmetric multipole field lines
"""

# ╔═╡ 39516f5e-fdfb-4d51-b150-f063ce14d0ce
dipole(r₀, θ₀, θ) = r₀ * sin(θ)^2 / sin(θ₀)^2   # P(n=1,m=0)

# ╔═╡ bbcd300b-9340-4646-a747-bf19351c043d
quadrupole_helper(θ) = sin(θ)^2 * cos(θ)

# ╔═╡ 573b2a15-f5f0-413e-bf31-be3ab3eb8e9e
quadrupole(r₀, θ₀, θ) =                         # P(n=2,m=0)
    √(abs(r₀^2 / quadrupole_helper(θ₀)) *
      abs(quadrupole_helper(θ)))

# ╔═╡ 390b8de8-1737-40f8-b2c5-6fef3a36c92f
octupole_helper(θ) = sin(θ)^2 * (5 * cos(θ)^2 - 1)

# ╔═╡ e9f18041-946f-4094-b8a3-47a26cedcc08
octupole(r₀, θ₀, θ) =                           # P(n=3,m=0)
    ∛(
        abs(r₀^3 / octupole_helper(θ₀)) *
        abs(octupole_helper(θ))
    )

# ╔═╡ 56a14b35-aaeb-43ce-b5e1-f5564029d3de
hexadecapole_helper(θ) = (7 * cos(θ)^3 - 3cos(θ)) * sin(θ)^2

# ╔═╡ 2259e5a3-8794-42d6-ad20-47ef47e12ba8
function hexadecapole(r₀, θ₀, θ)                # P(n=4,m=0)
    k = r₀^4 / hexadecapole_helper(θ₀)
    P = hexadecapole_helper(θ)
    return (abs(k) * abs(P))^(1//4)
end

# ╔═╡ 390a438d-93ca-474b-b832-00b85207c056
md"""
## Plotting
"""

# ╔═╡ af24d4ae-3fe4-4933-afa3-42157912ca98
pole_type = 2

# ╔═╡ c2ff0217-722d-4947-9f39-1cc766ee70b9
theta = [5, 10, 15, 20, 30, 40, 45, 50];

# ╔═╡ 15f9b6b2-2339-45e4-abfb-bbc67f660745
axpoles = Dict(1 => dipole, 2 => quadrupole, 3 => octupole, 4 => hexadecapole);

# ╔═╡ dc0d7df2-b4f9-40e4-a7f9-ca1bce87964d
md"""
Define the Earth
"""

# ╔═╡ a1bf7eff-157c-446f-802d-04dcb4d4d153
r0 = 1;

# ╔═╡ f86a1c63-a1dc-474a-a945-bbd053f1ba3f
θe = range(0, π, length = 1000);

# ╔═╡ 619d040d-bb80-474d-8b9e-9b1e0cc4211b
clines = [:red, :blue, :grey, :purple, :brown, :purple, :pink, :orange, :magenta, :olive, :cyan];

# ╔═╡ 347ae72c-0d1a-4a65-b019-94957e5cdcf9
let fig = Figure()
    ax = Axis(
        fig[1, 1];
        xlabel = "Earth radii",
        ylabel = "Earth radii",
        xminorgridvisible = true,
        yminorgridvisible = true,
        xminorticksvisible = true,
        yminorticksvisible = true,
    )

    poly!(Circle(Point2f(0, 0), r0), color = :lightgrey) # Plot the Earth

    ic = -1
    for i in theta
        ic += 1
        θ₀ = deg2rad(i)
        θ = range(θ₀, π, length = 1000)
        rad = axpoles[pole_type].(r0, θ₀, θ)
        xb = rad .* sin.(θ)
        yb = rad .* cos.(θ)
        inside_earth_idx = rad .< 1
        xb[inside_earth_idx] .= NaN
        yb[inside_earth_idx] .= NaN
        color = clines[ic % 10 + 1]
        lines!(ax, xb, yb; color)
        # Assume a symmetrical distribution in the four quadrants
        lines!(ax, xb, -yb; color)
        lines!(ax, -xb, yb; color)
        lines!(ax, -xb, -yb; color)
    end

    fig
end

# ╔═╡ 7b106b3f-93d2-4e65-8bf5-090906423470
md"""
## Notebook setup 🔧
"""

# ╔═╡ 89bde001-5138-4818-87fe-8f4bb276614b
import PlutoUI: TableOfContents

# ╔═╡ 56732420-62f7-4b43-af9e-e843774c5cf4
TableOfContents()

# ╔═╡ Cell order:
# ╟─6555358a-0e7e-11f0-0a63-e716c46957ec
# ╟─bf129681-b331-43f7-969a-a8b281e692e6
# ╟─4bea5744-631e-4cf0-b3cf-0a0dde25bdca
# ╠═39516f5e-fdfb-4d51-b150-f063ce14d0ce
# ╠═573b2a15-f5f0-413e-bf31-be3ab3eb8e9e
# ╠═bbcd300b-9340-4646-a747-bf19351c043d
# ╠═e9f18041-946f-4094-b8a3-47a26cedcc08
# ╠═390b8de8-1737-40f8-b2c5-6fef3a36c92f
# ╠═2259e5a3-8794-42d6-ad20-47ef47e12ba8
# ╠═56a14b35-aaeb-43ce-b5e1-f5564029d3de
# ╟─390a438d-93ca-474b-b832-00b85207c056
# ╠═af24d4ae-3fe4-4933-afa3-42157912ca98
# ╠═c2ff0217-722d-4947-9f39-1cc766ee70b9
# ╠═15f9b6b2-2339-45e4-abfb-bbc67f660745
# ╟─dc0d7df2-b4f9-40e4-a7f9-ca1bce87964d
# ╠═a1bf7eff-157c-446f-802d-04dcb4d4d153
# ╠═f86a1c63-a1dc-474a-a945-bbd053f1ba3f
# ╠═619d040d-bb80-474d-8b9e-9b1e0cc4211b
# ╠═f8722ee8-4a4e-41ee-8acb-b1c373eed434
# ╠═347ae72c-0d1a-4a65-b019-94957e5cdcf9
# ╟─7b106b3f-93d2-4e65-8bf5-090906423470
# ╠═89bde001-5138-4818-87fe-8f4bb276614b
# ╠═56732420-62f7-4b43-af9e-e843774c5cf4
