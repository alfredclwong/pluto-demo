### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 5266b77a-7779-11eb-331e-a5a301a484ad
begin
using PlutoUI: Slider
using Plots: plot, plot!, xlims!, ylims!, @layout

html"""
<h1>MLA Design Tool Mockup</h1>
<style>
main{
		margin-left:-15%;
		max-width:60%;
}
"""
end

# ╔═╡ 2c2ea610-796f-11eb-12d5-03770b36ae3a
module Wigner include("../Wigner.jl") end

# ╔═╡ 62eac400-795f-11eb-208f-87ee3dfc5356
let δ = .01,
	sliders = Dict(
		"device x" => (@bind xdevicewidth Slider(0:δ:28, default=28)),
		"eyebox x" => (@bind xeyeboxwidth Slider(0:δ:2, default=2)),
		"lens x" => (@bind xlenswidth Slider(.5:δ:3, default=3)),
		"fov x" => (@bind xfovmax Slider(0:δ:40, default=40)),
		"eye z" => (@bind zeye Slider(0:δ:20, default=20))
	)


md"""
$((pairs(sliders)...)...)
"""
end

# ╔═╡ 04e6e752-7839-11eb-181e-e5ecae56e424
begin
xdevice = (-xdevicewidth/2, xdevicewidth/2)
xeyebox = (-xeyeboxwidth/2, xeyeboxwidth/2)
xfov = (-xfovmax, xfovmax) .|> float

device = Wigner.box([xdevice, tand.(xfov)])
eyebox = Wigner.propagate(-zeye, Wigner.box([xeyebox, tand.(xfov)]))
lenses = [Wigner.cover(eyebox, (xi, xi+xlenswidth), xfov) for xi in Wigner.centredrange(xdevice..., xlenswidth)[1:end-1]]
filter!(l -> l != nothing, lenses)

nothing
end

# ╔═╡ e88f55e8-7a8c-11eb-32de-6b4ce50d5908
begin

elements = [device, eyebox, lenses...]

lplot = plot()
plot!.(elements)
xlims!(-32, 32)
ylims!(-1, 1)

rplot = plot()
plot!.(Wigner.propagate.(zeye, e) for e in elements)
xlims!(-32, 32)
ylims!(-1, 1)

plot(lplot, rplot, layout=@layout[a b], size=(852, 480))

end

# ╔═╡ Cell order:
# ╟─5266b77a-7779-11eb-331e-a5a301a484ad
# ╠═2c2ea610-796f-11eb-12d5-03770b36ae3a
# ╟─62eac400-795f-11eb-208f-87ee3dfc5356
# ╠═04e6e752-7839-11eb-181e-e5ecae56e424
# ╠═e88f55e8-7a8c-11eb-32de-6b4ce50d5908
