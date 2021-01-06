### A Pluto.jl notebook ###
# v0.12.18

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

# ╔═╡ 8b0a2bd0-4ffa-11eb-36d8-1d012eacd12a
using Plots, LinearAlgebra#, StructArrays

# ╔═╡ f9940280-4f8a-11eb-3271-59c3c3601353
md"""
# Neuroscience
Ported from some coursework I did in Python/``\LaTeX`` last year.

#### Physiological model of a spiking neuron
##### Background
Neuron spikes are modelled as 4D systems ``\mathbf{x} = (v, m, h, n)`` where v denotes the membrane potential (in mV) and ``n,\ m,\ h \\\in [0, 1]`` are dimensionless state variables representing ion channel activation. The Hodgkin-Huxley model is as follows

``\dot{v} = -g_{N_a} m^3 h (v-e_{N_a}) - g_K n^4 (v-e_K) - g_L (v-e_L) + I_{ext}``

``\dot{m} = \alpha_m(v) (1-m) - \beta_m(v) m``

``\dot{h} = \alpha_h(v) (1-h) - \beta_h(v) h``

``\dot{n} = \alpha_n(v) (1-n) - \beta_n(v) n``
"""

# ╔═╡ 842eccd2-4ffa-11eb-25cd-4f80f90f6d38
α(v::T) where {T <: AbstractFloat} = [
	(2.5-.1(v+65)) / (exp(2.5-.1(v+65))-1);  # α_m
	(.07exp(-(v+65)/20));                    # α_h
	(.1-.01(v+65)) / (exp(1-.1(v+65))-1)]    # α_n

# ╔═╡ 57507b90-4ffb-11eb-0a6f-6357944096f6
β(v::T) where {T <: AbstractFloat} = [
	4exp(-(v+65)/18);                        # β_m
	1/(exp(3-.1(v+65))+1);                   # β_h
	.125exp(-(v+65)/80)]                     # β_n

# ╔═╡ 17690c60-4ffe-11eb-031a-ab50ed1623e9
#=
begin
	local vs = collect(-65.0:.1:0.0)
	plot(vs, vcat(map(f -> hcat(map(f, vs)...), [α, β])...)',
		label=hcat(("$(a)_$(b)" for a in split("αβ","") for b in split("mhn",""))...),
		size=(2500, 1000), thickness_scaling=4, legend=:outertopright)
end
=#

# ╔═╡ b8bd00be-501d-11eb-0204-1998b4ca14df
begin
	local vs = collect(-100.0:.1:50.0)
	plot(vs, hcat([α(v)./(α(v).+β(v)) for v in vs]...)',
		xlim=(vs[1], vs[end]), ylim=(0, 1),
		label=hcat(["$(c)_∞" for c in split("mhn", "")]...),
		size=(2500, 1000), thickness_scaling=4, legend=:outertopright)
	plot!([-65], seriestype=:vline, label=nothing, line=(:dash, :black))
end

# ╔═╡ d0150aa0-5014-11eb-0a0a-9b5e244bd725
md"""
We simulate a step input applied to a neuron at rest (approx. -65mV). The slider below allows us to adjust the amplitude of the input current.
"""

# ╔═╡ df675930-4fb1-11eb-02f4-41b64aa8084d
struct State{T}
	v::T
	m::T
	h::T
	n::T
end

# ╔═╡ 5f2ce9c0-4fa1-11eb-2fea-4ffd865ff528
function tick(x::State, I_ext::T1; δ::T2) where {T1 <: AbstractFloat, T2 <: AbstractFloat}
	g = [120; 36; .3]
	e = [50; -77; -54.4]
	mhn = [x.m; x.h; x.n]
	nlin = [x.m^3 * x.h; x.n^4; 1]

	vdot = -g .* nlin ⋅ (x.v .- e) + I_ext
	mhndot = (α(x.v) .* (1 .- mhn)) .- (β(x.v) .* mhn)

	State(x.v + vdot*δ, (mhn .+ mhndot*δ)...)
end

# ╔═╡ bbca554e-4fb0-11eb-0a77-af54bd299739
function constant_pulse(I_ext::T1, t1::T2, t0::T2=50.; δ::T2=1e-3) where {T1 <: AbstractFloat, T2 <: AbstractFloat}
	n0 = Int(t0/δ); n1 = Int(t1/δ); n = n0 + n1
	I = vcat(zeros(n0), ones(n1)*I_ext)
	#X = StructArray(x for _ in 1:n)
	X = [State(0., 0., 0., 0.) for _ in 1:n]
	X[1] = State(-65., .05, .60, .32)
	for i = 2:n
		X[i] = tick(X[i-1], I[i-1], δ=δ)
	end
	(0:δ:t1, X[n0:end])
end

# ╔═╡ c7d9db40-5019-11eb-377f-499b5fe6e762
@bind I_ext html"<input type=\"range\" value=\"2.3\" min=\"0\" max=\"10\" step=\".02\" style=\"width:30%\" oninput=\"this.nextElementSibling.value=this.value\"><output>2.3</output>"

# ╔═╡ 47c69dae-501b-11eb-0a7e-9bcf8f023cf1
@bind t html"<input type=\"number\" value=\"100\">"

# ╔═╡ baaf2c70-5020-11eb-0716-0512c15562c7
begin
	local T, X = constant_pulse(Float64(I_ext), Float64(t), δ=1e-2)
	local l = @layout [a ; b]
	p1 = plot(T, [x.v for x in X],
		xlims=(0, t), ylims=(-100, 50),
		size=(2000, 2000), thickness_scaling=4, legend=false)
	p2 = plot(T, hcat([[x.m, x.h, x.n] for x in X]...)',
		xlims=(0, t),
		size=(2000, 2000), thickness_scaling=4, legend=false)
	plot(p1, p2, layout=l)
end

# ╔═╡ Cell order:
# ╠═8b0a2bd0-4ffa-11eb-36d8-1d012eacd12a
# ╟─f9940280-4f8a-11eb-3271-59c3c3601353
# ╟─842eccd2-4ffa-11eb-25cd-4f80f90f6d38
# ╟─57507b90-4ffb-11eb-0a6f-6357944096f6
# ╟─17690c60-4ffe-11eb-031a-ab50ed1623e9
# ╟─b8bd00be-501d-11eb-0204-1998b4ca14df
# ╟─d0150aa0-5014-11eb-0a0a-9b5e244bd725
# ╠═df675930-4fb1-11eb-02f4-41b64aa8084d
# ╟─5f2ce9c0-4fa1-11eb-2fea-4ffd865ff528
# ╟─bbca554e-4fb0-11eb-0a77-af54bd299739
# ╟─c7d9db40-5019-11eb-377f-499b5fe6e762
# ╠═47c69dae-501b-11eb-0a7e-9bcf8f023cf1
# ╠═baaf2c70-5020-11eb-0716-0512c15562c7
