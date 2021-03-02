module WignerTools
export propagate, box, cover, centredrange

using StaticArrays: SVector, SMatrix, @SVector, @SMatrix
using Colors: HSVA
import Polyhedra
using Polyhedra: polyhedron, vrep, MixedMatVRep, Polyhedron, HalfSpace, vrepiscomputed, DefaultLibrary
import GLPK
lib = DefaultLibrary{Float64}(GLPK.Optimizer)

"""
	propagate(dz::Real, x̄::SVector{N,T}; θ̄::SVector{N,<:Real}, tanθ̄::SVector{N,<:Real}) where {N,T<:Real}

Propagate a ray along the z-axis by `dz`.

`x̄` gives the initial position in the hyperplane.
`θ̄` specifies the angle of the ray from the z-axis within each dimension in degrees.
`tanθ̄` is a convenience option for when tan(θ̄) has already been computed.

Either `θ̄` or `tanθ̄` should be provided, but not both.
"""
function propagate(
	dz::Real,
	x̄::SVector{N,<:Real};
	θ̄::SVector{N,<:Real} = (@SVector zeros(N)),
	tanθ̄::SVector{N,<:Real} = (@SVector zeros(N))
) where N
	@assert θ̄ .* tanθ̄ == zeros(N)
	x̄ + (tand.(θ̄) + tanθ̄) * dz
end

## POLYHEDRA ##

"""
	propagate(dz::Real, rays::Polyhedron)

A polyhedron represents a collection of rays. Assuming angles are given as tanθ̄, we apply
a shear matrix to propagate the entire collection efficiently.
"""
propagate(dz::Real, rays::Polyhedron) = (@SMatrix [1 dz;0 1]) * rays

"""
	vertices(poly::Polyhedron)

Force computation of, and extract, vrep from polyhedron.
"""
vertices(poly::Polyhedron) = MixedMatVRep(vrep(poly)).V

"""
	box(intervals::SMatrix{N,2}) where N

N-dim axis-aligned hyperrectangle.
"""
function box(intervals::SMatrix{N,2}) where N
	halfspaces = reduce(vcat,
		[HalfSpace(s(e(i,N)), s(bound)) for (s, bound) in zip((-,+), intervals[i,:])]
		for i in 1:N
	)
	polyhedron(reduce(∩, halfspaces), lib)
end

box(intervals::Array{<:Real,2}) = box(SMatrix{size(intervals)...}(intervals))

function box(intervals::Vector{<:Tuple{Real,Real}})
	box(reduce(vcat, [[interval[1] interval[2]] for interval in intervals]))
end

"""
	boundingbox(poly::Polyhedron)

N-dim axis-oriented bounding box for a polyhedron. Returns nothing if polyhedron is empty.
"""
function boundingbox(poly::Polyhedron)
	vertexlist = vertices(poly)
	if isempty(vertexlist)
		return nothing
	end
	mins, maxs = (f(vertexlist, dims=1) for f in [minimum, maximum])
	box(reduce(vcat, zip(mins, maxs)))
end

"""
	cover(eyebox::Polyhedron, xlens::Tuple{<:Real,<:Real}, xfov::Tuple{<:Real,<:Real})

Cover the eyebox across a given range in the device plane and fov.
"""
function cover(eyebox::Polyhedron, xlens::Tuple{<:Real,<:Real}, xfov::Tuple{<:Real,<:Real})
    fulllens = box([xlens, xfov])
	intersection = intersect(eyebox, fulllens)
	boundingbox(intersection)
end

## OTHER ##

"""
	e(i::Integer, n::Integer)

Canonical basis vector utility function.
"""
function e(i::Integer, n::Integer)
	ē = zeros(n)
	ē[i] = 1
	SVector{n}(ē)
end

function colorrange(n::Integer; s::Real=1, v::Real=1, a::Real=1)
	n == 0 ? [] : range(HSVA(0,s,v,a); stop=HSVA(359,s,v,a), length=n+1)[1:n]
end

function centredrange(start::Real, stop::Real, step::Real)
	range = stop - start
	n = max(0, floor(range / step))
	start+(range-n*step)/2:step:stop
end

end # module
