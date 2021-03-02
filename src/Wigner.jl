using StaticArrays: SVector, SMatrix, @SVector, @SMatrix
using Colors: HSVA
import Polyhedra
using Polyhedra: polyhedron, vrep, MixedMatVRep, Polyhedron, HalfSpace, vrepiscomputed

"""
	propagate(dz::Real, x̄::SVector{N,T}; θ̄::SVector{N,<:Real}) where {N,T<:Real}

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

# """
# 	Ray{N,T<:Real}

# `N` is the dimension of the hyperplane within which rays are sourced.
# `T` is the numeric type used to represent positions and angles.

# Angles (degrees) are taken as deviations from the generalized z-axis within their
# respective dimensions (field of view).
# """
# struct Ray{N,T<:Real}
# 	x̄::SVector{N,T}
# 	θ̄::SVector{N,T}
# end

# """
# 	Ray(x̄::SVector{N,T}, θ̄::SVector{N,T}, z::Real) where {N,T<:Real}

# Find a ray's source from its propagated position (e.g. on the eyebox).
# """
# function Ray(x̄::SVector{N,T}, θ̄::SVector{N,T}, z::Real) where {N,T<:Real}
# 	Ray(propagate(-z, x̄; θ̄), θ̄)
# end

# """
# 	propagate(z::Real, ray::Ray)

# Propagate a Ray.
# """
# propagate(z::Real, ray::Ray) = Ray(ray.x̄, ray.θ̄, -z)

# """
# 	point(ray::Ray; tan::Bool = false)

# Express a ray as a 2N-dim vector. If `tan` is true, take tand.(θ̄) instead of θ̄.
# """
# function point(ray::Ray{N}; tan::Bool = false) where N
# 	SVector{2N}(ray.x̄..., (tan ? tand : identity).(ray.θ̄)...)
# end

## POLYHEDRA ##

propagate(dz::Real, rays::Polyhedron) = (@SMatrix [1 dz;0 1]) * rays

"""
	e(i::Integer, n::Integer)
"""
function e(i::Integer, n::Integer)
	ē = zeros(n)
	ē[i] = 1
	SVector{n}(ē)
end

"""
	box(intervals::SMatrix{N,2}) where N

N-dim axis-aligned hyperrectangle.
"""
function box(intervals::SMatrix{N,2}) where N
	halfspaces = reduce(vcat,
		[HalfSpace(s(e(i,N)), s(bound)) for (s, bound) in zip((-,+), intervals[i,:])]
		for i in 1:N
	)
	polyhedron(reduce(∩, halfspaces))
end

"""
	box(intervals::Array{<:Real,2})

Array convenience conversion.
"""
box(intervals::Array{<:Real,2}) = box(SMatrix{size(intervals)...}(intervals))

"""
	box(intervals::Vector{<:Tuple{Real,Real}})

Vector{Tuple} convenience conversion.
"""
function box(intervals::Vector{<:Tuple{Real,Real}})
	box(reduce(vcat, [[interval[1] interval[2]] for interval in intervals]))
end

function boundingbox(poly::Polyhedron)
	vertexlist = vertices(poly)
	if isempty(vertexlist)
		return nothing
	end
	mins, maxs = (f(vertexlist, dims=1) for f in [minimum, maximum])
	box(reduce(vcat, zip(mins, maxs)))
end

function cover(eyebox::Polyhedron, xlens::Tuple{<:Real,<:Real}, xfov::Tuple{<:Real,<:Real})
    fulllens = box([xlens, xfov])
	intersection = intersect(eyebox, fulllens)
	boundingbox(intersection)
end

vertices(poly::Polyhedron) = MixedMatVRep(vrep(poly)).V

# function poly(rays::Vector{<:Ray}; tan::Bool = false)
# 	f = tan ? tand : identity
# 	polyhedron(vrep([[p[1], f(p[2])] for p in point.(rays)]))
# end

## OTHER ##

function colorrange(n::Integer; s::Real=1, v::Real=1, a::Real=1)
	n == 0 ? [] : range(HSVA(0,s,v,a); stop=HSVA(359,s,v,a), length=n+1)[1:n]
end

function centredrange(start::Real, stop::Real, step::Real)
	range = stop - start
	n = max(0, floor(range / step))
	start+(range-n*step)/2:step:stop
end
