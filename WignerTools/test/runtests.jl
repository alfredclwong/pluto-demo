using Test
using WignerTools
using StaticArrays

@test propagate(0.1, (@SVector [0, 1]); θ̄ = (@SVector [0, 45])) == [0, 1.1]
