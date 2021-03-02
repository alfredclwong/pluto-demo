using Test
using WignerTools  # warning: need to ] dev WignerTools before running tests
using StaticArrays

@test propagate(0.1, (@SVector [0, 1]); θ̄ = (@SVector [0, 45])) == [0, 1.1]
