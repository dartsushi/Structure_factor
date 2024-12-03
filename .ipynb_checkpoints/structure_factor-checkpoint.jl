using MPSKit, TensorKit, KrylovKit
T = ComplexF64
#=
    This is an implementation of the structure factor calculation based on PHYSICAL REVIEW B 105, 195140 (2022)
=#

χ = 8
ψtop = InfiniteMPS(rand,T,ℂ^2,ℂ^χ)
