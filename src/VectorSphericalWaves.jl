module VectorSphericalWaves

# Vector spherical harmonics calculated using equations from:
# Mishchenko, M.I., Travis, L.D., and Lacis, A.A. (2002). Scattering, absorption, and emission of light by small particles (Cambridge University Press).

# returned vector spherical wave functions (VSWF) are matrices have the following indices:
# - index for component direction (e.g., 0,1,2 represent r,θ,ϕ)
# - index for order and rank (mn)
# - index for point in space

# inputs
# -- r,θ,ϕ : each is a 1D array, representing spactial spehrical coordinates of points in space, where VSWF will be evaluated
# -- n_max : integer, maximum rank of VSWF
# --

include("utils.jl")
include("vectorSphericalWaves_complex.jl")
include("vectorSphericalWaves_complex_all_indices.jl")
include("vectorSphericalWaves_separateRealImag.jl")

end # module
