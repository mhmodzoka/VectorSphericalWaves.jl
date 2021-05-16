# Importing
using LinearAlgebra
using SpecialFunctions
using VectorSphericalWaves
using Zygote
using ChainRulesCore


# calculate B(θ), C(θ), P(θ)
function BB_mn_of_θ(m, n, θ)
    """
    I assume each of m, n, θ is a single number
    """
    # TODO you can use literal syntax for this
    return vcat(
        0,                  # r-component
        VectorSphericalWaves.τₘₙ(m, n, θ),      # θ-component
        im .* VectorSphericalWaves.πₘₙ(m, n, θ)  # ϕ-component      
    ) # equation C.19
end
# jacobian(BB_mn_of_θ, 1, 1, 1.2)

function BB_cross_BB_mn_of_θ(m, n, θ)   
    # TODO you can use literal syntax for this
    z = cross(BB_mn_of_θ(m, n, θ), 2 .* im .* BB_mn_of_θ(m, n, θ))
    return vcat(real(z), imag(z))
end
# jacobian(BB_cross_BB_mn_of_θ, 1, 1, 1.2)

function B_mn_of_θ_real(m, n, θ)
    x = real(VectorSphericalWaves.B_mn_of_θ(m, n, θ)) .+ imag(VectorSphericalWaves.B_mn_of_θ(m, n, θ))
    return x
end
jacobian(B_mn_of_θ_real, 1, 1, 1.2)



function frule(
    (NO_FIELDS, Δm, Δn, Δθ),
    BB_mn_of_θ,
    m, n, θ
)
    y = BB_mn_of_θ(m, n, θ)
    Δy = vcat(
        (0, 0, 0,),                  # r-component
        gradient(VectorSphericalWaves.τₘₙ, m, n, θ),      # θ-component
        im .* gradient(VectorSphericalWaves.πₘₙ, m, n, θ)  # ϕ-component      
    )
    print("Δy")
    print(Δy)
    return y, Δy
end

function BB_cross_BB_mn_of_θ(m, n, θ)
    return cross(BB_mn_of_θ(m, n, θ), BB_mn_of_θ(m, n, θ + 1e-3))
end

function BB_dot_BB_mn_of_θ(m, n, θ)
    return dot(VectorSphericalWaves.B_mn_of_θ(m, n, θ), VectorSphericalWaves.B_mn_of_θ(m, n, θ + 1e-3))
end

# trying to auto differentiate from lower to upper level
gradient(besselj, 1, 1.2)
gradient(VectorSphericalWaves.wignerdjmn, 1, 1, 1, 1.2)
gradient(VectorSphericalWaves.πₘₙ, 1, 1, 1.2)
gradient(VectorSphericalWaves.τₘₙ, 1, 1, 1.2)
jacobian(BB_mn_of_θ, 1, 1, 1.2)
jacobian(BB_cross_BB_mn_of_θ, 1, 1, 1.2)
gradient(BB_dot_BB_mn_of_θ, 1, 1, 1.2)
jacobian(VectorSphericalWaves.B_mn_of_θ, 1, 1, 1.2) # ERROR: ArgumentError: jacobian does not accept complex output

# test Jacobiam
m, n, θ = 1, 1, 0.2
jacobian(BB_mn_of_θ, m, n, θ)
small = 1e-10
(BB_mn_of_θ(m, n, θ + small) - BB_mn_of_θ(m, n, θ - small)) / (2 * small)

gradient(VectorSphericalWaves.wignerdjmn_ELZOUKA, 1, 1, 1, 0.1)