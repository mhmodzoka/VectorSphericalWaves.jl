using SpecialFunctions
using DataStructures
using StaticArrays
using ComplexOperations

import WignerD
import ChainRulesCore


export γ_mn
export γ_mn_dash
export wignerdjmn_ELZOUKA
export ∂wignerdjmn_by_∂θ
export single_index_from_m_n
export get_max_single_index_from_n_max


#############################################################################################
# import "WignerD.jl" code to calculate wigner-d (https://github.com/jishnub/WignerD.jl)
# [TODO] I have already installed "WignerD" package using the following commands: # include("WignerD.jl/src/WignerD.jl") ;  import Pkg; Pkg.add("WignerD")
# include("WignerD.jl/src/WignerD.jl")
# using .WignerD

#############################################################################################
"""
    Dirac delta function
returns 1 if the two inputs are equal. Otherwise, return Zero
"""
δ(x,y) = ==(x, y)

"""
    Fuse the two indices (m,n) into a single index.
"""
function single_index_from_m_n(m, n)
    return n * (n + 1) + m
end

function get_max_single_index_from_n_max(n_max)
    return single_index_from_m_n(n_max, n_max)
end


#############################################################
"""
    defining derivative for Bessel functions.
Inside "SpecialFunctions.jl" package, the definition of derivative of Bessel functions exists. However, it returns an error when called by Zygote.
The error because the derivative of Bessel function with respect to the integer rank is not defined, and the code return error "ChainRulesCore.@thunk(error("not implemented"))"
Since adjoint calculation will never require calculation gradient with respect to the integer rank, I am defining this integer as zero or `ChainRulesCore.Zero()`
"""

#
ChainRulesCore.@scalar_rule(
    besselj(ν, x),
    (
        ChainRulesCore.Zero(), # ChainRulesCore.@thunk(error("not implemented")),
        (besselj(ν - 1, x) - besselj(ν + 1, x)) / 2
    ),
)
ChainRulesCore.@scalar_rule(
    besseli(ν, x),
    (
        ChainRulesCore.Zero(), # ChainRulesCore.@thunk(error("not implemented")),
        (besseli(ν - 1, x) + besseli(ν + 1, x)) / 2,
    ),
)
ChainRulesCore.@scalar_rule(
    bessely(ν, x),
    (
        ChainRulesCore.Zero(), # ChainRulesCore.@thunk(error("not implemented")),
        (bessely(ν - 1, x) - bessely(ν + 1, x)) / 2,
    ),
)

# TODO: I need to define scalar_rule for `wignerdjmn`, like following. The problem is that `wignerdjmn` is not recogized by `ChainRulesCore`
"""
ChainRulesCore.@scalar_rule(
    wignerdjmn(s, m, n, θ),
    (
        ChainRulesCore.Zero(), # ChainRulesCore.@thunk(error("not implemented")),
        ChainRulesCore.Zero(), # ChainRulesCore.@thunk(error("not implemented")),
        ChainRulesCore.Zero(), # ChainRulesCore.@thunk(error("not implemented")),
        ∂wignerdjmn_by_∂θ(s, m, n, θ),
    ),
)
"""

# Wigner-d
"""
    Wigner-d function calculated from eq. B.1
This is calculated directly, not with recurrence relation.
The Wigner-d calculated by recurrence is more numerically stable.
Use this function only as a validation for the recurrence relation, or when automatic differentiation is needed.

TODO: for large s,m,n, this is not stable. Use the one with recurrence instead
"""
function wignerdjmn_ELZOUKA(s::Int, m::Int, n::Int, θ::R) where R <: Real
    # println("s=$s, m=$m, n=$n, θ=$θ")
    if θ == zero(θ) # TODO: make the zero the same type of θ. e.g., zero(θ)
        d = δ(m, n)
    elseif θ == π
        d = (-1)^(s - n) * δ(-n, m)
    else
        d = zero(θ)
        k_min = max(0, m - n)
        k_max = min(s + m, s - n)
        if k_max >= k_min
            for k in k_min:k_max
                d += (-1)^k *
                        (cos(θ / 2)^(2s - 2k + m - n) * sin(θ / 2)^(2k - m + n)) /
                        (factorial(k) * factorial(s + m - k) * factorial(s - n - k) * factorial(n - m + k))
            end
            d *= sqrt(factorial(s + m) * factorial(s - m) * factorial(s + n) * factorial(s - n))
        else # wigner-d is zero if there is any negative factorial
            return zero(θ)
        end
    end

    return d
end
wignerdjmn = wignerdjmn_ELZOUKA # I did it to make it work with auto-diff, although "wignerdjmn_ELZOUKA" is not efficient.
# I may need to define "ChainRulesCore.@scalar_rule" for "WignerD.wignerdjmn"
# wignerdjmn = WignerD.wignerdjmn

# derivative of wigner-D
function ∂wignerdjmn_by_∂θ(s::Int, m::Int, n::Int, θ::R; numerical_derivative=false, verysmallnumber=1e-30) where R <: Real
    """
    derivative of wigner-d with resepect to θ. Adopted from eq. B.25 from Mishchenko, M.I., Travis, L.D., and Lacis, A.A. (2002). Scattering, absorption, and emission of light by small particles (Cambridge University Press).
    """
    if numerical_derivative
        return (wignerdjmn(s, m, n, θ + verysmallnumber) - wignerdjmn(s, m, n, θ - verysmallnumber)) / (verysmallnumber * 2)
    else
        # try
        # return    (m - n * cos(θ)) / sin(θ) * wignerdjmn(s, m, n, θ) + sqrt((s + n) * (s - n + 1)) * wignerdjmn(s, m, n - 1, θ)
        # catch
        return -1 * (m - n * cos(θ)) / sin(θ) * wignerdjmn(s, m, n, θ) - sqrt((s - n) * (s + n + 1)) * wignerdjmn(s, m, n + 1, θ)
        # end
    end
end



#############################################################################################
# calculate π(θ) and τ(θ)
function πₘₙ(m::Int, n::Int, θ::R) where R <: Real
    return m / sin(θ) * wignerdjmn(n, 0, m, θ)
end

function τₘₙ(m::Int, n::Int, θ::R) where R <: Real
    return ∂wignerdjmn_by_∂θ(n, 0, m, θ)
end


function γ_mn(m::Int, n::Int)
    return sqrt(
        ((2n + 1) * factorial(n - m)) / (4π * n * (n + 1) * factorial(n + m))
    ) # equation C.22
end

function γ_mn_dash(m::Int, n::Int)
    return sqrt(
        ((2n + 1) * factorial(n - m)) / (4π * factorial(n + m))
    ) # equation C.25
end

function get_radial_function_and_special_derivative_given_kind(kind)
    """
    return the radial function that is appropriate with the type of VSWF
    """

    if kind in ["regular", "incoming", 1]
        radial_function = spherical_Bessel_j_n
        radial_function_special_derivative = one_over_x_by_∂_x_j_n_by_∂x
    elseif kind in ["irregular", "outgoing", 2]
        radial_function = spherical_Hankel_h1_n
        radial_function_special_derivative = one_over_x_by_∂_x_h_n_by_∂x
    else
        throw(DomainError(""" 'kind' has to be one of the following: ["regular", "incoming", "irregular", "outgoing"] """))
    end

    return radial_function, radial_function_special_derivative
end
