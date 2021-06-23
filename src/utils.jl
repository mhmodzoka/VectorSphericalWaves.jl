using SpecialFunctions
using DataStructures
using StaticArrays
using ComplexOperations
using Memoize
using LoopVectorization

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

TODO: for large s,m,n, this function is not stable due to large integers resulted from factorial calclulation and their multiplication. Here, I am using big() to avoid the overflow error. I need to implement when I should apply big(), while not affecting performance and autodifferentiation

"""
function wignerdjmn_ELZOUKA(s::I, m::I, n::I, θ::R) where {R <: Real,I <: Integer}
    # println("s=$s, m=$m, n=$n, θ=$θ")
    if θ == zero(θ) # TODO: make the zero the same type of θ. e.g., zero(θ)
        d = δ(m, n);
    elseif θ == π
        d = (-1)^(s - n) * δ(-n, m);
    else
        d = zero(θ)
        
        
        k_max = min(s + m, s - n)
        # TODO: find a better way to detect when do we need to use big integers
        if (
                (max(k_max, s + m - k_max, s - n - k_max, n - m + k_max, s + abs(m), s + abs(n)) >= 21) ||
                ((factorial(s + m) * factorial(s - m) * factorial(s + n) * factorial(s - n)) < 0) ||
                ((factorial(k_max) * factorial(s + m - k_max) * factorial(s - n - k_max) * factorial(n - m + k_max)) < 0)
            )
            s = big(s); m = big(m); n = big(n)            
            k_max = min(s + m, s - n)
        end
        k_min = max(0, m - n)

        if k_max >= k_min
            for k in k_min:k_max
                d += (-1)^k *
                        (cos(θ / 2)^(2s - 2k + m - n) * sin(θ / 2)^(2k - m + n)) /
                        (factorial(k) * factorial(s + m - k) * factorial(s - n - k) * factorial(n - m + k))
            end
            d *= sqrt(factorial(s + m) * factorial(s - m) * factorial(s + n) * factorial(s - n))
            return convert(typeof(θ), d)
        else # wigner-d is zero if there is any negative factorial
            return zero(θ)
        end
    end

    return convert(typeof(θ), d)
end



"""
    Calculate Wigner-d function, using recurrence and memoize
"""
@memoize function wignerdjmn_recurrence_memoize(s::I, m::I, n::I, θ::R) where {R <: Real,I <: Integer}

    s_min = max(abs(m), abs(n))

    if typeof(s) == Int && 2s_min >= 21
        s = big(s); m = big(m); n = big(n)
        s, s_min = promote(s, s_min)                        
    end

    x = cos(θ)

    if m > n # I need to memoize only m < n, because it is relevant to calculating angular functions
        return (-1)^(m - n) * wignerdjmn_recurrence_memoize(s, n, m, θ) # eq. B.5

    elseif θ < 0 # I need to memoize only θ > 0
        return (-1)^(m - n) * wignerdjmn_recurrence_memoize(s, m, n, -θ)

    elseif n < 0 # I need to memoize only n > 0
        return (-1)^(m - n) * wignerdjmn_recurrence_memoize(s, -m, -n, θ)

    elseif s < s_min
        return zero(θ)

    elseif s == s_min
        if typeof(s) == Int && (factorial((abs(m - n))) * factorial((abs(m + n)))) < 0 # check if the denominator may cause overflow
            s = big(s); m = big(m); n = big(n)
            s, s_min = promote(s, s_min)
        end
        # calculate ξ from eq. B.16
        if n >= m
            ξ_mn = one(s)
        else
            ξ_mn = (-1)^(m - n)
        end

        # calculate d^s_min__m_n(θ) from eq. B.24
        return convert(typeof(θ),
            ξ_mn * 2.0^(-s_min) * sqrt(
                factorial((2s_min)) / (factorial((abs(m - n))) * factorial((abs(m + n))))
            ) *
            (1 - x)^(abs(m - n) / 2) *
            (1 + x)^(abs(m + n) / 2)
        )

    else
        d_s_here_plus_1 = zero(θ)
        for s_here = s_min:s - 1
            d_s_here_plus_1 = 1 / (s_here * sqrt((s_here + 1)^2 - m^2) * sqrt((s_here + 1)^2 - n^2)) * (
                (2s_here + 1) * (s_here * (s_here + 1) * x - m * n) * wignerdjmn_recurrence_memoize(s_here, m, n, θ)
                - 1 * (s_here + 1) * sqrt(s_here^2 - m^2) * sqrt(s_here^2 - n^2) * wignerdjmn_recurrence_memoize(s_here - 1, m, n, θ)
            ) # eq. B.22
        end
        return convert(typeof(θ), d_s_here_plus_1)
    end
end


"""
    Calculate Wigner-d function, using recurrence
"""
function wignerdjmn_recurrence(s::I, m::I, n::I, θ::R) where {R <: Real,I <: Integer}

    s_min = max(abs(m), abs(n))

    if typeof(s) == Int && 2s_min >= 21
        s = big(s); m = big(m); n = big(n)
        s, s_min = promote(s, s_min)            
    end

    x = cos(θ)

    if m > n # I need to memoize only m < n, because it is relevant to calculating angular functions
        return (-1)^(m - n) * wignerdjmn_recurrence(s, n, m, θ) # eq. B.5

    elseif θ < 0 # I need to memoize only θ > 0
        return (-1)^(m - n) * wignerdjmn_recurrence(s, m, n, -θ)

    elseif n < 0 # I need to memoize only n > 0
        return (-1)^(m - n) * wignerdjmn_recurrence(s, -m, -n, θ)

    elseif s < s_min
        return zero(θ)

    elseif s == s_min
        if typeof(s) == Int && (factorial((abs(m - n))) * factorial((abs(m + n)))) < 0 # check if the denominator may cause overflow
            s = big(s); m = big(m); n = big(n)
            s, s_min = promote(s, s_min)
        end
        # calculate ξ from eq. B.16
        if n >= m
            ξ_mn = one(s)
        else
            ξ_mn = (-1)^(m - n)
        end

        # calculate d^s_min__m_n(θ) from eq. B.24
        return convert(typeof(θ),
            ξ_mn * 2.0^(-s_min) * sqrt(
                factorial((2s_min)) / (factorial((abs(m - n))) * factorial((abs(m + n))))
            ) *
            (1 - x)^(abs(m - n) / 2) *
            (1 + x)^(abs(m + n) / 2)
        )

    else
        d_s_here_plus_1 = zero(θ)
        s_here = s_min
        d_s_here_minus_1 = wignerdjmn_recurrence(s_here - 1, m, n, θ)
        d_s_here = wignerdjmn_recurrence(s_here, m, n, θ)
        for s_here = s_min:s - 1
            d_s_here_plus_1 = 1 / (s_here * sqrt((s_here + 1)^2 - m^2) * sqrt((s_here + 1)^2 - n^2)) * (
                (2s_here + 1) * (s_here * (s_here + 1) * x - m * n) * d_s_here
                - 1 * (s_here + 1) * sqrt(s_here^2 - m^2) * sqrt(s_here^2 - n^2) * d_s_here_minus_1
            ) # eq. B.22
            d_s_here_minus_1 = d_s_here
            d_s_here = d_s_here_plus_1            
        end
        return convert(typeof(θ), d_s_here_plus_1)
    end
end




wignerdjmn = wignerdjmn_ELZOUKA # I did it to make it work with auto-diff, although "wignerdjmn_ELZOUKA" is not efficient.
# wignerdjmn = wignerdjmn_recurrence_memoize
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
        try
            return    (m - n * cos(θ)) / sin(θ) * wignerdjmn(s, m, n, θ) + sqrt((s + n) * (s - n + 1)) * wignerdjmn(s, m, n - 1, θ)
        catch
            return -1 * (m - n * cos(θ)) / sin(θ) * wignerdjmn(s, m, n, θ) - sqrt((s - n) * (s + n + 1)) * wignerdjmn(s, m, n + 1, θ)
        end
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
