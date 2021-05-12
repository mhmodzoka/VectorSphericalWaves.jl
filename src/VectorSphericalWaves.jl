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




#############################################################################################
# import "WignerD.jl" code to calculate wigner-d (https://github.com/jishnub/WignerD.jl)
# [TODO] I have already installed "WignerD" package using the following commands: # include("WignerD.jl/src/WignerD.jl") ;  import Pkg; Pkg.add("WignerD")
# include("WignerD.jl/src/WignerD.jl")
# using .WignerD


using SpecialFunctions
using DataStructures
using StaticArrays

import WignerD

# exporting main functions
export M_mn_wave
export N_mn_wave

export B_mn_of_θ
export C_mn_of_θ
export P_mn_of_θ

export B_mn_of_θ_ϕ
export C_mn_of_θ_ϕ
export P_mn_of_θ_ϕ

export γ_mn
export γ_mn_dash

# exporting functions that calculate for all m,n combinations
export M_N_wave_all_m_n
export B_C_P_mn_of_θ_ϕ_for_all_m_n



export wignerdjmn_ELZOUKA

# export wignerdjmn
export ∂wignerdjmn_by_∂θ
export single_index_from_m_n
export get_max_single_index_from_n_max
export πₘₙ_τₘₙ_all_all_m_n

#############################################################################################
δ(x,y) = ==(x, y)

function single_index_from_m_n(m, n)
    return n * (n + 1) + m
end

function get_max_single_index_from_n_max(n_max)
    return single_index_from_m_n(n_max, n_max)
end


#############################################################
"""
    defining derivative for Bessel functions.
Inside `SpecialFunctions.jl` package, the definition of derivative of Bessel functions exists. However, it returns an error when called by Zygote.
The error because the derivative of Bessel function with respect to the integer rank is not defined, and the code return error `ChainRulesCore.@thunk(error("not implemented"))`
Since adjoint calculation will never require calculation gradient with respect to the integer rank, I am defining this integer as zero or `ChainRulesCore.Zero()`
"""
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

# Wigner-d
"""
Wigner-d function calculated from eq. B.1
"""
function wignerdjmn_ELZOUKA(s, m, n, θ)
    # println("s=$s, m=$m, n=$n, θ=$θ")
    if θ == 0
        d = δ(m, n)
    elseif θ == π
        d = (-1)^(s - n) * δ(-n, m)
    else
        d = 0
        for k in max(0, m - n):min(s + m, s - n)
            d += (-1)^k * 
                    (cos(θ / 2)^(2s - 2k + m - n) * sin(θ / 2)^(2k - m + n)) / 
                    (factorial(k) * factorial(s + m - k) * factorial(s - n - k) * factorial(n - m + k))
        end
        d *= sqrt(factorial(s + m) * factorial(s - m) * factorial(s + n) * factorial(s - n))
    end
    
    return d
end
wignerdjmn = wignerdjmn_ELZOUKA # I did it to make it work with auto-diff
# wignerdjmn = WignerD.wignerdjmn

# derivative of wigner-D
function ∂wignerdjmn_by_∂θ(s, m, n, θ; numerical_derivative=false, verysmallnumber=1e-30)
    """
    derivative of wigner-d with resepect to θ. Adopted from eq. B.25 from Mishchenko, M.I., Travis, L.D., and Lacis, A.A. (2002). Scattering, absorption, and emission of light by small particles (Cambridge University Press).
    """    
    if numerical_derivative
        return (wignerdjmn(s, m, n, θ + verysmallnumber) - wignerdjmn(s, m, n, θ - verysmallnumber)) / (verysmallnumber * 2)
    else
        # try
        return    (m - n * cos(θ)) / sin(θ) * wignerdjmn(s, m, n, θ) + sqrt((s + n) * (s - n + 1)) * wignerdjmn(s, m, n - 1, θ)
        # catch
        #    return -1 * (m - n * cos(θ)) / sin(θ) * wignerdjmn(s, m, n, θ) - sqrt((s - n) * (s + n + 1)) * wignerdjmn(s, m, n + 1, θ)
        # end
    end
end



#############################################################################################
# calculate π(θ) and τ(θ)
function πₘₙ(m, n, θ)
    return m / sin(θ) * wignerdjmn(n, 0, m, θ)
end

function τₘₙ(m, n, θ)
    return ∂wignerdjmn_by_∂θ(n, 0, m, θ)
end


#############################################################################################
# calculate B(θ), C(θ), P(θ)
function B_mn_of_θ(m, n, θ)
    """
    I assume each of m, n, θ is a single number
    """
    # TODO you can use literal syntax for this
    return vcat(
        0,                  # r-component
        τₘₙ(m, n, θ),      # θ-component
        im * πₘₙ(m, n, θ)  # ϕ-component      
    ) # equation C.19
end

function C_mn_of_θ(m, n, θ)
    """
    I assume each of m, n, θ is a single number
    """
    return vcat(
        0,                  # r-component
        im * πₘₙ(m, n, θ), # θ-component
        -1 * τₘₙ(m, n, θ),    # ϕ-component      
    ) # equation C.20
end

# TODO use partial application for M,N, fold them into the type sig of struct P
function P_mn_of_θ(m, n, θ)
    return vcat(
        wignerdjmn(n, 0, m, θ), # r-component
        0,                      # θ-component
        0,                      # ϕ-component      
    ) # equation C.21
end


#############################################################################################
# calculate B(θ,ϕ), C(θ,ϕ), P(θ,ϕ)
function B_mn_of_θ_ϕ(m, n, θ, ϕ)
    return (-1)^m * sqrt(factorial(n + m) / factorial(n - m)) * B_mn_of_θ(m, n, θ) * exp(im * m * ϕ) # equation C.16
end

function C_mn_of_θ_ϕ(m, n, θ, ϕ)
    return (-1)^m * sqrt(factorial(n + m) / factorial(n - m)) * C_mn_of_θ(m, n, θ) * exp(im * m * ϕ) # equation C.17
end

function P_mn_of_θ_ϕ(m, n, θ, ϕ)
    return (-1)^m * sqrt(factorial(n + m) / factorial(n - m)) * P_mn_of_θ(m, n, θ) * exp(im * m * ϕ) # equation C.17
end


#############################################################################################
# calculate spherical Bessel and Hankel functions and their derivarive
function spherical_Bessel_j_n(n, x)
    """
    Spherical Bessel function of the first kind.
    It can be calculated from ordinary Bessel function of the first kind "besselj" as in the code.
    check https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions
    """
    return √(π / 2x) * besselj(n + 1 / 2, x)
end

function spherical_Bessel_y_n(n, x)
    """
    Spherical Bessel function of the second kind.
    It can be calculated from ordinary Bessel function of the second kind "bessely" as in the code.
    check https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions
    """
    return √(π / 2x) * bessely(n + 1 / 2, x)
end

function spherical_Hankel_h1_n(n, x)
    """
    Spherical Hankel function of the first kind. It can be calculated from spherical Bessel functions of the first and second kinds as in the code:
    """
    return spherical_Bessel_j_n(n, x) + im * spherical_Bessel_y_n(n, x)
end

function one_over_x_by_∂_x_j_n_by_∂x(n, x)
    """
    Derivative of (spherical Bessel of first kind * x) divided by x
    """
    return (spherical_Bessel_j_n(n - 1, x) - n / x * spherical_Bessel_j_n(n, x))
end

function one_over_x_by_∂_x_y_n_by_∂x(n, x)
    """
    Derivative of (spherical Bessel of second kind * x) divided by x
    """
    return (spherical_Bessel_y_n(n - 1, x) - n / x * spherical_Bessel_y_n(n, x))
end

function one_over_x_by_∂_x_h_n_by_∂x(n, x)
    """
    Derivative of (spherical Hankel of first kind * x) divided by x
    """
    return one_over_x_by_∂_x_j_n_by_∂x(n, x) + im * one_over_x_by_∂_x_y_n_by_∂x(n, x)
end


#############################################################################################
# calculate (Rg)M(kr,θ,ϕ), (Rg)N(kr,θ,ϕ)
function γ_mn(m, n)
    return sqrt(
        ((2n + 1) * factorial(n - m)) / (4π * n * (n + 1) * factorial(n + m))
    ) # equation C.22
end

function γ_mn_dash(m, n)
    return sqrt(
        ((2n + 1) * factorial(n - m)) / (4π * factorial(n + m))
    ) # equation C.25
end

function get_radial_function_and_special_derivative_given_kind(kind)
    """
    return the radial function that is appropriate with the type of VSWF
    """

    if kind in ["regular", "incoming"]
        radial_function = spherical_Bessel_j_n
        radial_function_special_derivative = one_over_x_by_∂_x_j_n_by_∂x
    elseif kind in ["irregular", "outgoing"]
        radial_function = spherical_Hankel_h1_n
        radial_function_special_derivative = one_over_x_by_∂_x_h_n_by_∂x
    else
        throw(DomainError(""" 'kind' has to be one of the following: ["regular", "incoming", "irregular", "outgoing"] """))
    end

    return radial_function, radial_function_special_derivative

end

function M_mn_wave(m, n, kr, θ, ϕ; kind="regular")
    """
    Parameters
    ==========
    kind: string, either ["regular" or "incoming"] or ["irregular" or "outgoing"]
    """

    radial_function, _ = get_radial_function_and_special_derivative_given_kind(kind)
    return γ_mn(m, n) * radial_function(n, kr) * C_mn_of_θ_ϕ(m, n, θ, ϕ)    
end

function N_mn_wave(m, n, kr, θ, ϕ; kind="regular")
    """
    Parameters
    ==========
    kind: string, either ["regular" or "incoming"] or ["irregular" or "outgoing"]
    """

    radial_function, radial_function_special_derivative  = get_radial_function_and_special_derivative_given_kind(kind)
    return γ_mn(m, n) * (
        n * (n + 1) / kr * radial_function(n, kr) * P_mn_of_θ_ϕ(m, n, θ, ϕ)
        + (radial_function_special_derivative(n, kr) * B_mn_of_θ_ϕ(m, n, θ, ϕ))
    )
end


#############################################################################################
# This version is different from "VectorSphericalHarmonics.jl", as it calculate VSWF for all m,n and all kr, θ, ϕ in one call
# This should be more stable and faster, as it employs recurrence relations in calculating πₘₙ and τₘₙ
# using recurrence relations in calculating πₘₙ and τₘₙ

#############################################################################################
# calculate π(θ) and τ(θ) using recurrence relations
function πₘₙ_τₘₙ_all_all_m_n(n_max, θ; verbose=false)
    # calculate A with recurrence relation
    A = SortedDict(0 => 1.0)
    for m = 0:n_max - 1
        A[m + 1] = A[m] * sqrt((2m + 1) / (2 * (m + 1)))
    end

    # calculate πₘₙ with recurrence relation
    πₘₙ_all = zeros(get_max_single_index_from_n_max(n_max))
    for m = 1:n_max   # TODO: I think I need to start m from 0, not 1.     
        n = m        
        πₘₙ_all[single_index_from_m_n(m, n)] = m * A[m] * sin(θ)^(m - 1)
        if verbose; println("m=$m, n=$n, πₘₙ_all=$(πₘₙ_all[single_index_from_m_n(m, n)])"); end

        # calculate π₋ₘₙ from πₘₙ
        if m != 0
            πₘₙ_all[single_index_from_m_n(-m, n)] = (-1)^(m + 1) * πₘₙ_all[single_index_from_m_n(m, n)]
            if verbose; println("m=$(-m), n=$n, πₘₙ_all=$(πₘₙ_all[single_index_from_m_n(-m, n)])"); end
        end

        for n = (m + 1):n_max
            if n == m + 1
                πₘₙ_all[single_index_from_m_n(m, n)] = 1 / sqrt(n^2 - m^2) * ((2n - 1) * cos(θ) * πₘₙ_all[single_index_from_m_n(m, n - 1)])
                if verbose; println("m=$m, n=$n, πₘₙ_all=$(πₘₙ_all[single_index_from_m_n(m, n)])"); end
            else
                πₘₙ_all[single_index_from_m_n(m, n)] = 1 / sqrt(n^2 - m^2) * ((2n - 1) * cos(θ) * πₘₙ_all[single_index_from_m_n(m, n - 1)]) - sqrt((n - 1)^2 - m^2) * πₘₙ_all[single_index_from_m_n(m, n - 2)]
                if verbose; println("m=$m, n=$n, πₘₙ_all=$(πₘₙ_all[single_index_from_m_n(m, n)])"); end
            end

            # calculate π₋ₘₙ from πₘₙ
            if m != 0
                πₘₙ_all[single_index_from_m_n(-m, n)] = (-1)^(m + 1) * πₘₙ_all[single_index_from_m_n(m, n)]
                if verbose; println("m=$(-m), n=$n, πₘₙ_all=$(πₘₙ_all[single_index_from_m_n(-m, n)])"); end
            end
        end
    end
    τₘₙ_all = 0 # TODO: add the code for it
    return πₘₙ_all, τₘₙ_all 
end





#############################################################################################
# Legendre and Associated Legendre
function Legendre_polynomials_Pn_array(n_max, x)
    """
    Calculate all Legendre polynomials Pₙ(x) from n = 1 up to n=n_max using recurrence relation
    https://en.wikipedia.org/wiki/Legendre_polynomials

    returns
    dictionary. TODO: find a better way if this cause adjoint calculation trouble
    """
    P = SortedDict(
        0 => 1,
        1 => x,
    )

    for n = 1:(n_max - 1)
        P[n + 1] = ((2n + 1) * x * P[n] - n * P[n - 1]) / (n + 1)
    end

    return P
end

function Legendre_polynomials_Pn(n, x)
    """
    Calculate Legendre polynomials Pₙ(x) at a given n
    """    
    return Legendre_polynomials_Pn_array(n, x)[n]
end


function Associated_Legendre_polynomials_Pmn_array(m, n_max, x)
    """
    Associated Legendre polynomials Pᵐₙ(x) from n = m up to n=n_max using recurrence relation
    https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Recurrence_formula

    returns
    dictionary. TODO: find a better way if this cause adjoint calculation trouble
    """
    n = m
    P = SortedDict(
        n => (-1)^n * factorial(factorial(2n - 1)) * (1 - x^2)^(n / 2),
    )
    P[n + 1] = x * (2n + 1) * P[n]

    for n = (m + 1):(n_max - 1)
        P[n + 1] = ( (2n + 1) * x * P[n] - (n + m) * P[n - 1] ) / (n - m + 1)
    end

    return P
end

function Associated_Legendre_polynomials_Pmn_array(m, n, x)
    """
    Associated Legendre polynomials Pᵐₙ(x) at a given n
    """
    return Associated_Legendre_polynomials_Pmn_array(m, n, x)[n]
end

#############################################################################################
# Wigner-d, using recurrence
function wignerd_and_∂wignerd_for_all_s(s_max, m, n, θ; get_derivatives=true, verbose=false)
    """
    Calculate dˢₘₙ(θ) and ∂(dˢₘₙ(θ))/∂θ for all values of s, where s starts from s_min up to s_max. s_min is the maximum of |m| and |n|

    """
    # TODO: special case of m=0, n=0 (eq. B.27)
    s_min = max(abs(m), abs(n))

    # calculate ξ from eq. B.16
    if n >= m
        ξ_mn = 1
    else
        ξ_mn = (-1)^(m - n)
    end
    
    x = cos(θ)

    # calculate d^s_min__m_n(θ) from eq. B.124
    d_smin_m_n = 
        ξ_mn * 2.0^(-s_min) * sqrt(
            factorial(BigInt(2s_min)) / (factorial(BigInt(abs(m - n))) * factorial(BigInt(abs(m + n))))
        ) *
        (1 - x)^(abs(m - n) / 2) *
        (1 + x)^(abs(m + n) / 2)

    d = SortedDict(
        s_min - 1 => 0.0,
        s_min   => d_smin_m_n
    )    

    # applying the recurrence relation
    if n == 0
        if m == 0
            if verbose; println("m=$m, n=$n, I will use Legendre_polynomials_Pn_array"); end
            P_s = Legendre_polynomials_Pn_array(s_max + 1, x)
            for s = s_min:s_max + 1
                d[s] = P_s[s]
            end
        else
            if verbose; println("m=$m, n=$n, I will use Associated_Legendre_polynomials_Pmn_array"); end
            P_m_s = Associated_Legendre_polynomials_Pmn_array(m, s_max + 1, x)
            for s = s_min:s_max + 1
                d[s] = sqrt(factorial(s - m) / factorial(s + m)) * P_m_s[s]
            end            
        end
    else
        for s = s_min:s_max + 1 
            if verbose; println("m=$m, n=$n, I will use the general recurrence"); end
            d[s + 1] = 1 / (s * sqrt((s + 1)^2 - m^2) * sqrt((s + 1)^2 - n^2)) * (
                (2s + 1) * (s * (s + 1) * x - m * n) * d[s]
                - 1 * (s + 1) * sqrt(s^2 - m^2) * sqrt(s^2 - n^2) * d[s - 1]
            ) # eq. B.22
        end
    end

    # calculate the derivative ∂(dˢₘₙ(θ))/∂θ
    if get_derivatives        
        ∂d_∂θ = SortedDict(
            s_min - 1 => 0.0,
            s_min   => 0.0
        )

        sin_theta = sin(θ)

        for s = s_min:s_max
            if verbose; println("s=$s"); end
            ∂d_∂θ[s] = 1 / sin_theta * (
                -1 * ((s + 1) * sqrt((s^2 - m^2) * (s^2 - n^2))) / (s * (2s + 1)) * d[s - 1]
                - 1 * (m * n) / (s * (s + 1)) * d[s]
                + 1 * (
                        s *
                        sqrt((s + 1)^2 - m^2) *
                        sqrt((s + 1)^2 - n^2)
                    ) / 
                    ((s + 1) * (2s + 1)) * d[s + 1]
            )
        end
        return d, ∂d_∂θ
    else
        return d
    end

    
end




#############################################################################################
# calculate π(θ) and τ(θ) using Wigner-d that was calculated using recurrence relations
function πₘₙ_τₘₙ_all_all_m_n_using_wigner(n_max, θ; verbose=false)
    πₘₙ_all = zeros(get_max_single_index_from_n_max(n_max))
    τₘₙ_all = zeros(get_max_single_index_from_n_max(n_max))

    for m = 0:n_max  # TODO: special case of m=0      
        d_rec_all_n, ∂d_∂θ_rec_all_n  = wignerd_and_∂wignerd_for_all_s(n_max, 0, m, θ)
        for n = m:n_max
            if n != 0
                if verbose; println("m=$m, n=$n, calculate πₘₙ_all, τₘₙ_all"); end
                πₘₙ_all[single_index_from_m_n(m, n)] = m / sin(θ) * d_rec_all_n[n]
                τₘₙ_all[single_index_from_m_n(m, n)] = ∂d_∂θ_rec_all_n[n]
                
                if m != 0
                    πₘₙ_all[single_index_from_m_n(-m, n)] = (-1)^(m + 1) * πₘₙ_all[single_index_from_m_n(m, n)]
                    τₘₙ_all[single_index_from_m_n(-m, n)] = (-1)^(m)   * τₘₙ_all[single_index_from_m_n(m, n)]
                end
            end
        end
    end
    return πₘₙ_all, τₘₙ_all
end

function B_C_mn_of_θ_for_all_m_n(n_max, θ)
    """
    Calculate Bₙₘ(θ), Cₙₘ(θ), Pₙₘ(θ) equations C.19, C.20, C.21
    The order of m,n is according to the function "single_index_from_m_n"
    """
    πₘₙ_all, τₘₙ_all = πₘₙ_τₘₙ_all_all_m_n_using_wigner(n_max, θ)
    B = (_ -> zero(SVector{3,Complex})).(πₘₙ_all)
    C = (_ -> zero(SVector{3,Complex})).(πₘₙ_all)    
    for idx in eachindex(πₘₙ_all)
        B[idx] = [0,      τₘₙ_all[idx], im * πₘₙ_all[idx] ]
        C[idx] = [0, im * πₘₙ_all[idx], -1 * τₘₙ_all[idx] ]    
    end

    return B, C
end

function P_mn_of_θ_for_all_m_n(n_max, θ)
    """
    Calculate Pₙₘ(θ) using equations C.19, C.20
    The order of m,n is according to the function "single_index_from_m_n"
    """
    P = fill(zero(SVector{3,Complex}), get_max_single_index_from_n_max(n_max))
    for m = 0:n_max
        d_rec_all_n  = wignerd_and_∂wignerd_for_all_s(n_max, 0, m, θ; get_derivatives=false)
        for n = m:n_max
            if n != 0                
                P[single_index_from_m_n(m, n)] = [d_rec_all_n[n], 0, 0]                
                P[single_index_from_m_n(-m, n)] = [(-1)^m * d_rec_all_n[n], 0, 0]  # using symmetry relation B.5, we can get d_(0,-m) from d_(0,m)
            end
        end
    end
    return P
end

function B_C_P_mn_of_θ_ϕ_for_all_m_n(n_max, θ, ϕ)
    """
    Calculate Bₙₘ(θ,ϕ), Cₙₘ(θ,ϕ), Pₙₘ(θ,ϕ) for all m and n
    """
    B_of_θ, C_of_θ = B_C_mn_of_θ_for_all_m_n(n_max, θ)
    P_of_θ = P_mn_of_θ_for_all_m_n(n_max, θ)

    for n = 1:n_max
        for m = -n:n
            factor = (-1)^m * sqrt(factorial(n + m) / factorial(n - m)) * exp(im * m * ϕ)
            B_of_θ[single_index_from_m_n(m, n)] *= factor
            C_of_θ[single_index_from_m_n(m, n)] *= factor
            P_of_θ[single_index_from_m_n(m, n)] *= factor
        end
    end

    return B_of_θ, C_of_θ, P_of_θ    
end

function M_N_wave_all_m_n(n_max, kr, θ, ϕ; kind="regular")
    """
    Parameters
    ==========
    kind: string, either ["regular" or "incoming"] or ["irregular" or "outgoing"]
    """
    radial_function, radial_function_special_derivative  = get_radial_function_and_special_derivative_given_kind(kind)

    B_of_θ_ϕ, C_of_θ_ϕ, P_of_θ_ϕ  = B_C_P_mn_of_θ_ϕ_for_all_m_n(n_max, θ, ϕ)

    M = 0 .* B_of_θ_ϕ
    N = 0 .* B_of_θ_ϕ

    for n = 1:n_max
        for m = -n:n
            M[single_index_from_m_n(m, n)] = γ_mn(m, n) * radial_function(n, kr) * C_of_θ_ϕ[single_index_from_m_n(m, n)]    
            N[single_index_from_m_n(m, n)] = γ_mn(m, n) * (
                n * (n + 1) / kr * radial_function(n, kr)    * P_of_θ_ϕ[single_index_from_m_n(m, n)]
                + (radial_function_special_derivative(n, kr) * B_of_θ_ϕ[single_index_from_m_n(m, n)])
            )
        end
    end
    return M, N
end




end # module
