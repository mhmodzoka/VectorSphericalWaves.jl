# TODO: I cound that "utils.jl" has been `include`ed in the main `VectorSphericalWaves` module. Is there a more neat way? when I `include("utils.jl")` here, I get lots of warnings.

export M_mn_wave
export N_mn_wave
export M_mn_wave_SVector
export N_mn_wave_SVector


export B_mn_of_θ
export C_mn_of_θ
export P_mn_of_θ

export B_mn_of_θ_ϕ
export C_mn_of_θ_ϕ
export P_mn_of_θ_ϕ

#############################################################################################
# calculate B(θ), C(θ), P(θ)
"""
    Calculation of B(θ), equation C.19, returns array
I assume each of m, n, θ is a single number
"""
function B_mn_of_θ(m::Int, n::Int, θ::R) where R <: Real
    return vcat(
        zero(θ),                  # r-component
        τₘₙ(m, n, θ),      # θ-component
        im * πₘₙ(m, n, θ)  # ϕ-component
    ) # equation C.19
end

"""
    Calculation of B(θ), equation C.19, returns SVector
I assume each of m, n, θ is a single number
"""
function B_mn_of_θ_SVector(m::Int, n::Int, θ::R) where R <: Real
    """
    I assume each of m, n, θ is a single number
    """
    return SVector(
        zero(θ),                  # r-component
        τₘₙ(m, n, θ),      # θ-component
        im * πₘₙ(m, n, θ)  # ϕ-component
    ) # equation C.19
end

"""
    Calculation of C(θ), equation C.20, returns array
I assume each of m, n, θ is a single number
"""
function C_mn_of_θ(m::Int, n::Int, θ::R) where R <: Real    
    return vcat(
        zero(θ),                  # r-component
        im * πₘₙ(m, n, θ), # θ-component
        -1 * τₘₙ(m, n, θ),    # ϕ-component
    ) # equation C.20
end

"""
    Calculation of C(θ), equation C.20, returns SVector
I assume each of m, n, θ is a single number
"""
function C_mn_of_θ_SVector(m::Int, n::Int, θ::R) where R <: Real   
    return SVector(
        zero(θ),                  # r-component
        im * πₘₙ(m, n, θ), # θ-component
        -1 * τₘₙ(m, n, θ),    # ϕ-component
    ) # equation C.20
end

# TODO use partial application for M,N, fold them into the type sig of struct P
"""
    equation C.21, returns array
"""
function P_mn_of_θ(m::Int, n::Int, θ::R) where R <: Real
    # TODO @Alok, should I replace arrays with SMatrix? what are the drawbacks?
    # TODO: replace "0" with zero(type)
    return vcat(
        wignerdjmn(n, 0, m, θ), # r-component
        zero(θ),                # θ-component
        zero(θ),                # ϕ-component
    ) # equation C.21
end

"""
    equation C.21, returns SVector
"""
function P_mn_of_θ_SVector(m::Int, n::Int, θ::R) where R <: Real
    # TODO @Alok, should I replace arrays with SMatrix? what are the drawbacks?
    # TODO: replace "0" with zero(type)
    return SVector(
        wignerdjmn(n, 0, m, θ), # r-component
        zero(θ),                # θ-component
        zero(θ),                # ϕ-component
    ) # equation C.21
end


#############################################################################################
# calculate B(θ,ϕ), C(θ,ϕ), P(θ,ϕ), returns Array
function B_mn_of_θ_ϕ(m::Int, n::Int, θ::R, ϕ::R) where R <: Real
    return (-1)^m * sqrt(factorial(n + m) / factorial(n - m)) * B_mn_of_θ(m, n, θ) * exp(im * m * ϕ) # equation C.16
end

function C_mn_of_θ_ϕ(m::Int, n::Int, θ::R, ϕ::R) where R <: Real
    return (-1)^m * sqrt(factorial(n + m) / factorial(n - m)) * C_mn_of_θ(m, n, θ) * exp(im * m * ϕ) # equation C.17
end

function P_mn_of_θ_ϕ(m::Int, n::Int, θ::R, ϕ::R) where R <: Real
    return (-1)^m * sqrt(factorial(n + m) / factorial(n - m)) * P_mn_of_θ(m, n, θ) * exp(im * m * ϕ) # equation C.18
end

# calculate B(θ,ϕ), C(θ,ϕ), P(θ,ϕ), returns SVector
function B_mn_of_θ_ϕ_SVector(m::Int, n::Int, θ::R, ϕ::R) where R <: Real
    return (-1)^m * sqrt(factorial(n + m) / factorial(n - m)) * B_mn_of_θ_SVector(m, n, θ) * exp(im * m * ϕ) # equation C.16
end

function C_mn_of_θ_ϕ_SVector(m::Int, n::Int, θ::R, ϕ::R) where R <: Real
    return (-1)^m * sqrt(factorial(n + m) / factorial(n - m)) * C_mn_of_θ_SVector(m, n, θ) * exp(im * m * ϕ) # equation C.17
end

function P_mn_of_θ_ϕ_SVector(m::Int, n::Int, θ::R, ϕ::R) where R <: Real
    return (-1)^m * sqrt(factorial(n + m) / factorial(n - m)) * P_mn_of_θ_SVector(m, n, θ) * exp(im * m * ϕ) # equation C.18
end


#############################################################################################
# calculate spherical Bessel and Hankel functions and their derivarive
#function spherical_Bessel_j_n(n::Int, x::NN) where NN <: Number
    """
    Spherical Bessel function of the first kind.
    It can be calculated from ordinary Bessel function of the first kind "besselj" as in the code.
    check https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions
    """
#    return √(π / 2x) * besselj(n + 1 / 2, x)
#end

#function spherical_Bessel_y_n(n::Int, x::NN) where NN <: Number
    """
    Spherical Bessel function of the second kind.
    It can be calculated from ordinary Bessel function of the second kind "bessely" as in the code.
    check https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions
    """
#    return √(π / 2x) * bessely(n + 1 / 2, x)
#end

# overriding Bessel functions, not using SpecialFunctions
include("bessel.jl")
spherical_Bessel_j_n = spherical_Bessel_j_n_ELZOUKA
spherical_Bessel_y_n = spherical_Bessel_y_n_ELZOUKA

function spherical_Hankel_h1_n(n::Int, x::NN) where NN <: Number
    """
    Spherical Hankel function of the first kind. It can be calculated from spherical Bessel functions of the first and second kinds as in the code:
    """
    return spherical_Bessel_j_n(n, x) + im * spherical_Bessel_y_n(n, x)
end

function one_over_x_by_∂_x_j_n_by_∂x(n::Int, x::NN) where NN <: Number
    """
    Derivative of (spherical Bessel of first kind * x) divided by x
    """
    return (spherical_Bessel_j_n(n - 1, x) - n / x * spherical_Bessel_j_n(n, x))
end

function one_over_x_by_∂_x_y_n_by_∂x(n::Int, x::NN) where NN <: Number
    """
    Derivative of (spherical Bessel of second kind * x) divided by x
    """
    return (spherical_Bessel_y_n(n - 1, x) - n / x * spherical_Bessel_y_n(n, x))
end

function one_over_x_by_∂_x_h_n_by_∂x(n::Int, x::NN) where NN <: Number
    """
    Derivative of (spherical Hankel of first kind * x) divided by x
    """
    return one_over_x_by_∂_x_j_n_by_∂x(n, x) + im * one_over_x_by_∂_x_y_n_by_∂x(n, x)
end


#############################################################################################
# calculate (Rg)M(kr,θ,ϕ), (Rg)N(kr,θ,ϕ), returns Array
"""
    Parameters
    ==========
    kind: string, either ["regular" or "incoming"] or ["irregular" or "outgoing"]
"""
function M_mn_wave(m::Int, n::Int, kr::NN, θ::R, ϕ::R; kind="regular") where {R <: Real,NN <: Number}    
    radial_function, _ = get_radial_function_and_special_derivative_given_kind(kind)
    return γ_mn(m, n) * radial_function(n, kr) * C_mn_of_θ_ϕ(m, n, θ, ϕ)
end


"""
    Parameters
    ==========
    kind: string, either ["regular" or "incoming"] or ["irregular" or "outgoing"]
"""
function N_mn_wave(m::Int, n::Int, kr::NN, θ::R, ϕ::R; kind="regular") where {R <: Real,NN <: Number}    
    radial_function, radial_function_special_derivative  = get_radial_function_and_special_derivative_given_kind(kind)
    return γ_mn(m, n) * (
        n * (n + 1) / kr * radial_function(n, kr) * P_mn_of_θ_ϕ(m, n, θ, ϕ)
        + (radial_function_special_derivative(n, kr) * B_mn_of_θ_ϕ(m, n, θ, ϕ))
    )
end

# calculate (Rg)M(kr,θ,ϕ), (Rg)N(kr,θ,ϕ), returns SVector
"""
    Parameters
    ==========
    kind: string, either ["regular" or "incoming"] or ["irregular" or "outgoing"]
"""
function M_mn_wave_SVector(m::Int, n::Int, kr::NN, θ::R, ϕ::R; kind="regular") where {R <: Real,NN <: Number}    
    radial_function, _ = get_radial_function_and_special_derivative_given_kind(kind)
    return convert.(typeof(Complex(θ,θ)),
        γ_mn(m, n) * radial_function(n, kr) * C_mn_of_θ_ϕ_SVector(m, n, θ, ϕ)
    ) # make sure the output is of the same type as the input. # TODO: find a better way
end


"""
    Parameters
    ==========
    kind: string, either ["regular" or "incoming"] or ["irregular" or "outgoing"]
"""
function N_mn_wave_SVector(m::Int, n::Int, kr::NN, θ::R, ϕ::R; kind="regular") where {R <: Real,NN <: Number}    
    radial_function, radial_function_special_derivative  = get_radial_function_and_special_derivative_given_kind(kind)
    return convert.(typeof(Complex(θ,θ)),
            γ_mn(m, n) * (
            n * (n + 1) / kr * radial_function(n, kr) * P_mn_of_θ_ϕ_SVector(m, n, θ, ϕ)
            + (radial_function_special_derivative(n, kr) * B_mn_of_θ_ϕ_SVector(m, n, θ, ϕ))
        )
    )
end
