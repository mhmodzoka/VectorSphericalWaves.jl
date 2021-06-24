# TODO: I cound that "utils.jl" has been `include`ed in the main `VectorSphericalWaves` module. Is there a more neat way? when I `include("utils.jl")` here, I get lots of warnings.

using StaticArrays

# exporting main functions, separating real and imaginary parts
export M_mn_wave_SeparateRealImag
export N_mn_wave_SeparateRealImag

export M_mn_wave_SeparateRealImag_SMatrix
export N_mn_wave_SeparateRealImag_SMatrix

#############################################################################################
# The following functions are the same as before, but are desinged to be automatically differentiable (e.g., with Zygote).
# The main difference is that we don't use any complex numbers here. Complex numbers and vectors are represented as 1x2 and nx2 Matrices, respectively. The first and second columns represent the real and imaginary parts.
#############################################################################################
#############################################################################################
# calculate B(θ), C(θ), P(θ)
"""
    equation C.19, returns an array
It returns a 3x2 Matrix, the first and second columns are the real and imaginary parts, respectively. 
"""
function B_mn_of_θ_SeparateRealImag(m::I, n::I, θ::R) where {R <: Real, I <: Integer}
    # equation C.19    
    return hcat(
        vcat(0, τₘₙ(m, n, θ), 0),
        vcat(0, 0, πₘₙ(m, n, θ))
    )
end

"""
    Same as B_mn_of_θ_SeparateRealImag, but return SMatrix
It returns a 3x2 SMatrix, the first and second columns are the real and imaginary parts, respectively.
"""
function B_mn_of_θ_SeparateRealImag_SMatrix(m::I, n::I, θ::R) where {R <: Real, I <: Integer}
    # equation C.19    
    return hcat(
        SVector(0, τₘₙ(m, n, θ), 0),
        SVector(0, 0, πₘₙ(m, n, θ))
    )
end

"""
    equation C.20, returns an array
It returns a 3x2 Matrix, the first and second columns are the real and imaginary parts, respectively.
"""
function C_mn_of_θ_SeparateRealImag(m::I, n::I, θ::R) where {R <: Real, I <: Integer}
    # equation C.20    
    return hcat(
        vcat(0, 0, -1 * τₘₙ(m, n, θ)),
        vcat(0, πₘₙ(m, n, θ), 0)
    )
end

"""
    Same as C_mn_of_θ_SeparateRealImag, but return SMatrix
It returns a 3x2 SMatrix, the first and second columns are the real and imaginary parts, respectively.
"""
function C_mn_of_θ_SeparateRealImag_SMatrix(m::I, n::I, θ::R) where {R <: Real, I <: Integer}
    # equation C.20    
    return hcat(
        SVector(0, 0, -1 * τₘₙ(m, n, θ)),
        SVector(0, πₘₙ(m, n, θ), 0)
)
end

"""
    equation C.21, returns an array
It returns a 3x2 Matrix, the first and second columns are the real and imaginary parts, respectively.
"""
function P_mn_of_θ_SeparateRealImag(m::I, n::I, θ::R) where {R <: Real, I <: Integer}
    # equation C.21
    return hcat(
        P_mn_of_θ(m, n, θ),
        vcat(0, 0, 0)
)
end

"""
    Same as P_mn_of_θ_SeparateRealImag, returns SMatrix
It returns a 3x2 SMatrix, the first and second columns are the real and imaginary parts, respectively.
"""
function P_mn_of_θ_SeparateRealImag_SMatrix(m::I, n::I, θ::R) where {R <: Real, I <: Integer}
    # equation C.21
    return hcat(
        P_mn_of_θ_SVector(m, n, θ),
        SVector(zero(θ), zero(θ), zero(θ))
    )
end

function convert_from_fun_of_θ_to_fun_of_θ_ϕ(fun_tobe_converted::Function, m::I, n::I, θ::R, ϕ::R) where R <: Real
    # equation C.16, C.17, and C.18
    coeff = (-1)^m * sqrt_factorial_n_plus_m_over_factorial_n_minus_m(m, n)
    B_of_θ_coef = convert(R, coeff) .* fun_tobe_converted(m, n, θ)
    B_of_θ_coef_real = B_of_θ_coef[:,1]
    B_of_θ_coef_imag = B_of_θ_coef[:,2]
    exp_imϕ_real = cos(m * ϕ)
    exp_imϕ_imag = sin(m * ϕ)
    
    return hcat(
        B_of_θ_coef_real .* exp_imϕ_real - B_of_θ_coef_imag .* exp_imϕ_imag,
        B_of_θ_coef_real .* exp_imϕ_imag + B_of_θ_coef_imag .* exp_imϕ_real
    )
end

#############################################################################################
# calculate B(θ,ϕ), C(θ,ϕ), P(θ,ϕ)
function B_mn_of_θ_ϕ_SeparateRealImag(m::I, n::I, θ::R, ϕ::R) where {R <: Real, I <: Integer}
    # equation C.16
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(B_mn_of_θ_SeparateRealImag, m, n, θ, ϕ)
end

function C_mn_of_θ_ϕ_SeparateRealImag(m::I, n::I, θ::R, ϕ::R) where {R <: Real, I <: Integer}
    # equation C.17
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(C_mn_of_θ_SeparateRealImag, m, n, θ, ϕ)
end

function P_mn_of_θ_ϕ_SeparateRealImag(m::I, n::I, θ::R, ϕ::R) where {R <: Real, I <: Integer}
    # equation C.18
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(P_mn_of_θ_SeparateRealImag, m, n, θ, ϕ)
end

# same as above, but returns *_SMatrix
function B_mn_of_θ_ϕ_SeparateRealImag_SMatrix(m::I, n::I, θ::R, ϕ::R) where {R <: Real, I <: Integer}
    # equation C.16
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(B_mn_of_θ_SeparateRealImag_SMatrix, m, n, θ, ϕ)
end

function C_mn_of_θ_ϕ_SeparateRealImag_SMatrix(m::I, n::I, θ::R, ϕ::R) where {R <: Real, I <: Integer}
    # equation C.17
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(C_mn_of_θ_SeparateRealImag_SMatrix, m, n, θ, ϕ)
end

function P_mn_of_θ_ϕ_SeparateRealImag_SMatrix(m::I, n::I, θ::R, ϕ::R) where {R <: Real, I <: Integer}
    # equation C.18
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(P_mn_of_θ_SeparateRealImag_SMatrix, m, n, θ, ϕ)
end

#############################################################################################
# calculate spherical Bessel and Hankel functions and their derivarive
function SeparateRealImag_for_Bessel_Hankel_and_derivatives(fun_tobe_converted::Function, n::I, x_r::R, x_i::R) where {R <: Real, I <: Integer}
    bessel_complex = fun_tobe_converted(n, complex(x_r, x_i))
    return hcat(real(bessel_complex), imag(bessel_complex))
end

function spherical_Bessel_j_n_SeparateRealImag(n::I, x_r::R, x_i::R) where {R <: Real, I <: Integer}
    return SeparateRealImag_for_Bessel_Hankel_and_derivatives(spherical_Bessel_j_n, n, x_r, x_i)
end

function spherical_Bessel_y_n_SeparateRealImag(n::I, x_r::R, x_i::R) where {R <: Real, I <: Integer}
    return SeparateRealImag_for_Bessel_Hankel_and_derivatives(spherical_Bessel_y_n, n, x_r, x_i)
end

function spherical_Hankel_h1_n_SeparateRealImag(n::I, x_r::R, x_i::R) where {R <: Real, I <: Integer}
    return SeparateRealImag_for_Bessel_Hankel_and_derivatives(spherical_Hankel_h1_n, n, x_r, x_i)
end

function one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag(n::I, x_r::R, x_i::R) where {R <: Real, I <: Integer}
    # n_over_x = complex_divide(n * ones(size(x_r)), zeros(size(x_r)), x_r, x_i) # TODO: find a way to make sure all inputs have the same size
    n_over_x = complex_divide(n, zero(x_r), x_r, x_i)
    return (spherical_Bessel_j_n_SeparateRealImag(n - 1, x_r, x_i) - complex_multiply(n_over_x, spherical_Bessel_j_n_SeparateRealImag(n, x_r, x_i)))
end

function one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag(n::I, x_r::R, x_i::R) where {R <: Real, I <: Integer}
    # n_over_x = complex_divide(n * ones(size(x_r)), zeros(size(x_r)), x_r, x_i) # TODO: find a way to make sure all inputs have the same size
    n_over_x = complex_divide(n, zero(x_r), x_r, x_i)
    return (spherical_Bessel_y_n_SeparateRealImag(n - 1, x_r, x_i) - complex_multiply(n_over_x, spherical_Bessel_y_n_SeparateRealImag(n, x_r, x_i)))
end

function one_over_x_by_∂_x_h_n_by_∂x_SeparateRealImag(n::I, x_r::R, x_i::R) where {R <: Real, I <: Integer}
    der_j = one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag(n, x_r, x_i)
    der_y = one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag(n, x_r, x_i)
    return der_j + complex_multiply(0, 1, der_y[1], der_y[2])
    # return (one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag(n, x_r, x_i) + complex_multiply(0, 1, one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag(n, x_r, x_i)...))
end

function get_radial_function_and_special_derivative_given_kind_SeparateRealImag(kind)
    """
    return the radial function that is appropriate with the type of VSWF
    """

    if kind in ["regular", "incoming", 1]
        radial_function = spherical_Bessel_j_n_SeparateRealImag
        radial_function_special_derivative = one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag
    elseif kind in ["irregular", "outgoing", 2]
        radial_function = spherical_Hankel_h1_n_SeparateRealImag
        radial_function_special_derivative = one_over_x_by_∂_x_h_n_by_∂x_SeparateRealImag
    else
        throw(DomainError(""" 'kind' has to be one of the following: ["regular", "incoming", "irregular", "outgoing"] """))
    end

    return radial_function, radial_function_special_derivative
end

function M_mn_wave_SeparateRealImag(m::I, n::I, kr_r::R, kr_i::R, θ::R, ϕ::R, kind="regular") where {R <: Real, I <: Integer}
    radial_function, _ = get_radial_function_and_special_derivative_given_kind_SeparateRealImag(kind)
    gamma_by_radial = convert(R, γ_mn(m, n)) .* radial_function(n, kr_r, kr_i)
    return complex_multiply(vcat(gamma_by_radial, gamma_by_radial, gamma_by_radial), C_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ))
    # TODO: write this in a more elegant way: vcat(gamma_by_radial, gamma_by_radial, gamma_by_radial)
end

"""
    Same as `M_mn_wave_SeparateRealImag`, but returns SMatrix
"""
function M_mn_wave_SeparateRealImag_SMatrix(m::I, n::I, kr_r::R, kr_i::R, θ::R, ϕ::R, kind="regular") where {R <: Real, I <: Integer}
    radial_function, _ = get_radial_function_and_special_derivative_given_kind_SeparateRealImag(kind)
    gamma_by_radial = convert(R, γ_mn(m, n)) .* radial_function(n, kr_r, kr_i)
    return complex_multiply_SMatrix(gamma_by_radial .* SMatrix{3,2}(ones(3, 2)), C_mn_of_θ_ϕ_SeparateRealImag_SMatrix(m, n, θ, ϕ))
    # TODO: write this in a more elegant way: vcat(gamma_by_radial, gamma_by_radial, gamma_by_radial)
end

function N_mn_wave_SeparateRealImag(m::I, n::I, kr_r::R, kr_i::R, θ::R, ϕ::R, kind="regular") where {R <: Real, I <: Integer}
    radial_function, radial_function_special_derivative  = get_radial_function_and_special_derivative_given_kind_SeparateRealImag(kind)
    ceoff = complex_multiply(complex_divide(n * (n + 1), zero(kr_r), kr_r, kr_i), radial_function(n, kr_r, kr_i))
    rad_fun_der = radial_function_special_derivative(n, kr_r, kr_i)
    return convert(R, γ_mn(m, n)) .* (
        complex_multiply([ceoff; ceoff; ceoff], P_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ))
        + complex_multiply([rad_fun_der; rad_fun_der; rad_fun_der], B_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ))
    )
    # TODO: write this in a more elegant way: [ceoff; ceoff; ceoff]
end

"""
    Same as `N_mn_wave_SeparateRealImag`, but returns SMatrix
"""
function N_mn_wave_SeparateRealImag_SMatrix(m::I, n::I, kr_r::R, kr_i::R, θ::R, ϕ::R, kind="regular") where {R <: Real, I <: Integer}
    radial_function, radial_function_special_derivative  = get_radial_function_and_special_derivative_given_kind_SeparateRealImag(kind)
    ceoff = complex_multiply(complex_divide(n * (n + 1), zero(kr_r), kr_r, kr_i), radial_function(n, kr_r, kr_i))
    rad_fun_der = radial_function_special_derivative(n, kr_r, kr_i)
    return convert(R, γ_mn(m, n)) .* (
        complex_multiply_SMatrix(ceoff .* SMatrix{3,2}(ones(3, 2)), P_mn_of_θ_ϕ_SeparateRealImag_SMatrix(m, n, θ, ϕ))
        + complex_multiply_SMatrix(rad_fun_der .* SMatrix{3,2}(ones(3, 2)), B_mn_of_θ_ϕ_SeparateRealImag_SMatrix(m, n, θ, ϕ))
    )
    # TODO: write this in a more elegant way: [ceoff; ceoff; ceoff]
end
