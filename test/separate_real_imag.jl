# TODO:
# 1- I need jacobian and gradient return Zero for any argument of type "Int"

using LinearAlgebra
using SpecialFunctions
using VectorSphericalWaves
using Zygote
using ChainRulesCore

import FiniteDifferences

#############################################################################################
"""
    Multiplication of two vectors Z1 and Z2.
We assume the first and second columns of each vector are the real and imaginary parts, respectively
Each of Z1 and Z2 has two columns, and arbitrary number of rows.
Multiplication algorithm for complex numbers can be found here: https://en.wikipedia.org/wiki/Multiplication_algorithm
"""
function complex_multiply(Z1, Z2)
    return vcat(complex_multiply.(Z1[:,1], Z1[:,2], Z2[:,1], Z2[:,2])...)
end
"""
    Multiplication of two vectors Z1=a+im*b and Z2=c+im*d.
a, b, c, and d are all real numbers or arrays of real numbers. All shoud have the same size.
Multiplication algorithm for complex numbers can be found here: https://en.wikipedia.org/wiki/Multiplication_algorithm
"""
function complex_multiply(a, b, c, d)
    # return hcat(a * c - b * d, b * c + a * d)
    return hcat(a .* c - b .* d, b .* c + a .* d)
end

"""
    Division of two vectors Z1 and Z2.
We assume the first and second columns of each vector are the real and imaginary parts, respectively
Each of Z1 and Z2 has two columns, and arbitrary number of rows.
Division algorithm for complex numbers can be found here: https://en.wikipedia.org/wiki/Multiplication_algorithm
"""
function complex_divide(Z1, Z2)
    return vcat(complex_divide.(Z1[:,1], Z1[:,2], Z2[:,1], Z2[:,2])...)
end
"""
    Division of two vectors Z1=a+im*b and Z2=c+im*d.
a, b, c, and d are all real numbers or arrays of real numbers. All shoud have the same size.
Division algorithm for complex numbers can be found here: https://en.wikipedia.org/wiki/Multiplication_algorithm
"""
function complex_divide(a, b, c, d)    
    return hcat((a * c + b * d) / (c^2 + d^2), (b * c - a * d) / (c^2 + d^2))
end


m, n, kr_r, kr_i, θ, ϕ = 1, 1, 1e5, 0, 0.2, 0.3
kr = complex(kr_r, kr_i)

gradient(VectorSphericalWaves.πₘₙ, m, n, θ)

#############################################################################################
# calculate B(θ), C(θ), P(θ)
function B_mn_of_θ_SeparateRealImag(m::Int, n::Int, θ::Real)
    # equation C.19
    B_real = vcat(0, VectorSphericalWaves.τₘₙ(m, n, θ), 0)
    B_imag = vcat(0, 0, VectorSphericalWaves.πₘₙ(m, n, θ))
    return hcat(B_real, B_imag)
end
jacobian(B_mn_of_θ_SeparateRealImag, m, n, θ)

function C_mn_of_θ_SeparateRealImag(m::Int, n::Int, θ::Real)
    # equation C.20  
    C_real = vcat(0, 0, -1 * VectorSphericalWaves.τₘₙ(m, n, θ))
    C_imag = vcat(0, VectorSphericalWaves.πₘₙ(m, n, θ), 0)
    return hcat(C_real, C_imag)
end
jacobian(C_mn_of_θ_SeparateRealImag, m, n, θ)


function P_mn_of_θ_SeparateRealImag(m::Int, n::Int, θ::Real)
    # equation C.21
    return hcat(P_mn_of_θ(m, n, θ), vcat(0, 0, 0))
end
jacobian(P_mn_of_θ_SeparateRealImag, m, n, θ)

jacobian(VectorSphericalWaves.P_mn_of_θ, m, n, θ) # and this works too

function convert_from_fun_of_θ_to_fun_of_θ_ϕ(fun_tobe_converted::Function, m::Int, n::Int, θ::Real, ϕ::Real)
    # equation C.16, C.17, and C.18
    coeff = (-1)^m * sqrt(factorial(n + m) / factorial(n - m))
    B_of_θ_coef = coeff .* fun_tobe_converted(m, n, θ)
    B_of_θ_coef_real = B_of_θ_coef[:,1]
    B_of_θ_coef_imag = B_of_θ_coef[:,2]
    exp_imϕ_real = cos(m * ϕ)
    exp_imϕ_imag = sin(m * ϕ)
    """
    # did't work
    B_real = B_of_θ_coef_real .* real(exp(im * m * ϕ)) - B_of_θ_coef_imag .* real(exp(im * m * ϕ))
    B_imag = B_of_θ_coef_real .* imag(exp(im * m * ϕ)) + B_of_θ_coef_imag .* real(exp(im * m * ϕ))
    """
    B_real = B_of_θ_coef_real .* exp_imϕ_real - B_of_θ_coef_imag .* exp_imϕ_imag
    B_imag = B_of_θ_coef_real .* exp_imϕ_imag + B_of_θ_coef_imag .* exp_imϕ_real

    return hcat(B_real, B_imag)
end

#############################################################################################
# calculate B(θ,ϕ), C(θ,ϕ), P(θ,ϕ)
function B_mn_of_θ_ϕ_SeparateRealImag(m::Int, n::Int, θ::Real, ϕ::Real)
    # equation C.16
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(B_mn_of_θ_SeparateRealImag, m, n, θ, ϕ)
end
jacobian(B_mn_of_θ_ϕ_SeparateRealImag, m, n, θ, ϕ)

function C_mn_of_θ_ϕ_SeparateRealImag(m::Int, n::Int, θ::Real, ϕ::Real)
    # equation C.17
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(C_mn_of_θ_SeparateRealImag, m, n, θ, ϕ)
end
jacobian(C_mn_of_θ_ϕ_SeparateRealImag, m, n, θ, ϕ)

function P_mn_of_θ_ϕ_SeparateRealImag(m::Int, n::Int, θ::Real, ϕ::Real)
    # equation C.18
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(P_mn_of_θ_SeparateRealImag, m, n, θ, ϕ)
end
jacobian(P_mn_of_θ_ϕ_SeparateRealImag, m, n, θ, ϕ)

#############################################################################################
# calculate spherical Bessel and Hankel functions and their derivarive
function SeparateRealImag_for_Bessel_Hankel_and_derivatives(fun_tobe_converted::Function, n::Int, x_r::Real, x_i::Real)
    bessel_complex = fun_tobe_converted(n, complex(x_r, x_i))
    return hcat(real(bessel_complex), imag(bessel_complex))
end

function spherical_Bessel_j_n_SeparateRealImag(n::Int, x_r::Real, x_i::Real)    
    return SeparateRealImag_for_Bessel_Hankel_and_derivatives(VectorSphericalWaves.spherical_Bessel_j_n, n, x_r, x_i)
end
jacobian(spherical_Bessel_j_n_SeparateRealImag, n, kr_r, kr_i)

function spherical_Bessel_y_n_SeparateRealImag(n::Int, x_r::Real, x_i::Real)    
    return SeparateRealImag_for_Bessel_Hankel_and_derivatives(VectorSphericalWaves.spherical_Bessel_y_n, n, x_r, x_i)
end
jacobian(spherical_Bessel_y_n_SeparateRealImag, n, kr_r, kr_i)

function spherical_Hankel_h1_n_SeparateRealImag(n::Int, x_r::Real, x_i::Real)
    return SeparateRealImag_for_Bessel_Hankel_and_derivatives(VectorSphericalWaves.spherical_Hankel_h1_n, n, x_r, x_i)
end
jacobian(spherical_Hankel_h1_n_SeparateRealImag, n, kr_r, kr_i)

function one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag(n::Int, x_r::Real, x_i::Real)
    # n_over_x = complex_divide(n * ones(size(x_r)), zeros(size(x_r)), x_r, x_i) # TODO: find a way to make sure all inputs have the same size
    n_over_x = complex_divide(n, 0, x_r, x_i)
    return (spherical_Bessel_j_n_SeparateRealImag(n - 1, x_r, x_i) - complex_multiply(n_over_x, spherical_Bessel_j_n_SeparateRealImag(n, x_r, x_i)))
end
jacobian(one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag, n, kr_r, kr_i)

function one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag(n::Int, x_r::Real, x_i::Real)
    # n_over_x = complex_divide(n * ones(size(x_r)), zeros(size(x_r)), x_r, x_i) # TODO: find a way to make sure all inputs have the same size
    n_over_x = complex_divide(n, 0, x_r, x_i)
    return (spherical_Bessel_y_n_SeparateRealImag(n - 1, x_r, x_i) - complex_multiply(n_over_x, spherical_Bessel_y_n_SeparateRealImag(n, x_r, x_i)))
end
jacobian(one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag, n, kr_r, kr_i)

function one_over_x_by_∂_x_h_n_by_∂x_SeparateRealImag(n::Int, x_r::Real, x_i::Real)
    der_j = one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag(n, x_r, x_i)
    der_y = one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag(n, x_r, x_i)
    return der_j + complex_multiply(0, 1, der_y[1], der_y[2])
    # return (one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag(n, x_r, x_i) + complex_multiply(0, 1, one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag(n, x_r, x_i)...))
end
jacobian(one_over_x_by_∂_x_h_n_by_∂x_SeparateRealImag, n, kr_r, kr_i)

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


function M_mn_wave_SeparateRealImag(m::Int, n::Int, kr_r::Real, kr_i::Real, θ::Real, ϕ::Real, kind="regular")
    # NOTE: I should remove all kwargs to make gradients and jacobian calculations smooth
    radial_function, _ = get_radial_function_and_special_derivative_given_kind_SeparateRealImag(kind)
    gamma_by_radial = γ_mn(m, n) .* radial_function(n, kr_r, kr_i)
    return complex_multiply(vcat(gamma_by_radial, gamma_by_radial, gamma_by_radial), C_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ))
    # TODO: write this in a more elegant way: vcat(gamma_by_radial, gamma_by_radial, gamma_by_radial)
end
jacobian(M_mn_wave_SeparateRealImag,m, n, kr_r, kr_i, θ, ϕ, "regular") # I can't add kwarg "kind". How can I add it?

# quick test
M_wave_calc_using_complex_numbers = M_mn_wave(m, n, complex(kr_r, kr_i), θ, ϕ; kind="regular")
M_mn_wave_SeparateRealImag(m, n, kr_r, kr_i, θ, ϕ; kind="regular") == hcat(real(M_wave_calc_using_complex_numbers), imag(M_wave_calc_using_complex_numbers))


function N_mn_wave_SeparateRealImag(m::Int, n::Int, kr_r::Real, kr_i::Real, θ::Real, ϕ::Real, kind="regular")
    radial_function, radial_function_special_derivative  = get_radial_function_and_special_derivative_given_kind_SeparateRealImag(kind)
    ceoff = complex_multiply(complex_divide(n * (n + 1), 0, kr_r, kr_i), radial_function(n, kr_r, kr_i))
    rad_fun_der = radial_function_special_derivative(n, kr_r, kr_i)
    return γ_mn(m, n) .* (
        complex_multiply([ceoff; ceoff; ceoff], P_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ))
        + complex_multiply([rad_fun_der; rad_fun_der; rad_fun_der], B_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ))
    )
    # TODO: write this in a more elegant way: [ceoff; ceoff; ceoff]
end
jacobian(N_mn_wave_SeparateRealImag,m, n, kr_r, kr_i, θ, ϕ, "regular") # I can't add kwarg "kind". How can I add it?
# quick test
N_wave_calc_using_complex_numbers = N_mn_wave(m, n, complex(kr_r, kr_i), θ, ϕ; kind="regular")
N_mn_wave_SeparateRealImag(m, n, kr_r, kr_i, θ, ϕ, "regular") == hcat(real(N_wave_calc_using_complex_numbers), imag(N_wave_calc_using_complex_numbers))