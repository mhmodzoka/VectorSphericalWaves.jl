using LinearAlgebra
using SpecialFunctions
using VectorSphericalWaves
using Zygote
using ChainRulesCore

#############################################################################################
"""
    Multiplication of two vectors Z1 and Z2.
We assume the first and second columns of each vector are the real and imaginary parts, respectively
https://en.wikipedia.org/wiki/Multiplication_algorithm
"""
function complex_multiply(Z1, Z2)
    return vcat(complex_multiply.(Z1[:,1], Z1[:,2], Z2[:,1], Z2[:,2])...)
end
"""
    Multiplication of two vectors Z1=a+im*b and Z2=c+im*d.
https://en.wikipedia.org/wiki/Multiplication_algorithm
"""
function complex_multiply(a, b, c, d)
    # return hcat(a * c - b * d, b * c + a * d)
    return hcat(a * c - b * d, b * c + a * d)
end

"""
    Division of two vectors Z1 and Z2.
We assume the first and second columns of each vector are the real and imaginary parts, respectively
https://en.wikipedia.org/wiki/Multiplication_algorithm
"""
function complex_divide(Z1, Z2)
    return vcat(complex_divide.(Z1[:,1], Z1[:,2], Z2[:,1], Z2[:,2])...)
end
"""
    Division of two vectors Z1=a+im*b and Z2=c+im*d.
https://en.wikipedia.org/wiki/Multiplication_algorithm
"""
function complex_divide(a, b, c, d)    
    return hcat((a * c + b * d) / (c^2 + d^2), (b * c - a * d) / (c^2 + d^2))
end


m, n, kr_r, kr_i, θ, ϕ = 1, 1, 1e5, 0, 0.2, 0.3
kr = complex(kr_r, kr_i)

gradient(VectorSphericalWaves.πₘₙ, m, n, θ)

#############################################################################################
# calculate B(θ), C(θ), P(θ)

# this works
function B_mn_of_θ_SeparateRealImag(m, n, θ)
    # equation C.19
    B_real = vcat(0, VectorSphericalWaves.τₘₙ(m, n, θ), 0)
    B_imag = vcat(0, 0, VectorSphericalWaves.πₘₙ(m, n, θ))
    return hcat(B_real, B_imag)
end
jacobian(B_mn_of_θ_SeparateRealImag, m, n, θ)

# this works
function C_mn_of_θ_SeparateRealImag(m, n, θ)  
    # equation C.20  
    C_real = vcat(0, 0, -1 * VectorSphericalWaves.τₘₙ(m, n, θ))
    C_imag = vcat(0, VectorSphericalWaves.πₘₙ(m, n, θ), 0)
    return hcat(C_real, C_imag)
end
jacobian(C_mn_of_θ_SeparateRealImag, m, n, θ)

# this works
function P_mn_of_θ_SeparateRealImag(m, n, θ)    
    # equation C.21
    return hcat(P_mn_of_θ(m, n, θ), vcat(0, 0, 0))
end
jacobian(P_mn_of_θ_SeparateRealImag, m, n, θ)

# and this works too
jacobian(VectorSphericalWaves.P_mn_of_θ, m, n, θ)

function convert_from_fun_of_θ_to_fun_of_θ_ϕ(fun_tobe_converted, m, n, θ, ϕ)
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
function B_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ)
    # equation C.16
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(B_mn_of_θ_SeparateRealImag, m, n, θ, ϕ)
end
jacobian(B_mn_of_θ_ϕ_SeparateRealImag, m, n, θ, ϕ)

function C_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ)
    # equation C.17
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(C_mn_of_θ_SeparateRealImag, m, n, θ, ϕ)
end
jacobian(C_mn_of_θ_ϕ_SeparateRealImag, m, n, θ, ϕ)

function P_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ)
    # equation C.18
    return convert_from_fun_of_θ_to_fun_of_θ_ϕ(P_mn_of_θ_SeparateRealImag, m, n, θ, ϕ)
end
jacobian(P_mn_of_θ_ϕ_SeparateRealImag, m, n, θ, ϕ)

#############################################################################################
# calculate spherical Bessel and Hankel functions and their derivarive
function SeparateRealImag_for_Bessel_Hankel_and_derivatives(fun_tobe_converted, n, x_r, x_i)
    bessel_complex = fun_tobe_converted(n, complex(x_r, x_i))
    return hcat(real(bessel_complex), imag(bessel_complex))
end

function spherical_Bessel_j_n_SeparateRealImag(n, x_r, x_i)    
    return SeparateRealImag_for_Bessel_Hankel_and_derivatives(VectorSphericalWaves.spherical_Bessel_j_n, n, x_r, x_i)
end
jacobian(spherical_Bessel_j_n_SeparateRealImag, n, kr_r, kr_i)

function spherical_Bessel_y_n_SeparateRealImag(n, x_r, x_i)    
    return SeparateRealImag_for_Bessel_Hankel_and_derivatives(VectorSphericalWaves.spherical_Bessel_y_n, n, x_r, x_i)
end
jacobian(spherical_Bessel_y_n_SeparateRealImag, n, kr_r, kr_i)

function spherical_Hankel_h1_n_SeparateRealImag(n, x_r, x_i)
    return SeparateRealImag_for_Bessel_Hankel_and_derivatives(VectorSphericalWaves.spherical_Hankel_h1_n, n, x_r, x_i)
end
jacobian(spherical_Hankel_h1_n_SeparateRealImag, n, kr_r, kr_i)

function one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag(n, x_r, x_i)
    # n_over_x = complex_divide(n * ones(size(x_r)), zeros(size(x_r)), x_r, x_i) # TODO: find a way to make sure all inputs have the same size
    n_over_x = complex_divide(n, 0, x_r, x_i)
    return (spherical_Bessel_j_n_SeparateRealImag(n - 1, x_r, x_i) - complex_multiply(n_over_x, spherical_Bessel_j_n_SeparateRealImag(n, x_r, x_i)))
end
jacobian(one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag, n, kr_r, kr_i)

function one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag(n, x_r, x_i)
    # n_over_x = complex_divide(n * ones(size(x_r)), zeros(size(x_r)), x_r, x_i) # TODO: find a way to make sure all inputs have the same size
    n_over_x = complex_divide(n, 0, x_r, x_i)
    return (spherical_Bessel_y_n_SeparateRealImag(n - 1, x_r, x_i) - complex_multiply(n_over_x, spherical_Bessel_y_n_SeparateRealImag(n, x_r, x_i)))
end
jacobian(one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag, n, kr_r, kr_i)

function one_over_x_by_∂_x_h_n_by_∂x_SeparateRealImag(n, x_r, x_i)
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

    if kind in ["regular", "incoming"]
        radial_function = spherical_Bessel_j_n_SeparateRealImag
        radial_function_special_derivative = one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag
    elseif kind in ["irregular", "outgoing"]
        radial_function = spherical_Hankel_h1_n_SeparateRealImag
        radial_function_special_derivative = one_over_x_by_∂_x_h_n_by_∂x_SeparateRealImag
    else
        throw(DomainError(""" 'kind' has to be one of the following: ["regular", "incoming", "irregular", "outgoing"] """))
    end

    return radial_function, radial_function_special_derivative
end


function M_mn_wave_SeparateRealImag(m, n, kr_r, kr_i, θ, ϕ; kind="regular")
    radial_function, _ = get_radial_function_and_special_derivative_given_kind_SeparateRealImag(kind)
    gamma_by_radial = γ_mn(m, n) .* radial_function(n, kr_r, kr_i)
    return complex_multiply(vcat(gamma_by_radial, gamma_by_radial, gamma_by_radial), C_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ))
end
jacobian(M_mn_wave_SeparateRealImag,m, n, kr_r, kr_i, θ, ϕ) # I can't add kwarg "kind". How can I add it?



"""
# this didn't work
function C_mn_of_θ_SeparateRealImag(m, n, θ) 
    C = vcat(
        0,                  # r-component
        im * VectorSphericalWaves.πₘₙ(m, n, θ), # θ-component
        -1 * VectorSphericalWaves.τₘₙ(m, n, θ),    # ϕ-component      
    ) # equation C.20   
    C_real = real(C)
    C_imag = imag(C)
    return hcat(C_real, C_imag)
end
jacobian(C_mn_of_θ_SeparateRealImag, m, n, θ)
"""

"""
# this didn't work
function C_mn_of_θ_SeparateRealImag(m, n, θ)
    C = VectorSphericalWaves.C_mn_of_θ(m, n, θ)
    return hcat(real(C), imag(C))    
end
jacobian(C_mn_of_θ_SeparateRealImag, m, n, θ)
"""



function C_mn_of_θ_ϕ_SeparateRealImag(m, n, θ, ϕ)
    C = VectorSphericalWaves.C_mn_of_θ_ϕ(m, n, θ, ϕ)
    return [real(C), imag(C)]
end
jacobian(C_mn_of_θ_ϕ_SeparateRealImag, m, n, θ, ϕ)


function M_mn_wave_SeparateRealImag(m, n, kr_r, kr_i, θ, ϕ; kind="regular")
    kr = complex(kr_r, kr_i)
    M = VectorSphericalWaves.M_mn_wave(m, n, kr, θ, ϕ; kind=kind)
    return [real(M), imag(M)]
end
m, n, kr_r, kr_i, θ, ϕ = 1, 1, 1e5, 0, 0.2, 0.3
jacobian(M_mn_wave_SeparateRealImag, m, n, kr_r, kr_i, θ, ϕ)






m, n, θ = 1, 2, 1.2
gradient(B_mn_of_θ_separate_realimage, m, n, θ)
ForwardDiff.jacobian(B_mn_of_θ_separate_realimage, m, n, θ)


VectorSphericalWaves.spherical_Bessel_j_n(1, 1.3 + im)

gradient(VectorSphericalWaves.spherical_Bessel_j_n, 1, 1.3)

function spherical_Bessel_j_n_realimag(n, x_r, x_i)
    """
    Spherical Bessel function of the first kind.
    It can be calculated from ordinary Bessel function of the first kind "besselj" as in the code.
    check https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions
    """
    x_complex = complex(x_r, x_i)
    z_complex = √(π / 2x_complex) * besselj(n + 1 / 2, x_complex)
    return [real(z_complex), imag(z_complex)]
end 
jacobian(spherical_Bessel_j_n_realimag, 1, 1.3, 0.1)



function spherical_Bessel_j_n_SeparateRealImag(n, x_r, x_i)
    bessel_complex = VectorSphericalWaves.spherical_Bessel_j_n(n, complex(x_r, x_i))
    return [real(bessel_complex), imag(bessel_complex)]
end
jacobian(spherical_Bessel_j_n_realimag, 1, 1.3, 0.1)


