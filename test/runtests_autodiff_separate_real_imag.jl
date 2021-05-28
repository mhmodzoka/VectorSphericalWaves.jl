# TODO:
# 1- I need Zygote.jacobian and gradient return Zero for any argument of type "Int"

# NOTE: Lessons learned:
# 1- I should remove all kwargs to make gradients and Zygote.jacobian calculations smooth

#############################################################################################
using LinearAlgebra
using SpecialFunctions
using VectorSphericalWaves
using ChainRulesCore
using ComplexOperations
using Traceur

import FiniteDifferences
import Zygote

#############################################################################################
# testing automatic differenitation
m, n, kr_r, kr_i, θ, ϕ = 1, 1, 1e5, 1e2, 0.2, 0.3
m, n, kr_r, kr_i, θ, ϕ = 1, 1, BigFloat(1e5), BigFloat(1e2), BigFloat(0.2), BigFloat(0.3)
kr = complex(kr_r, kr_i)

Zygote.gradient(VectorSphericalWaves.πₘₙ, m, n, θ)
Zygote.jacobian(VectorSphericalWaves.B_mn_of_θ_SeparateRealImag, m, n, θ)
Zygote.jacobian(VectorSphericalWaves.C_mn_of_θ_SeparateRealImag, m, n, θ)
Zygote.jacobian(VectorSphericalWaves.P_mn_of_θ_SeparateRealImag, m, n, θ)
Zygote.jacobian(VectorSphericalWaves.P_mn_of_θ, m, n, θ) # and this works too

Zygote.jacobian(VectorSphericalWaves.B_mn_of_θ_ϕ_SeparateRealImag, m, n, θ, ϕ)
Zygote.jacobian(VectorSphericalWaves.C_mn_of_θ_ϕ_SeparateRealImag, m, n, θ, ϕ)
Zygote.jacobian(VectorSphericalWaves.P_mn_of_θ_ϕ_SeparateRealImag, m, n, θ, ϕ)

Zygote.jacobian(VectorSphericalWaves.spherical_Bessel_j_n_SeparateRealImag, n, kr_r, kr_i)

Zygote.jacobian(VectorSphericalWaves.spherical_Bessel_y_n_SeparateRealImag, n, kr_r, kr_i)

Zygote.jacobian(VectorSphericalWaves.spherical_Hankel_h1_n_SeparateRealImag, n, kr_r, kr_i)

Zygote.jacobian(VectorSphericalWaves.one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag, n, kr_r, kr_i)

Zygote.jacobian(VectorSphericalWaves.one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag, n, kr_r, kr_i)

Zygote.jacobian(VectorSphericalWaves.one_over_x_by_∂_x_h_n_by_∂x_SeparateRealImag, n, kr_r, kr_i)

Zygote.jacobian(VectorSphericalWaves.M_mn_wave_SeparateRealImag,m, n, kr_r, kr_i, θ, ϕ, "regular") # I can't add kwarg "kind". How can I add it?

M_wave_calc_using_complex_numbers = M_mn_wave(m, n, complex(kr_r, kr_i), θ, ϕ; kind="regular")
VectorSphericalWaves.M_mn_wave_SeparateRealImag(m, n, kr_r, kr_i, θ, ϕ, "regular") == hcat(real(M_wave_calc_using_complex_numbers), imag(M_wave_calc_using_complex_numbers))
Zygote.jacobian(VectorSphericalWaves.M_mn_wave_SeparateRealImag,m, n, kr_r, kr_i, θ, ϕ, "regular") # I can't add kwarg "kind". How can I add it?

N_wave_calc_using_complex_numbers = N_mn_wave(m, n, complex(kr_r, kr_i), θ, ϕ; kind="regular")
VectorSphericalWaves.N_mn_wave_SeparateRealImag(m, n, kr_r, kr_i, θ, ϕ, "regular") == hcat(real(N_wave_calc_using_complex_numbers), imag(N_wave_calc_using_complex_numbers))
Zygote.jacobian(VectorSphericalWaves.N_mn_wave_SeparateRealImag,m, n, kr_r, kr_i, θ, ϕ, "regular") # I can't add kwarg "kind". How can I add it?
