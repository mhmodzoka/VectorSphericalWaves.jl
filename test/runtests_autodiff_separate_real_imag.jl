# TODO:
# 1- I need Zygote.jacobian and gradient return Zero for any argument of type "Int"

# NOTE: Lessons learned:
# 1- I should remove all kwargs to make gradients and Zygote.jacobian calculations smooth

#############################################################################################
using LinearAlgebra
using SpecialFunctions
using ChainRulesCore
using ComplexOperations
using Traceur

import FiniteDifferences
import Zygote
import VectorSphericalWaves


#function fr
#############################################################################################
# testing automatic differenitation

m, n, kr_r, kr_i, θ, ϕ = 1, 1, 1e5, 1e2, 0.2, 0.3
s = 2 * n

# m, n, kr_r, kr_i, θ, ϕ = 1, 1, BigFloat(1e5), BigFloat(1e2), BigFloat(0.2), BigFloat(0.3)
kr = complex(kr_r, kr_i)

finite_difference_method = FiniteDifferences.central_fdm(5, 1; factor=1e-8);

Zygote.jacobian(VectorSphericalWaves.C_mn_of_θ_SeparateRealImag, m, n, θ)

# testing wignerdjmn
wignerdjmn_lean(θ) = VectorSphericalWaves.wignerdjmn(s, m, n, θ)
@time wignerdjmn_jacobian = Zygote.gradient(VectorSphericalWaves.wignerdjmn, s,m,n, θ)
@time wignerdjmn_lean_jacobian = Zygote.gradient(wignerdjmn_lean, θ)
@time wignerdjmn_lean_jacobian_FD = FiniteDifferences.jacobian(finite_difference_method, wignerdjmn_lean, θ)

# testing πₘₙ
m, n = -1, 1
πₘₙ_lean(θ) = VectorSphericalWaves.πₘₙ(m, n, θ)
@time πₘₙ_jacobian = Zygote.gradient(VectorSphericalWaves.πₘₙ,m, n, θ)
@time πₘₙ_lean_jacobian = Zygote.gradient(πₘₙ_lean, θ)
@time πₘₙ_lean_jacobian_FD = FiniteDifferences.jacobian(finite_difference_method, πₘₙ_lean, θ)

# testing τₘₙ
m, n = -1, 1
τₘₙ_lean(θ) = VectorSphericalWaves.τₘₙ(m, n, θ)
@time τₘₙ_jacobian = Zygote.gradient(VectorSphericalWaves.τₘₙ,m, n, θ)
@time τₘₙ_lean_jacobian = Zygote.gradient(τₘₙ_lean, θ)
@time τₘₙ_lean_jacobian_FD = FiniteDifferences.jacobian(finite_difference_method, τₘₙ_lean, θ)

n_max = 3
for n = 1:n_max
    for m = -n:n
        try
            πₘₙ_lean(θ) = VectorSphericalWaves.πₘₙ(m, n, θ)            
            πₘₙ_jacobian_FD = FiniteDifferences.jacobian(finite_difference_method, πₘₙ_lean, θ)

            y, back = Zygote.pullback(πₘₙ_lean, θ)
            πₘₙ_jacobian_pullback = back(θ)

            πₘₙ_lean_jacobian = Zygote.gradient(πₘₙ_lean, θ)
            πₘₙ_jacobian = Zygote.gradient(VectorSphericalWaves.πₘₙ, m, n, θ)
            println("n=$n, m=$m, πₘₙ_lean_jacobian=$πₘₙ_lean_jacobian, πₘₙ_jacobian=$(πₘₙ_jacobian[3]), πₘₙ_jacobian_pullback=$πₘₙ_jacobian_pullback")
        
            if abs(πₘₙ_jacobian_FD[1][1] - πₘₙ_lean_jacobian[1]) > 1e-5
                println("------- this causes large error: n=$n, m=$m")
            end
        catch
            println("------- this causes runtime error: n=$n, m=$m")
        end
    end
end

function ff(m::I, n::I, θ::R) where {R <: Real, I <: Integer}
    return m / sin(θ) * VectorSphericalWaves.wignerdjmn(n, 0, m, θ)
end
ChainRulesCore.@scalar_rule(
    ff(m::Integer, n::Integer, θ::Real),
    (       
        ChainRulesCore.ZeroTangent(),
        ChainRulesCore.ZeroTangent(),
        Zygote.gradient(
            (θ) -> ff(m,n,θ), θ
        )[1]
    )
)

######################################
# Zygote example to illustrate my point of Zygote ignoring differentiation with respect to Integers
import Zygote
# naiively defining my function, and using Zygote to calculate the gradient
f(n::Integer,x::Real) = sin(n*x)
Zygote.gradient(f, 100, π)
# Zygote performs automatic differentiation with respect to `n`, although `n` is an integer and `f` can't be differentiated with respect to `n`

# I hope Zygote can apply `Zygote.dropgrad()` to all arguments that are Integers, like `n`
f(n::Integer,x::Real) = sin(Zygote.dropgrad(n) * x)
Zygote.gradient(f, 100, π)
######################################

max_permissible_error = 1e-6
n_max = 7
for n = 1:n_max
    for m = -n:n     

        # TODO: try to use `isapprox`
        
        # πₘₙ, τₘₙ -------------------------
        for f in [
                VectorSphericalWaves.πₘₙ, VectorSphericalWaves.τₘₙ
            ]
            try
                error = abs(Zygote.gradient(f, m, n, θ)[3] - FiniteDifferences.jacobian(finite_difference_method, θ -> f(m, n, θ), θ)[1][1])                
                if error > max_permissible_error
                    println("error between Zygote and FiniteDifferences is high, error=$error at <<$f>> at n=$n, m=$m")
                end
            catch
                println("runtime error with <<$f>> at n=$n, m=$m")
            end
        end
        
        # B_mn_of_θ_SeparateRealImag, C_mn_of_θ_SeparateRealImag, P_mn_of_θ_SeparateRealImag ----------------------
        for f in [
                VectorSphericalWaves.B_mn_of_θ_SeparateRealImag,
                VectorSphericalWaves.C_mn_of_θ_SeparateRealImag,
                VectorSphericalWaves.P_mn_of_θ_SeparateRealImag
            ]
            try
                error = maximum(abs.(Zygote.jacobian(f, m, n, θ)[3] - FiniteDifferences.jacobian(finite_difference_method, θ -> f(m, n, θ), θ)[1]))
                if error > max_permissible_error
                    println("error between Zygote and FiniteDifferences is high, error=$error at <<$f>> at n=$n, m=$m")
                end
            catch
                println("runtime error with <<$f>> at n=$n, m=$m")
            end
        end

        # B_mn_of_θ_ϕ_SeparateRealImag, C_mn_of_θ_ϕ_SeparateRealImag, P_mn_of_θ_ϕ_SeparateRealImag ----------------------
        for f in [
                VectorSphericalWaves.B_mn_of_θ_ϕ_SeparateRealImag,
                VectorSphericalWaves.C_mn_of_θ_ϕ_SeparateRealImag,
                VectorSphericalWaves.P_mn_of_θ_ϕ_SeparateRealImag
            ]
            try
                error = maximum(abs.(Zygote.jacobian(f, m, n, θ, ϕ)[3] - FiniteDifferences.jacobian(finite_difference_method, (θ, ϕ) -> f(m, n, θ, ϕ), θ, ϕ)[1]))
                if error > max_permissible_error
                    println("error between Zygote and FiniteDifferences is high, error=$error at <<$f>> at n=$n, m=$m")
                end
            catch
                println("runtime error with <<$f>> at n=$n, m=$m")
            end
        end

        # B_mn_of_θ_ϕ_SeparateRealImag, C_mn_of_θ_ϕ_SeparateRealImag, P_mn_of_θ_ϕ_SeparateRealImag ----------------------
        for f in [
            VectorSphericalWaves.spherical_Bessel_j_n_SeparateRealImag,
            VectorSphericalWaves.spherical_Bessel_y_n_SeparateRealImag,
            VectorSphericalWaves.spherical_Hankel_h1_n_SeparateRealImag,
            VectorSphericalWaves.one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag,
            VectorSphericalWaves.one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag,
            VectorSphericalWaves.one_over_x_by_∂_x_h_n_by_∂x_SeparateRealImag
        ]
        try
            error = maximum(abs.(vcat(Zygote.jacobian(f, n, kr_r, kr_i)[2:3]...) - vcat(FiniteDifferences.jacobian(finite_difference_method, (kr_r, kr_i) -> f(n, kr_r, kr_i), θ, ϕ)...)))
            if error > max_permissible_error
                println("error between Zygote and FiniteDifferences is high, error=$error at <<$f>> at n=$n, m=$m")
            end
        catch
            println("runtime error with <<$f>> at n=$n, m=$m")
        end
    end
    end
end

f = VectorSphericalWaves.spherical_Bessel_j_n_SeparateRealImag
f(n, kr_r, kr_i)


Zygote.gradient(VectorSphericalWaves.πₘₙ, m, n, θ)
Zygote.gradient(VectorSphericalWaves.τₘₙ, m, n, θ)
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



function B_mn_of_θ_SeparateRealImag(m::Int, n::Int, θ::R) where R <: Real
    # equation C.19
    B_real = vcat(0, VectorSphericalWaves.τₘₙ(m, n, θ), 0)
    B_imag = vcat(0, 0, VectorSphericalWaves.πₘₙ(m, n, θ))
    return hcat(B_real, B_imag)
end

m, n, t_ = 5, 7, LinRange(0, π, 10000);
@time B_mn_of_θ_SeparateRealImag(m, n, t_[1]);
@time B_mn_of_θ_SeparateRealImag.(m, n, t_);
@time B_mn_of_θ_SeparateRealImag.(m, n, t_);


function B_mn_of_θ_SeparateRealImag(m::Int, n::Int, θ::R) where R <: Real
    # equation C.19    
    return hcat(vcat(0, VectorSphericalWaves.τₘₙ(m, n, θ), 0), vcat(0, 0, VectorSphericalWaves.πₘₙ(m, n, θ)))
end

m, n, t_ = 5, 7, LinRange(0, π, 10000);
@time B_mn_of_θ_SeparateRealImag(m, n, t_[1]);
@time B_mn_of_θ_SeparateRealImag.(m, n, t_);
@time B_mn_of_θ_SeparateRealImag.(m, n, t_);


