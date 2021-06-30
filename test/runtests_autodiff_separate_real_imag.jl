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
import Flux
import VectorSphericalWaves

function calculating_error(array_tobe_evaluated, array_groundtruth; max_permissible_error = 1e-6, max_permissible_relative_error = 1e-6)
    if 0 in array_groundtruth                    
        error = maximum(abs.(array_tobe_evaluated - array_groundtruth))       
        if error > max_permissible_error
            println("error between Zygote and manual_numerical_diff is high, error=$error at <<$f>> at n=$n")
        end               
        
    
    else
        relative_error = maximum(abs.((array_groundtruth - array_tobe_evaluated)/array_groundtruth))       
    
        if relative_error > max_permissible_relative_error
            println("relative error between Zygote and manual_numerical_diff is high, error=$relative_error at <<$f>> at n=$n")
        end
        
    end

end

function ∂wignerdjmn(s,m,n,θ)
    Flux.gradient(VectorSphericalWaves.wignerdjmn, s,m,n,θ)
end
function ∂M_mn_wave_SeparateRealImag(m,n, kr_r, kr_i, θ, ϕ, kind)
    Flux.jacobian(VectorSphericalWaves.M_mn_wave_SeparateRealImag, m,n, kr_r, kr_i, θ, ϕ, kind)
end
function ∂N_mn_wave_SeparateRealImag(m,n, kr_r, kr_i, θ, ϕ, kind)
    Flux.jacobian(VectorSphericalWaves.N_mn_wave_SeparateRealImag, m,n, kr_r, kr_i, θ, ϕ, kind)
end


m, n, kr_r, kr_i, θ, ϕ = 1, 1, 20., 1., 0.2, 0.3
s = 2 * n

# m, n, kr_r, kr_i, θ, ϕ = 1, 1, BigFloat(1e5), BigFloat(1e2), BigFloat(0.2), BigFloat(0.3)
kr = complex(kr_r, kr_i)

finite_difference_method = FiniteDifferences.central_fdm(5, 1; factor=1e-8);

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

n_max = 8
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


##################################################################

n_max = 20
epss = big(1e-30); 
max_permissible_error = 1e-6
max_permissible_relative_error = 1e-6

for n = 1:n_max
    # TODO: try to use `isapprox`        
    
    # Bessel functions
    for f in [
        VectorSphericalWaves.spherical_Bessel_j_n_SeparateRealImag,
        VectorSphericalWaves.spherical_Bessel_y_n_SeparateRealImag,
        VectorSphericalWaves.spherical_Hankel_h1_n_SeparateRealImag,
        VectorSphericalWaves.one_over_x_by_∂_x_j_n_by_∂x_SeparateRealImag,
        VectorSphericalWaves.one_over_x_by_∂_x_y_n_by_∂x_SeparateRealImag,
        VectorSphericalWaves.one_over_x_by_∂_x_h_n_by_∂x_SeparateRealImag
    ]
        try
            println("Zygote($f)")
            @time jacobian_Zygote = vcat(Zygote.jacobian(f, n, kr_r, kr_i)[2:3]...)
            println("manual_numerical_diff($f)")
            @time jacobian_manual_numerical_diff = hcat(
                (f(n, big(kr_r+epss), big(kr_i)) - f(n, big(kr_r-epss), big(kr_i)))/ (2*epss),
                (f(n, big(kr_r), big(kr_i+epss)) - f(n, big(kr_r), big(kr_i-epss)))/ (2*epss)
            )'
            jacobian_FD = vcat(FiniteDifferences.jacobian(finite_difference_method, (kr_r, kr_i) -> f(n, kr_r, kr_i), θ, ϕ)...)

            if 0 in jacobian_manual_numerical_diff            
                error_FD = maximum(abs.(jacobian_Zygote - jacobian_FD))
                error_manual_numerical_diff = maximum(abs.(jacobian_Zygote - jacobian_manual_numerical_diff))
                if min(error_FD, error_manual_numerical_diff) > max_permissible_error
                    if error_FD > max_permissible_error
                        println("error between Zygote and FiniteDifferences is high, error=$error_FD at <<$f>> at n=$n")
                    end
                    if error_manual_numerical_diff > max_permissible_error
                        println("error between Zygote and manual_numerical_diff is high, error=$error_manual_numerical_diff at <<$f>> at n=$n")
                    end                    
                end   
            
            else
                relative_error_Zygote = maximum(abs.((jacobian_manual_numerical_diff - jacobian_Zygote)/jacobian_manual_numerical_diff))
                relative_error_FD = maximum(abs.((jacobian_manual_numerical_diff- jacobian_FD)/jacobian_manual_numerical_diff))
                if min(relative_error_Zygote, relative_error_FD) > max_permissible_relative_error
                    if relative_error_Zygote > max_permissible_relative_error
                        println("relative error between Zygote and manual_numerical_diff is high, error=$relative_error_Zygote at <<$f>> at n=$n")
                    end
                    if relative_error_FD > max_permissible_relative_error
                        println("relative error between FiniteDifference and manual_numerical_diff is high, error=$relative_error_FD at <<$f>> at n=$n")
                    end                
                end            
            end
        catch
            println("runtime error with <<$f>> at n=$n")
        end
    end
end


for n = 1:n_max
    for m = -n:n
        # TODO: try to use `isapprox`
        
        # πₘₙ, τₘₙ -------------------------
        for f in [
                VectorSphericalWaves.πₘₙ,
                VectorSphericalWaves.τₘₙ
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

        # M and N
        for f in [
            VectorSphericalWaves.M_mn_wave_SeparateRealImag,
            VectorSphericalWaves.N_mn_wave_SeparateRealImag,            
        ]
            for kind in ["regular", "irregular"]
                try
                    println("Zygote($f, kind=$kind)")
                    @time jacobian_Zygote = vcat(Zygote.jacobian(f, m,n, kr_r, kr_i, θ, ϕ, kind)[3:6]...)
                    println("Flux($f, kind=$kind)")
                    if f == VectorSphericalWaves.M_mn_wave_SeparateRealImag
                        @time zz = ∂M_mn_wave_SeparateRealImag(m,n, kr_r, kr_i, θ, ϕ, kind)
                    elseif f == VectorSphericalWaves.N_mn_wave_SeparateRealImag
                        @time zz = ∂N_mn_wave_SeparateRealImag(m,n, kr_r, kr_i, θ, ϕ, kind)
                    end
                    println("manual_numerical_diff($f, kind=$kind)")
                    @time jacobian_manual_numerical_diff = vcat(
                        vcat((f(m,n, big(kr_r+epss), big(kr_i), big(θ), big(ϕ), kind) - f(m,n, big(kr_r-epss), big(kr_i), big(θ), big(ϕ), kind)) / (2*epss)...),
                        vcat((f(m,n, big(kr_r), big(kr_i+epss), big(θ), big(ϕ), kind) - f(m,n, big(kr_r), big(kr_i-epss), big(θ), big(ϕ), kind)) / (2*epss)...),
                        vcat((f(m,n, big(kr_r), big(kr_i), big(θ+epss), big(ϕ), kind) - f(m,n, big(kr_r), big(kr_i), big(θ-epss), big(ϕ), kind)) / (2*epss)...),
                        vcat((f(m,n, big(kr_r), big(kr_i), big(θ), big(ϕ+epss), kind) - f(m,n, big(kr_r), big(kr_i), big(θ), big(ϕ-epss), kind)) / (2*epss)...),
                    )
                    println("FiniteDifference($f, kind=$kind)")
                    @time jacobian_FD = vcat(FiniteDifferences.jacobian(finite_difference_method, (kr_r, kr_i, θ, ϕ) -> f(m,n, kr_r, kr_i, θ, ϕ, kind), kr_r, kr_i, θ, ϕ)...)

                    if 0 in jacobian_manual_numerical_diff                    
                        error = maximum(
                            abs.(jacobian_Zygote - jacobian_manual_numerical_diff)
                        )
                        if error > max_permissible_error
                            println("error between Zygote and manual_numerical_diff is high, error=$error at <<$f>>, kind=$kind at n=$n, m=$m")
                            println("here is [Zygote, manual_numerical_diff, abs(Zygote-manual_numerical_diff)]:")
                            display(hcat(jacobian_Zygote, jacobian_manual_numerical_diff, abs.(jacobian_Zygote - jacobian_manual_numerical_diff)))
                            println()
                            println()
                        end

                    else
                        relative_error = maximum(
                            abs.(
                                (jacobian_manual_numerical_diff - jacobian_Zygote)./jacobian_manual_numerical_diff
                            )
                        )
                        if relative_error > max_permissible_relative_error
                            println("relative error between Zygote and manual_numerical_diff is high, relative_error=$relative_error at <<$f>>, kind=$kind at n=$n, m=$m")
                            println("here is [Zygote, manual_numerical_diff, abs(Zygote-manual_numerical_diff)]:")
                            display(hcat(jacobian_Zygote, jacobian_manual_numerical_diff, (jacobian_manual_numerical_diff - jacobian_Zygote)./jacobian_manual_numerical_diff))
                            println()
                            println()
                        end

                    end
                catch
                    println("runtime error with <<$f>>, kind=$kind at n=$n, m=$m")
                end
            end
        end
    end
end