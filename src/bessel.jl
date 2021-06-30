function besselj_ELZOUKA(n::Int, x::NN) where NN <: Number
end

function bessely_ELZOUKA(n::Int, x::NN) where NN <: Number
end

"""
Calculate Bessel functions using recurrence relation
"""
function recurrence_relation_bessel(function_here::Function, n::Int, x::NN) where NN <: Number
    return (2*(n-1)/x) * function_here(n-1, x) - function_here(n-2, x)
end

"""
    Recurrence formula for ALL spherical Bessel functions.
I mean by spherical Bessel functions the following:
jₙ, spherical Bessel of the first kind
yₙ, spherical Bessel of the second kind
hₙ⁽¹⁾
hₙ⁽²⁾
http://dlmf.nist.gov/10.51.E1
"""
function recurrence_relation_spherical_bessel(function_here::Function, n::Int, x::NN) where NN <: Number
    return ((2*(n-1)+1)/x) * function_here(n-1, x) - function_here(n-2, x)
end


function recurrence_relation_spherical_bessel(n::Int, x::NN, bessel_n_minus_1::NN, bessel_n_minus_2::NN) where NN <: Number
    return (2*(n-1)+1)/x * bessel_n_minus_1 - bessel_n_minus_2
end

"""
    Calculate spherical Bessel function, given the Bessel functions at n=0 and n=1
Works for ALL spherical Bessel functions
"""
function spherical_Bessel_n_ELZOUKA(n::Int, x::NN, bessel_0::NN, bessel_1::NN) where NN <: Number
    if n == 0
        return bessel_0
    elseif n == 1
        return bessel_1
    else
        # use the recurrence
        n_now = 2
        bessel_n_minus_2 = bessel_0
        bessel_n_minus_1 = bessel_1
        while n_now <= n
            bessel_n_now = recurrence_relation_spherical_bessel(n_now, x, bessel_n_minus_1, bessel_n_minus_2)
            if n_now == n
                return bessel_n_now
            else
                bessel_n_minus_2 = bessel_n_minus_1
                bessel_n_minus_1 = bessel_n_now
                n_now += 1
            end
        end
    end
end

# jₙ, spherical Bessel of the first kind =======================================
"""
    Spherical Bessel function of the first kind, at n=0
from http://dlmf.nist.gov/10.49.E3
"""
function spherical_Bessel_j_0_ELZOUKA(x::NN) where NN <: Number
    return sin(x)/x
end

"""
    Spherical Bessel function of the first kind, at n=1
from http://dlmf.nist.gov/10.49.E3
"""
function spherical_Bessel_j_1_ELZOUKA(x::NN) where NN <: Number
    return sin(x)/(x^2) - cos(x)/x
end

"""
    Spherical Bessel function of the first kind, using recurrence relation
"""
function spherical_Bessel_j_n_ELZOUKA(n::Int, x::NN) where NN <: Number
    if n == 0
        return spherical_Bessel_j_0_ELZOUKA(x)
    elseif n == 1
        return spherical_Bessel_j_1_ELZOUKA(x)
    else
        bessel_0 = spherical_Bessel_j_0_ELZOUKA(x)
        bessel_1 = spherical_Bessel_j_1_ELZOUKA(x)
        return spherical_Bessel_n_ELZOUKA(n, x, bessel_0, bessel_1)
    end
end


ChainRulesCore.@scalar_rule(
    spherical_Bessel_j_n_ELZOUKA(n::Integer, x::Real),
    (
        ChainRulesCore.ZeroTangent(),        
        (spherical_Bessel_j_n_ELZOUKA(n - 1, x) - spherical_Bessel_j_n_ELZOUKA(n + 1, x)) / 2,
    )
)


# yₙ, spherical Bessel of the second kind =======================================
"""
    Spherical Bessel function of the second kind, at n=0
from http://dlmf.nist.gov/10.49.E5
"""
function spherical_Bessel_y_0_ELZOUKA(x::NN) where NN <: Number
    return -cos(x)/x
end

"""
    Spherical Bessel function of the second kind, at n=1
from http://dlmf.nist.gov/10.49.E5
"""
function spherical_Bessel_y_1_ELZOUKA(x::NN) where NN <: Number
    return -cos(x)/(x^2) - sin(x)/x
end

"""
    Spherical Bessel function of the second kind, using recurrence relation
"""
function spherical_Bessel_y_n_ELZOUKA(n::Int, x::NN) where NN <: Number
    if n == 0
        return spherical_Bessel_y_0_ELZOUKA(x)
    elseif n == 1
        return spherical_Bessel_y_1_ELZOUKA(x)
    else
        bessel_0 = spherical_Bessel_y_0_ELZOUKA(x)
        bessel_1 = spherical_Bessel_y_1_ELZOUKA(x)
        return spherical_Bessel_n_ELZOUKA(n, x, bessel_0, bessel_1)
    end
end


ChainRulesCore.@scalar_rule(
    spherical_Bessel_y_n_ELZOUKA(n::Integer, x::Real),
    (
        ChainRulesCore.ZeroTangent(),        
        (spherical_Bessel_y_n_ELZOUKA(n - 1, x) - spherical_Bessel_y_n_ELZOUKA(n + 1, x)) / 2,
    )
)
