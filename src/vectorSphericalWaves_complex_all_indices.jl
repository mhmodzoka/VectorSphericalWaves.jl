# TODO: I cound that "utils.jl" has been `include`ed in the main `VectorSphericalWaves` module. Is there a more neat way? when I `include("utils.jl")` here, I get lots of warnings.

# exporting functions that calculate for all m,n combinations
export M_N_wave_all_m_n
export B_C_P_mn_of_θ_ϕ_for_all_m_n
export πₘₙ_τₘₙ_all_all_m_n

#############################################################################################
# This version is different from "VectorSphericalHarmonics.jl", as it calculate VSWF for all m,n and all kr, θ, ϕ in one call
# This should be more stable and faster, as it employs recurrence relations in calculating πₘₙ and τₘₙ
# using recurrence relations in calculating πₘₙ and τₘₙ

#############################################################################################
# calculate π(θ) and τ(θ) using recurrence relations
function πₘₙ_τₘₙ_all_all_m_n(n_max::I, θ::R; verbose=false) where {R <: Real, I <: Integer}
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
function Legendre_polynomials_Pn_array(n_max::I, x::NN) where {I <: Integer, NN <: Number}
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

function Legendre_polynomials_Pn(n::I, x::NN) where {I <: Integer, NN <: Number}
    """
    Calculate Legendre polynomials Pₙ(x) at a given n
    """
    return Legendre_polynomials_Pn_array(n, x)[n]
end

function Associated_Legendre_polynomials_Pmn_array(m::I, n_max::I, x) where {I <: Integer}
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

function Associated_Legendre_polynomials_Pmn(m::I, n::I, x) where {I <: Integer}
    """
    Associated Legendre polynomials Pᵐₙ(x) at a given n
    """
    return Associated_Legendre_polynomials_Pmn_array(m, n, x)[n]
end

#############################################################################################
# Wigner-d, using recurrence
# TODO: I think this will not work for large s,m,n values, try to fix it as in `wignerdjmn_ELZOUKA`
function wignerd_and_∂wignerd_for_all_s(s_max::I, m::I, n::I, θ::R; get_derivatives=true, verbose=false) where {R <: Real, I <: Integer}
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

    # calculate d^s_min__m_n(θ) from eq. B.24
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
function πₘₙ_τₘₙ_all_all_m_n_using_wigner(n_max::I, θ::R; verbose=false) where {R <: Real, I <: Integer}
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

function B_C_mn_of_θ_for_all_m_n(n_max::I, θ::R) where {R <: Real, I <: Integer}
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

function P_mn_of_θ_for_all_m_n(n_max::I, θ::R) where {R <: Real, I <: Integer}
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

function B_C_P_mn_of_θ_ϕ_for_all_m_n(n_max::I, θ::R, ϕ::R) where {R <: Real, I <: Integer}
    """
    Calculate Bₙₘ(θ,ϕ), Cₙₘ(θ,ϕ), Pₙₘ(θ,ϕ) for all m and n
    """
    B_of_θ, C_of_θ = B_C_mn_of_θ_for_all_m_n(n_max, θ)
    P_of_θ = P_mn_of_θ_for_all_m_n(n_max, θ)

    for n = 1:n_max
        for m = -n:n
            factor = (-1)^m * convert(R, sqrt_factorial_n_plus_m_over_factorial_n_minus_m(m,n)) * exp(im * m * ϕ)
            B_of_θ[single_index_from_m_n(m, n)] *= factor
            C_of_θ[single_index_from_m_n(m, n)] *= factor
            P_of_θ[single_index_from_m_n(m, n)] *= factor
        end
    end

    return B_of_θ, C_of_θ, P_of_θ
end

function M_N_wave_all_m_n(n_max::I, kr::NN, θ::R, ϕ::R; kind="regular") where {R <: Real, I <: Integer, NN <: Number}
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
            M[single_index_from_m_n(m, n)] = convert(R, γ_mn(m, n)) * radial_function(n, kr) * C_of_θ_ϕ[single_index_from_m_n(m, n)]
            N[single_index_from_m_n(m, n)] = convert(R, γ_mn(m, n)) * (
                n * (n + 1) / kr * radial_function(n, kr)    * P_of_θ_ϕ[single_index_from_m_n(m, n)]
                + (radial_function_special_derivative(n, kr) * B_of_θ_ϕ[single_index_from_m_n(m, n)])
            )
        end
    end
    return M, N
end
