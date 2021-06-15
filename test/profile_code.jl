using VectorSphericalWaves
using ProfileView
using Zygote

m, n, kr_r, kr_i, θ, ϕ = 6, 10, 1e5, 0, 0.2, 0.3
kr = complex(kr_r, kr_i)

θ = 0.4
@profiler for i=1:1e6; VectorSphericalWaves.πₘₙ(m, n, θ); end


@time for i=1:1e5; VectorSphericalWaves.πₘₙ(m, n, θ); end


Zygote.gradient(VectorSphericalWaves.πₘₙ, m, n, θ)

@time for i=1:1e2; N_mn_wave(m, n, complex(kr_r, kr_i), θ, ϕ; kind="regular"); end


using VectorSphericalWaves
m, n, kr_r, kr_i, θ, ϕ = 6, 10, 1e5, 0, 0.2, 0.3 ; kr = complex(kr_r, kr_i); θ = rand(1000);
N_mn_wave(m, n, kr, θ[1], ϕ; kind="regular")
M_mn_wave(m, n, kr, θ[1], ϕ; kind="regular")
@time for n = 1 : 6
    for m = -n:n
        N_mn_wave.(m, n, kr, θ, ϕ; kind="regular")
        M_mn_wave.(m, n, kr, θ, ϕ; kind="regular")
        N_mn_wave.(m, n, kr, θ, ϕ; kind="irregular")
        M_mn_wave.(m, n, kr, θ, ϕ; kind="irregular")
    end
end