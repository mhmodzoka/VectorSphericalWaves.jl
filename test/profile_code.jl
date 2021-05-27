using VectorSphericalWaves
using ProfileView
using Zygote

m, n, kr_r, kr_i, θ, ϕ = 1, 1, 1e5, 0, 0.2, 0.3
kr = complex(kr_r, kr_i)

θ = 0.4
@profiler for i=1:1e6; VectorSphericalWaves.πₘₙ(m, n, θ); end


@time for i=1:1e5; VectorSphericalWaves.πₘₙ(m, n, θ); end


Zygote.gradient(VectorSphericalWaves.πₘₙ, m, n, θ)

@time for i=1:1e2; N_mn_wave(m, n, complex(kr_r, kr_i), θ, ϕ; kind="regular"); end
