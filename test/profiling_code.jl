using VectorSphericalWaves
using ProfileView

m, n, kr_r, kr_i, θ, ϕ = 1, 1, 1e5, 0, 0.2, 0.3
kr = complex(kr_r, kr_i)

@profiler VectorSphericalWaves.πₘₙ( m, n, θ)
