using VectorSphericalWaves

"""
## Failure, doesn't support ^
using MultiFloat
MultiFloats.use_bigfloat_transcendentals()
"""
n_max, kr_r, kr_i, θ, ϕ = 2, 1e6, 1e2, 0.2, 0.3;
kr = complex(kr_r, kr_i);
M_N_wave_all_m_n(n_max, kr, θ, ϕ; kind="regular");


# BigFloat
n_max, kr_r, kr_i, θ, ϕ = 2, BigFloat(1e6), BigFloat(1e2), BigFloat(0.2), BigFloat(0.3);
kr = complex(kr_r, kr_i);
M_N_wave_all_m_n(n_max, kr, θ, ϕ; kind="regular");
