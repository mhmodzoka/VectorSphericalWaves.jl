using VectorSphericalWaves

m, n, kr_r, kr_i, θ, ϕ = 1, 1, 1e5, 1e2, 0.2, 0.3
kr = complex(kr_r, kr_i)

M_mn_wave_SeparateRealImag_SMatrix(m, n, kr_r, kr_i, θ, ϕ, "regular");
N_mn_wave_SeparateRealImag_SMatrix(m, n, kr_r, kr_i, θ, ϕ; "regular");

# BigFloat
m, n, kr_r, kr_i, θ, ϕ = 1, 1, BigFloat(1e6), BigFloat(1e2), BigFloat(0.2), BigFloat(0.3);
kr = complex(kr_r, kr_i);

M_mn_wave(m, n, kr, θ, ϕ; kind="regular");
N_mn_wave(m, n, kr, θ, ϕ; kind="regular");