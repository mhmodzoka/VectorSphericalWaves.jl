using VectorSphericalWaves

m, n, kr_r, kr_i, θ, ϕ = 1, 1, 1e6, 1e2, 0.2, 0.3;
kr = complex(kr_r, kr_i);

M_mn_wave(m, n, kr, θ, ϕ; kind="regular");
N_mn_wave(m, n, kr, θ, ϕ; kind="regular");

# BigFloat
m, n, kr_r, kr_i, θ, ϕ = 1, 1, BigFloat(1e6), BigFloat(1e2), BigFloat(0.2), BigFloat(0.3);
kr = complex(kr_r, kr_i);

M_mn_wave(m, n, kr, θ, ϕ; kind="regular");
N_mn_wave(m, n, kr, θ, ϕ; kind="regular");


# SVector comparison
m, n, kr_r, kr_i, θ, ϕ = 1, 1, 1e6, 1e2, 0.2, 0.3;
kr = complex(kr_r, kr_i);

θ_array, ϕ_array = LinRange(0, π, 1000), LinRange(0, π, 1000)

println("M and N returning Arrays")
@time M_mn_wave.(m, n, kr, θ_array, ϕ_array; kind="regular");
@time N_mn_wave.(m, n, kr, θ_array, ϕ_array; kind="regular");

println("M and N returning SVector, you should notice the sligth speedup")
@time VectorSphericalWaves.M_mn_wave_SVector.(m, n, kr, θ_array, ϕ_array; kind="regular");
@time VectorSphericalWaves.N_mn_wave_SVector.(m, n, kr, θ_array, ϕ_array; kind="regular");

