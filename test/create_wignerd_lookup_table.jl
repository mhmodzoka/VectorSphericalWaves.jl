import VectorSphericalWaves.wignerdjmn_recurrence_memoize
using HDF5
using JLD
using Interpolations

create_lookup_table = true
if create_lookup_table
    const wignerdjmn_VS_theta_cache = Dict()
    n_θ_points = Int(1e4);
    θ_array = LinRange(0, π, n_θ_points);
    for s = 1:6
        for m = -s:s
            for n = -s:s
                wignerdjmn_VS_theta_cache[s,m,n] = wignerdjmn_recurrence_memoize.(s, m, n, θ_array)
            end
        end
    end
    save("cach.jld", "wignerdjmn", wignerdjmn_VS_theta_cache, "θ_array", θ_array)
    
end
"""
    Calculate Wigner-d function, using recurrence and memoize, and saves the cache in a lookup table
"""
const wignerdjmn_VS_theta_cache = load("cach.jld")["wignerdjmn"]
const θ_array_interp = load("cach.jld")["θ_array"]
function wignerdjmn_recurrence_memoize_lookup(s::Int, m::Int, n::Int, θ::R) where R <: Real
    wignerdjmn_array = wignerdjmn_VS_theta_cache[s,m,n]        
    # itp = interpolate(θ_array_interp, wignerdjmn_array, BSpline(Linear()))
    itp = LinearInterpolation(θ_array_interp, wignerdjmn_array)
    return itp(θ)
end
s,m,n = 5, 2, 1
t_ = rand(10000);
@time VectorSphericalWaves.wignerdjmn_ELZOUKA.(s,m,n,t_);
@time VectorSphericalWaves.wignerdjmn_recurrence_memoize.(s,m,n,t_);
@time wignerdjmn_recurrence_memoize_lookup.(s,m,n,t_);

"""
Surprisingly, the direct calculation of WignerD is faster than memoizing!
"""