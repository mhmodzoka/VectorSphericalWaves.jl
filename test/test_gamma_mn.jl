using VectorSphericalWaves

n_max = 200

println("start testing `VectorSphericalWaves.γ_mn` ...")
for n = 1:n_max
    for m = -n:n
        try
            VectorSphericalWaves.γ_mn(m,n)
        catch e
            println("error for : n=$n, m=$m")
            """
            bt = backtrace()
            msg = sprint(showerror, e, bt)
            println(msg)
            """
            println(e)
            println()
            println()
        end
    end
end
println("end testing `VectorSphericalWaves.γ_mn` ...")
println()
println()
println("start testing `VectorSphericalWaves.γ_mn_dash` ...")
for n = 1:n_max
    for m = -n:n
        try
            VectorSphericalWaves.γ_mn_dash(m,n)
        catch e
            println("error for : n=$n, m=$m")
            """
            bt = backtrace()
            msg = sprint(showerror, e, bt)
            println(msg)
            """
            println(e)
            println()
            println()
        end
    end
end
println("end testing `VectorSphericalWaves.γ_mn_dash` ...")



