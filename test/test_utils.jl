using VectorSphericalWaves
n_max = 40
functions_tobe_tested = [VectorSphericalWaves.γ_mn, VectorSphericalWaves.γ_mn_dash]
for f in functions_tobe_tested
    println("now testing <<$f>>")
    for n = 1:n_max
        for m = -n:n
            try
                f(m, n);
            catch e
                println("here is the arguments that caused error: n=$n, m=$m")
                println(e)
                break
            end
        end
    end
end