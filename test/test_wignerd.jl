using VectorSphericalWaves
using WignerD

wignerd_functions_list = [
    WignerD.wignerdjmn,
    VectorSphericalWaves.wignerdjmn_ELZOUKA,
    VectorSphericalWaves.wignerdjmn_recurrence, 
    VectorSphericalWaves.wignerdjmn_recurrence_memoize,
    VectorSphericalWaves.∂wignerdjmn_by_∂θ
]

# checking the limit of s,m,n for wignerd
th = (LinRange(0, pi, Int(10)));

VectorSphericalWaves.∂wignerdjmn_by_∂θ(7,-7,7,th[2])

s_max = 15
for f in wignerd_functions_list
    println()
    println("==============================================================================")
    println("Now testing <<$f>>")
    @time for s = 1:s_max        
        for n = 1:s
            for m = -n:n
                try
                    wigner = f.(s, m, n, th);
                catch e
                    println("using $f, here is the arguments that caused error: s=$s, m=$m, n=$n")
                    bt = backtrace()
                    msg = sprint(showerror, e, bt)
                    println(msg)
                    #println(e)
                    println()
                    break
                end

            end
        end
    end
end

"""
sample results
==============================================================================
Now testing <<wignerdjmn_recurrence>>
  0.936349 seconds (6.94 M allocations: 212.850 MiB, 14.74% gc time)
==============================================================================
Now testing <<wignerdjmn_ELZOUKA>>
  1.108147 seconds (10.80 M allocations: 272.367 MiB, 30.52% gc time)
==============================================================================
Now testing <<wignerdjmn_recurrence_memoize>>
  0.050361 seconds (76.13 k allocations: 4.080 MiB)
"""
println()
println()
println("comparing wignerd_functions:")
println(wignerd_functions_list)
# TODO: WignerD.wignerdjmn is different from the functions we have in src/utils.jl near poles (θ=0 and θ=π)
# cross validation of Wignerd
for s = 1:s_max        
    for n = 1:s
        for m = -n:n
            try
                error_12 = filter(!isnan, abs.((wignerd_functions_list[1].(s, m, n, th) - wignerd_functions_list[2].(s, m, n, th)) ./ wignerd_functions_list[1].(s, m, n, th)));
                error_13 = filter(!isnan, abs.((wignerd_functions_list[1].(s, m, n, th) - wignerd_functions_list[3].(s, m, n, th)) ./ wignerd_functions_list[1].(s, m, n, th)));
                error_14 = filter(!isnan, abs.((wignerd_functions_list[1].(s, m, n, th) - wignerd_functions_list[4].(s, m, n, th)) ./ wignerd_functions_list[1].(s, m, n, th)));
                println("s=$s, n=$n, m=$m: error_12_max = $(max(error_12...)), error_13_max = $(max(error_13...)), error_14_max = $(max(error_14...))")
            catch
                println("using 00, here is the arguments that caused error: s=$s, n=$n, m=$m")
                println()
                break
            end

        end
    end
end
