using VectorSphericalWaves

wignerd_functions_list = [
    VectorSphericalWaves.wignerdjmn_recurrence, VectorSphericalWaves.wignerdjmn_ELZOUKA, 
    VectorSphericalWaves.wignerdjmn_recurrence_memoize
]

# checking the limit of s,m,n for wignerd
th = (LinRange(0, pi, Int(10)));
s_max = 20
for f in wignerd_functions_list
    println("==============================================================================")
    println("Now testing <<$f>>")
    @time for s = 1:s_max        
        for n = 1:s
            for m = -n:n
                try
                    wigner = f.(s, m, n, th);
                catch
                    println("using $f, here is the arguments that caused error: s=$s, n=$n, m=$m")
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

# cross validation of Wignerd
for s = 1:s_max        
    for n = 1:s
        for m = -n:n
            try
                error_12 = abs.((wignerd_functions_list[1].(s, m, n, th) - wignerd_functions_list[2].(s, m, n, th)) ./ wignerd_functions_list[1].(s, m, n, th));
                error_13 = abs.((wignerd_functions_list[1].(s, m, n, th) - wignerd_functions_list[3].(s, m, n, th)) ./ wignerd_functions_list[1].(s, m, n, th));
                println("s=$s, n=$n, m=$m: error_12 = $error_12, error_13 = $error_13, ")
            catch
                println("using 00, here is the arguments that caused error: s=$s, n=$n, m=$m")
                println()
                break
            end

        end
    end
end
