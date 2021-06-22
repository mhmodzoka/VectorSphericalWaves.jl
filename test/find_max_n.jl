using VectorSphericalWaves
using Memoize

th = (LinRange(0, pi, Int(10)));
"""
n,m = 7,-7
VectorSphericalWaves.wignerdjmn_recurrence.(n,0,m,th)
"""
s_max = 12
for f in [VectorSphericalWaves.wignerdjmn_ELZOUKA, VectorSphericalWaves.wignerdjmn_recurrence, VectorSphericalWaves.wignerdjmn_recurrence_memoize]
    println("==============================================================================")
    println("Now testing <<$f>>")
    for s = 1 : s_max
        for n = 1:s
            for m = -n:n
                try
                    wigner = f.(s,m,n,th);
                catch
                    println("using $f, here is the arguments that caused error: s=$s, n=$n, m=$m")
                    println()
                    break
                end

            end
        end
    end
end
