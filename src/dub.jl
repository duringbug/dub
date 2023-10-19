module dub
include("factorize/cholesky.jl")

function test()
    cholesky_result = cholesky([4 -1 -1; -1 17//4 17//4; 1 11//4 7//2])
    println(cholesky_result)
end

end # module dub