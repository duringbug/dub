module dub
include("factorize/cholesky.jl")
include("factorize/QR.jl")
include("factorize/QR.jl")
function test()
    cholesky_result = cholesky([4 -1 -1; -1 17//4 17//4; 1 11//4 7//2])
    println(cholesky_result)
    QR([4 -1 -1; -1 17//4 17//4; 1 11//4 7//2], Givens)
end
end # module dub
dub.test()