module dub
include("factorize/cholesky.jl")
include("factorize/QR.jl")
include("factorize/LU.jl")
include("test/test.jl")
test()
end
