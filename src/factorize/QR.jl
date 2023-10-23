"""
struct QR{T<:Number}

这个方法对A进行QR分解,返回矩阵Q与R

Arguments:
- `A`: 待分解的矩阵(m•n).
 
Returns:
- `A`: 待分解矩阵(m•n)
- `R`: R矩阵(m•m)
- `Q`: Q矩阵(n•n)

Examples:
```
julia> using dub

julia> qr = dub.QR{Float64}([
               1.0 1.0
               2.0 0.0
               2.0 1.0
           ])
dub.QR{Float64}([1.0 1.0; 2.0 0.0; 2.0 1.0], [0.33333333333333354 0.6666666666666664 0.6666666666666665; 0.6666666666666665 -0.6666666666666667 0.33333333333333376; 0.6666666666666665 0.33333333333333365 -0.6666666666666666], [3.000000000000001 1.0; -5.551115123125783e-16 0.9999999999999998; -6.661338147750938e-16 2.220446049250313e-16])

julia> println(qr.Q)
[0.33333333333333354 0.6666666666666664 0.6666666666666665; 0.6666666666666665 -0.6666666666666667 0.33333333333333376; 0.6666666666666665 0.33333333333333365 -0.6666666666666666]

julia> println(qr.R)
[3.000000000000001 1.0; -5.551115123125783e-16 0.9999999999999998; -6.661338147750938e-16 2.220446049250313e-16]
```


"""
struct QR{T<:Number}
    A::Matrix{T}
    Q::Matrix{T}
    R::Matrix{T}
    function QR{T}(A) where {T<:Number}
        m, n = check_QR_matrix(A)
        Q, R = qr_factorize(A, m, n, T)
        new{T}(A, Q, R)
    end
    function QR(A::Matrix{T}) where {T<:Number}
        m, n = check_QR_matrix(A)
        Q, R = qr_factorize(A, m, n, T)
        new{T}(A, Q, R)
    end
    function check_QR_matrix(A::Matrix{T}) where {T<:Number}
        m, n = size(A)
        if m < n
            throw(ArgumentError("Matrix A must be m  n: $(size(A, 1)) < $(size(A, 2))"))
        end
        return m, n
    end
    function qr_factorize(A::Matrix{T}, m, n, S) where {T<:Number}
        Q = zeros(S, m, m)
        for i in 1:m
            Q[i, i] = one(S)
        end
        R = zeros(S, m, n)
        for i in 1:n
            I = zeros(S, m, m)
            for i in 1:m
                I[i, i] = one(S)
            end
            x = A[i:m, i]
            alpha = sqrt(x' * x)
            x[1] = x[1] - alpha
            w = x / alpha
            W = zeros(m, m)
            W[i:m, i:m] = w * w'
            H = I - 2 * W
            A = H * A
            Q = H * Q
        end
        R = A
        return inv(Q), R
    end
end