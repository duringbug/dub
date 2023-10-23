"""
struct QR{T<:Number}

这个方法对A进行QR分解,返回矩阵Q与R

Arguments:
- `A`: 待分解的矩阵(m•n).
 
Returns:
- `A`: 分解后矩阵(n•n).经过了行变换
- `R`: R矩阵
- `Q`: Q矩阵

Examples:
```
qr = QR{Float64}([
1.0 1.0
2.0 0.0
2.0 1.0
])
println(qr.Q)
println(qr.R)
println(qr.Q * qr.R)
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
            w = x / sqrt(x' * x)
            W = zeros(m, m)
            W[i:3, i:3] = w * w'
            H = I - 2 * W
            A = H * A
            Q = inv(H) * Q
        end
        R = A
        return Q, R
    end
end