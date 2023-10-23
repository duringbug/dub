"""
struct Cholesky{T<:Number}

这个方法对A进行Cholesky分解,返回矩阵G

Arguments:
- `A`: 待分解的矩阵(n•n),正定且对称.
 
Returns:
- `A`: 待分解矩阵(n•n)
- `G`: A=G•G'

Examples:
```

```
"""
struct Cholesky{T<:Number}

    A::Matrix{T}
    n::Int
    R::Matrix{T}

    function Cholesky{T}(A) where {T<:Number}
        n, A = check_cholesky_matrix(A)
        R = cholesky_factorize(A, n, T)
        new{T}(A, n, R)
    end

    function Cholesky(A::Matrix{T}) where {T<:Number}
        n, A = check_cholesky_matrix(A)
        R = cholesky_factorize(A, n, T)
        new{T}(A, n, R)
    end

    function check_cholesky_matrix(A::Matrix{T}) where {T<:Number}
        if size(A, 1) != size(A, 2)
            throw(ArgumentError("Matrix A must be square: $(size(A, 1)) != $(size(A, 2))"))
        elseif A' != A
            throw(ArgumentError("Matrix A must be symmetric"))
        end
        n = size(A, 1)
        return n, A
    end

    function cholesky_factorize(A, n, S)
        R = zeros(S, n, n)
        v = zeros(n)
        for j in 1:n
            for i in j:n
                v[i] = A[i, j]
                for k in 1:j-1
                    v[i] = v[i] - R[j, k] * R[i, k]
                end
                if v[j] <= 0
                    throw(ArgumentError("Matrix A must be positive definite"))
                end
                R[i, j] = v[i] / sqrt(v[j])
            end
        end
        return R
    end
end