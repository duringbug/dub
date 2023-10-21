struct QR{T<:Number}
    A::Matrix{T}
    function QR{T}(A) where {T<:Number}
        m, n = check_QR_matrix(A)
        new{T}(A)
    end
    function QR(A::Matrix{T}) where {T<:Number}
        m, n = check_QR_matrix(A)
        new{T}(A)
    end
    function check_QR_matrix(A::Matrix{T}) where {T<:Number}
        m, n = size(A)
        if m < n
            throw(ArgumentError("Matrix A must be m  n: $(size(A, 1)) < $(size(A, 2))"))
        end
        return m, n
    end
end