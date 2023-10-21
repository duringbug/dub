"""
function QR(A, qr_type::QRType)

这个方法对A进行QR分解,返回矩阵Q与R

Arguments:
- `A`: 待分解的矩阵(M•N).
- `qr_type`: 分解的方法: @enum QRType Householder GramSchmidt Givens.
 
Returns:
QR分解的Q(M•M)与R(M•N)矩阵

Examples:
```
QR([4 -1 -1; -1 17//4 17//4; 1 11//4 7//2], Givens)
```
"""
struct QR{T<:Matrix{Number}}
end