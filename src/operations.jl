"""
  rep(ind::Tuple)

axiliary for vecnorm count how many times a block should be counted
```jldoctest
julia> rep((1,2,3))
6

julia> rep((1,2,2))
3

julia> rep((1,1,1))
1
```
"""
function rep(ind::Tuple)
  @inbounds c = counts([ind...])
  div(factorial(length(ind)), mapreduce(factorial, * , c))
end

"""
vecnorm(bt::SymmetricTensor{N}, k=2)

Returns float, a k-norm of bt, for k != 0

```jldoctest
julia> vecnorm(st)
0.5273572868359742

julia> vecnorm(st, k=1)
1.339089
```
"""
function vecnorm(bt::SymmetricTensor{T, N}, k::Union{Float64, Int}=2) where T<:AbstractFloat where N
  k != 0 || throw(AssertionError("0' th norm not supported"))
  ret = 0.
  for i in indices(N, bt.bln)
    @inbounds ret += rep(i) * vecnorm(getblockunsafe(bt, i), k)^k
  end
  ret^(1/k)
end

"""
  cnorms(c::Vector{SymmetricTensor{T}})

Returns vector of Floats of norms of cumulants of order 3, ..., m, normalsed
by √||C₂||ᵏ
.....

"""
function cnorms(c::Vector{SymmetricTensor{T}}) where T <: AbstractFloat
  n = vecnorm(c[2])
  [vecnorm(c[k])/(n)^(k/2) for k in 3:length(c)]
end
