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
  l = length(ind)
  c = counts([ind...])
  div(factorial(l), mapreduce(x -> factorial(x), * , c))
end

"""
vecnorm(bt::SymmetricTensor{N}; k=2)

Returns float, a k-norm of bt

```jldoctest
julia> vecnorm(st)
0.5273572868359742

julia> vecnorm(st ;k=1)
1.339089
```
"""
function vecnorm{T <: AbstractFloat, N}(bt::SymmetricTensor{T, N}; k=2)
  ret = 0.
  for i in indices(N, bt.bln)
    ret += rep(i) * vecnorm(bt[i...], k)^k
  end
  ret^(1/k)
end
