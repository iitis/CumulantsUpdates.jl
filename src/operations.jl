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
vecnorm(bt::SymmetricTensor{N}, k=2)

Returns float, a k-norm of bt, for k != 0

```jldoctest
julia> vecnorm(st)
0.5273572868359742

julia> vecnorm(st, k=1)
1.339089
```
"""
function vecnorm{T <: AbstractFloat, N}(bt::SymmetricTensor{T, N}, k::Union{Float64, Int}=2)
  k != 0 || throw(AssertionError("0' th norm not supported"))
  ret = 0.
  for i in indices(N, bt.bln)
    @inbounds ret += rep(i) * vecnorm(bt[i...], k)^k
  end
  ret^(1/k)
end

"""
  cnorms{T <: AbstractFloat}(c::Vector{SymmetricTensor{T}}, norm::Bool, k::Union{Float64, Int})

Returns vector of Floats of k norms of cumulants of order 1, ..., m. If norm = true
.....

"""
function cnorms{T <: AbstractFloat}(c::Vector{SymmetricTensor{T}}, norm::Bool = true,
                                                                      k::Union{Float64, Int}=2)
  if norm
    mcov = [vecnorm(cum, k) for cum in c[1:2]]
    return [mcov..., [vecnorm(c[i], k)/(mcov[2])^(i/2) for i in 3:length(c)]...]
  else
    return [vecnorm(cum, k) for cum in c]
  end
end

"""
  normperelementT <: AbstractFloat}(c::Vector{SymmetricTensor{T}}, norm::Bool, k::Union{Float64, Int})

Returns vector of Floats of k norms of cumulants of order 1, ..., m, each value is divided by a number of elements in a given cumulant's tensor)

"""

normperelement{T <: AbstractFloat}(c::Vector{SymmetricTensor{T}}, norm::Bool = true, k::Union{Float64, Int}=2) =
 [vecnorm(c[i], k)/(c[i].dats)^i for i in 1:length(c)]

"""
  cumnorms{T <: AbstractFloat}(X::Matrix{T}, m::Int = 4, norm::Bool = true, k::Union{Float64, Int}=2)

Given multivariate data X computes a vector of cumulants of order 1,...,m cash them and data
in /tmp/cumdata.jld and returns a vector of Floats of k norm of those cumulants.
If norm = true each value is normalised (divided by a number of elements in a given cumulant's tensor)
"""

function cumnorms{T <: AbstractFloat}(X::Matrix{T}, m::Int = 4, norm::Bool = true, k::Union{Float64, Int}=2, b::Int = 3)
  c = cumulants(X, m, b)
  cpath = "/tmp/cumdata.jld"
  save(cpath, "c", c, "x", X)
  cnorms(c, norm, k)
end

"""
  cumupdatnorms{T <: AbstractFloat}(X::Matrix{T}, norm::Bool = true, k::Union{Float64, Int}=2)

Given multivariate data X and vector of cumulants of order 1,...,m stored in /tmp/cumdata.jld
updates data by Xup and cumulants, store them in "/tmp/cumdata.jld"
and returns a vector of Floats of k norm of updated cumulants.
If norm = true each value is normalised (divided by a number of elements in a given cumulant's tensor)
"""

function cumupdatnorms{T <: AbstractFloat}(Xup::Matrix{T}, norm::Bool = true, k::Union{Float64, Int}=2)
  cpath = "/tmp/cumdata.jld"
  isfile(cpath) || throw(AssertionError("no cumulants cashed please run cumnorms first"))
  c = load(cpath, "c")
  X = load(cpath, "x")
  Xprim = vcat(X, Xup)[(size(Xup, 1)+1):end,:]
  cup = cumulantsupdat(c, X, Xup)
  save(cpath, "c", cup, "x", Xprim)
  cnorms(cup, norm, k)
end
