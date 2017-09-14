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
function vecnorm(bt::SymmetricTensor{T, N}, k::Union{Float64, Int}=2) where T<:AbstractFloat where N
  k != 0 || throw(AssertionError("0' th norm not supported"))
  ret = 0.
  for i in indices(N, bt.bln)
    @inbounds ret += rep(i) * vecnorm(getblockunsafe(bt, i), k)^k
  end
  ret^(1/k)
end

"""
  cnorms(c::Vector{SymmetricTensor{T}}, norm::Bool, k::Union{Float64, Int})

Returns vector of Floats of k norms of cumulants of order 1, ..., m. If norm = true
.....

"""
function cnorms(c::Vector{SymmetricTensor{T}}, norm::Bool = true, k::Union{Float64, Int}=2) where T <: AbstractFloat
  if norm
    mcov = [vecnorm(cum, k) for cum in c[1:2]]
    return [mcov..., [vecnorm(c[i], k)/(mcov[2])^(i/2) for i in 3:length(c)]...]
  else
    return [vecnorm(cum, k) for cum in c]
  end
end

"""
  normperelement(c::Vector{SymmetricTensor{T}}, norm::Bool, k::Union{Float64, Int})

Returns vector of Floats of k norms of cumulants of order 1, ..., m, each value is divided by a number of elements in a given cumulant's tensor)

"""

normperelement(c::Vector{SymmetricTensor{T}}, norm::Bool = true, k::Union{Float64, Int}=2) where T <: AbstractFloat =
 [vecnorm(c[i], k)/(c[i].dats)^i for i in 1:length(c)]

"""
  cumnorms(X::Matrix{T}, m::Int = 4, norm::Bool = true, k::Union{Float64, Int}=2
                                                        cache::Bool = true, b::Int = 3)

Returns a vector of Floats of k norm of cumulants of order 1, ..., m calculated for data X.
If norm = true norms for m > 2 are normalised by the norm of the secod cumulant rised to m/2.
If cache save cumulants in /tmp/cumdata.jld.
b - the block size for the block structure of cumulants.
"""

function cumnorms(X::Matrix{T}, m::Int = 4, norm::Bool = true,
                                            k::Union{Float64, Int}=2, b::Int = 3,
                                            cache::Bool = true) where T <: AbstractFloat
  c = cumulants(X, m, b)
  if cache
    cpath = "/tmp/cumdata.jld"
    try save(cpath, "c", c, "x", X) catch end
  end
  cnorms(c, norm, k)
end

"""
  cumupdatnorms(X::Matrix{T}, Xup::Matrix{T}, m::Int = 4, norm::Bool = true,
                                                          k::Union{Float64, Int}=2,
                                                          b::Int = 3)

Returns a vector of Floats of k norm of cumulants of order 1, ..., m calculated for data X,
updated by Xup. Save updates cumulants in /tmp/cumdata.jld.
If cumulants are stored in "/tmp/cumdata.jld" the function uses it.
If norm = true norms for m > 2 are normalised by the norm of the secod cumulant rised to m/2.
b - the block size for the block structure of cumulants.
"""

function cumupdatnorms(X::Matrix{T}, Xup::Matrix{T}, m::Int = 4,
                                                     norm::Bool = true,
                                                     k::Union{Float64, Int}=2,
                                                     b::Int = 3,
                                                     cache::Bool = true) where T <: AbstractFloat
  Xprim = vcat(X, Xup)[(size(Xup, 1)+1):end,:]
  cpath = "/tmp/cumdata.jld"
  if isfile(cpath)
    d = load(cpath)
    X1 = try(d["x"]) catch end
    c = try(d["c"]) catch end
    if ((length(c) == m) & (X1 == X) & (typeof(c[1]) <: SymmetricTensor))
      cup = cumulantsupdat(c, X, Xup)
    else
      cup = cumulants(Xprim, m, b)
    end
  else
    cup = cumulants(Xprim, m, b)
  end
  if cache
    try save(cpath, "c", cup, "x", Xprim) catch end
  end
  cnorms(cup, norm, k), Xprim
end
