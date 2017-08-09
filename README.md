# Cumupdates.jl

Updates moment tensors of any order for multivariate data aw well as a sequence of cumulant tensors of order `1,2,...,m`
Functions takes and return tensor or array of tensors in `SymmetricTensors` type. Requires [SymmetricTensors.jl](https://github.com/ZKSI/SymmetricTensors.jl). To convert to array, run:

```julia
julia> convert(Array, st::SymmetricTensors{T, m})
```

```julia
julia>  convert(SymmetricTensor, a::Array{T,m})
```
Requires [Cumulants.jl](https://github.com/ZKSI/Cumulants.jl). Advised to compute first moment tensors or a vector of cumulants using `Cumulants.jl` that returns
`SymmetricTensors{T,m}` for N'th moment or `[SymmetricTensors{T,1}, SymmetricTensors{T,2}, ..., SymmetricTensors{T,m}]` for cumulant series of order 1,2,..., N

As of 01/01/2017 [kdomino](https://github.com/kdomino) is the lead maintainer of this package.

## Instalation

Within Julia, run

```julia
julia> Pkg.clone("https://bison.iitis.pl/kdomino/cumulant-update/tree/master/Cumupdates")
```

to install the files.  Julia 0.5 is required.

## Functions


### Moment update

```julia
julia> momentupdat{T<:AbstractFloat, m}(M::SymmetricTensor{T, m}, X::Matrix{T}, Xup::Matrix{T})
```

Returns a `SymmetricTensor{T, m}` of the moment tensor of order `m` of updated multivariate data in the observation window of size `t`: `X' = vcat(X,Xup)[1+size(Xup,1):end, :]` such that `size(X') = (t, n)`. Input: moment tensor `M` in the `SymmetricTensors` type with `m` modes and size `M.dats = n`, calculated for `X::Matrix{T}` such that `size(X) = (t, n)`, i.e. data with `n` marginal variables and `t` realisations; the update `Xup::Matrix{T}` such that `size(X) = (tup, n)` and `tup < t`

```julia
julia> x = ones(6, 2);

julia> m = moment(x, 3)
SymmetricTensors.SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[1.0 1.0; 1.0 1.0]

[1.0 1.0; 1.0 1.0]],2,1,2,true)

julia> y = 2*ones(2,2);

julia> momentupdat(m, x, y)
SymmetricTensors.SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[3.33333 3.33333; 3.33333 3.33333]

[3.33333 3.33333; 3.33333 3.33333]],2,1,2,true)

julia> moment(vcat(x,y)[3:end,:],3)
SymmetricTensors.SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[3.33333 3.33333; 3.33333 3.33333]

[3.33333 3.33333; 3.33333 3.33333]],2,1,2,true)

```

### Cumulants update


```julia
julia> cumulantsupdat{T<:AbstractFloat}(cum::Vector{SymmetricTensor{T}}, X::Matrix{T}, Xup::Matrix{T})
```
Returns a vector `[SymmetricTensor{T, 1}, SymmetricTensor{T, 2}, ...,SymmetricTensor{T, m}]` of cumulant tensors of order `1,2,...,m` of updated multivariate data in the observation window of size `t`: `X' = vcat(X,Xup)[1+size(Xup,1):end, :]` such that `size(X') = (t, n)`.
Input: a vector `[SymmetricTensor{T, 1}, SymmetricTensor{T, 2}, ...,SymmetricTensor{T, m}]` of cumulant tensors of order `1,2,...,m`, calculated for `X::Matrix{T}` such that `size(X) = (t, n)`, i.e. data with `n` marginal variables and `t` realisations; the update `Xup::Matrix{T}` such that `size(X) = (tup, n)` and `tup < t`.
 
 
### Vector norm

```julia
julia> vecnorm{T <: AbstractFloat, m}(st::SymmetricTensor{T, m}; k::Union{Float64, Int})
```

Returns a vector norm of the `SymmetricTensors` type, `vecnorm(st; k=k) = vecnorn(convert(Array, st),k)`, for `k \neq 0`

```julia
julia> te = [-0.112639 0.124715 0.124715 0.268717 0.124715 0.268717 0.268717 0.046154];

julia> st = convert(SymmetricTensor, (reshape(te, (2,2,2))));

julia> vecnorm(st)
0.5273572868359742

julia> vecnorm(st, 2.5)
0.4468668679541424

julia> vecnorm(st, 1)
1.339089

```
