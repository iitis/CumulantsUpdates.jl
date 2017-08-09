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

## Installation

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
Input: a vector `[SymmetricTensor{T, 1}, SymmetricTensor{T, 2}, ...,SymmetricTensor{T, m}]` of cumulant tensors of order `1,2,...,m`, calculated for `X::Matrix{T}` such that `size(X) = (t, n)`, i.e. data with `n` marginal variables and `t` realisations; the update `Xup::Matrix{T}` such that `size(X) = (tup, n)` and `tup < t`. The input vector of cumulant tensors must by of sequent orders starting form `1`, otherwise error will be return. If cumulants in input vector are not computed for the same data, results are meaningless.

```julia
julia> x = ones(10,2);

julia> y = 2*ones(2,2);

julia> c = cumulants(x, 3);

julia> cumulantsupdat(c, x, y)
3-element Array{SymmetricTensors.SymmetricTensor{Float64,N},1}:
 SymmetricTensors.SymmetricTensor{Float64,1}(Nullable{Array{Float64,1}}[[1.2,1.2]],2,1,2,true)                                             
 SymmetricTensors.SymmetricTensor{Float64,2}(Nullable{Array{Float64,2}}[[0.16 0.16; 0.16 0.16]],2,1,2,true)                                
 SymmetricTensors.SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[0.096 0.096; 0.096 0.096] 
[0.096 0.096; 0.096 0.096]],2,1,2,true)
 
julia> cumulants(vcat(x,y)[3:end, :], 3)
3-element Array{SymmetricTensors.SymmetricTensor{Float64,N},1}:
 SymmetricTensors.SymmetricTensor{Float64,1}(Nullable{Array{Float64,1}}[[1.2,1.2]],2,1,2,true)                                             
 SymmetricTensors.SymmetricTensor{Float64,2}(Nullable{Array{Float64,2}}[[0.16 0.16; 0.16 0.16]],2,1,2,true)                                
 SymmetricTensors.SymmetricTensor{Float64,3}(Nullable{Array{Float64,3}}[[0.096 0.096; 0.096 0.096]
[0.096 0.096; 0.096 0.096]],2,1,2,true)

```
 
### Vector norm

```julia
julia> vecnorm{T <: AbstractFloat, m}(st::SymmetricTensor{T, m}; k::Union{Float64, Int})
```

Returns a vector norm of the `SymmetricTensors` type, `vecnorm(st, k) = vecnorn(convert(Array, st),k)`, for `k != 0`

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


# Performance analysis

To analyse the computational time of moments and cumulants updates vs a naive recalculation that uses `Cumulants.jl`, we supply the executable script `comptimes.jl`.
This script returns to a .jld file computational times, given following parameters:
* `-m (Int)`: cumulant's order, by default `m = 4`,
* `-n (vararg Int)`: numbers of marginal variables, by default `m = 30`,
* `-t (Int)`: number of realisations of random variable, by default `t = 2500000`,
* `-u (vararg Int)`: number of realisations of update, by default `u = 25000, 30000, 35000, 40000`,
* `-b (Int)`: blocks size, by default `b = 3`.

Computational times and parameters are saved in the .jld file in /res directory. All comparisons performed by this script use one core.

To analyse the computational time of cumulants updates for different block sizes `1 < b =< Int(sqrt(n))`, we supply the executable script `comptimesblocks.jl`.
This script returns to a .jld file computational times, given following parameters:
* `-m (Int)`: cumulant's order, by default `m = 4`,
* `-n (Int)`: numbers of marginal variables, by default `m = 48`,
* `-u (vararg Int)`: number of realisations of the update, by default `u = 20000, 40000`.

Computational times and parameters are saved in the .jld file in /res directory. All comparisons performed by this script use one core. 

To plot a graph run /res/plotcomptimes.jl followed by a `*.jld` file name


# Citing this work


Krzysztof Domino, Piotr Gawron, *On updates of high order cumulant tensors*, [arXiv:1701.06446](https://arxiv.org/abs/1701.06446)

