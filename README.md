# Cumupdates.jl

Updates moment tensors of any order for multivariate data aw well as a sequence of cumulant tensors of order `1,2,...,m`
Functions takes and return tensor or array of tensors in `SymmetricTensors` type. Requires [SymmetricTensors.jl](https://github.com/ZKSI/SymmetricTensors.jl). To convert to array, run:

```julia
julia> convert(Array, data::SymmetricTensors{T, N})
```


```julia
julia>  convert(SymmetricTensor, ca::Array{T,N})
```

As of 01/01/2017 [kdomino](https://github.com/kdomino) is the lead maintainer of this package.
