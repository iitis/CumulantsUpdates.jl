module CumulantsUpdates
  using Cumulants
  using SymmetricTensors
  using StatsBase
  using JLD
  import Cumulants: outerprodcum
  import SymmetricTensors: indices, getblockunsafe
  import Base: vecnorm

  include("updates.jl")
  include("operations.jl")

  export momentupdat, momentarray
  export cumulantsupdat
  export vecnorm, moms2cums!, cums2moms
end