module CumulantsUpdates
  using Cumulants
  using SymmetricTensors
  using StatsBase
  using JLD
  import Cumulants: outerprodcum
  import SymmetricTensors: pyramidindices, getblockunsafe
  import Base: vecnorm

  include("updates.jl")
  include("operations.jl")

  export dataupdat, momentupdat, momentarray
  export cumulantscache, cumulantsupdat
  export vecnorm, moms2cums!, cums2moms, cumulantsupdate!, DataMoments, savedm, loaddm
end
