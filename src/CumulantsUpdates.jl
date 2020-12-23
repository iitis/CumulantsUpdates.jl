module CumulantsUpdates
  using Cumulants
  using SymmetricTensors
  using StatsBase
  using FileIO
  using JLD2
  using Distributed
  import Cumulants: outerprodcum
  import SymmetricTensors: pyramidindices, getblockunsafe
  import LinearAlgebra: norm
  if VERSION >= v"1.3"
    using CompilerSupportLibraries_jll
  end

  include("updates.jl")
  include("operations.jl")

  export dataupdat, momentupdat, momentarray
  export norm, moms2cums!, cums2moms, cumulantsupdate!, DataMoments, savedm, loaddm
end
