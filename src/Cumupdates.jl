module Cumupdates
  using Cumulants
  using SymmetricTensors
  using Distributions
  using StatsBase
  import Cumulants: outerprodcum
  import SymmetricTensors: indices
  import Base: vecnorm

  include("momentupdate.jl")

  include("gendata.jl")

  include("operations.jl")

  export momentupdat, momentarray, moms2cums!, cums2moms, cumulantsupdat, vecnorm
end
