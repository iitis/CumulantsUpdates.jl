module Cumupdates
  using Cumulants
  using SymmetricTensors
  using Distributions
  using StatsBase
  using JLD
  import Cumulants: outerprodcum
  import SymmetricTensors: indices, getblockunsafe
  import Base: vecnorm

  include("updates.jl")

  include("gendata.jl")

  include("operations.jl")

  export momentupdat, cumulantsupdat, vecnorm, cumnorms, cumupdatnorms
end
