using Base.Test
using Cumulants
using SymmetricTensors
using CumulantsUpdates
using JLD
import CumulantsUpdates: rep, momentarray, moms2cums!, cums2moms, cnorms, normperelement

include("updatetest.jl")
include("operationstest.jl")
