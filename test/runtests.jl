using Base.Test
using Cumulants
using SymmetricTensors
using CumulantsUpdates
using JLD
import CumulantsUpdates: rep, momentarray, moms2cums!, cums2moms, cnorms, normperelement
import CumulantsUpdates: invers_gen, clcopulagen, cormatgen, tcopulagen, gcopulagen, u2stnormal, u2tdist
import CumulantsUpdates: gcopulatmarg, tdistdat, tcopulagmarg, gcopulatmarg, normdist

include("updatetest.jl")
include("gendattest.jl")
include("operationstest.jl")
