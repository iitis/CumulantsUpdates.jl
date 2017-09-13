using Base.Test
using Cumulants
using SymmetricTensors
using Cumupdates
using JLD
import Cumupdates: rep, momentarray, moms2cums!, cums2moms, cnorms, normperelement
import Cumupdates: invers_gen, clcopulagen, cormatgen, tcopulagen, gcopulagen, u2stnormal, u2tdist
import Cumupdates: gcopulatmarg, tdistdat, tcopulagmarg, gcopulatmarg, normdist

include("updatetest.jl")
include("gendattest.jl")
include("operationstest.jl")
