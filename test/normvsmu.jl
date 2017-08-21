#!/usr/bin/env julia

using Cumupdates
addprocs(6)
@everywhere using Cumupdates
using JLD
import Cumupdates: cormatgen, cumulants, tcopulagmarg, gcopulatmarg, tdistdat, normdist

function c4normvsmu(cormat::Matrix{Float64}, t::Int, wsize::Int, mu::Vector{Int})
  k = div(t, wsize)
  norms = zeros(k+1, length(mu))
  for j in 1:length(mu)
  x = normdist(cormat, t, mu[j])
  norms[1,:] .= cumnorms(x, 4, true, 2, 3)[4]
    for i in 1:k
      xup = tcopulagmarg(cormat, wsize, mu[j])
      x = vcat(x, xup)[(size(xup, 1)+1):end,:]
      norms[i+1,j] = cumupdatnorms(xup, true, 1)[4]
      println(wsize*i/t)
      println(mu[j])
    end
  end
  norms
end


function nvmu(t::Int=200000, tup::Int=50000, mu::Vector{Int} = [10, 12, 14], nar::Vector{Int} = [30])
  srand(43)
  for n in nar
    println(n)
    uv = [i*(tup/t) for i in 0:div(t, tup)]
    cormats = [cormatgen(n), 0.6*cormatgen(n)+0.4*eye(n), 0.4*cormatgen(n)+0.6*eye(n), cor(randn(2*n, n))]
    nor = zeros(length(cormats), length(mu), div(t, tup)+1)
    sl = zeros(length(cormats), length(mu))
    b  = zeros(length(cormats))
    for j in 1:length(cormats)
      println(vecnorm(cormats[j]))
      mod4 = c4normvsmu(cormats[j], t, tup, mu)
      slopes = zeros(Float64, mu)
      for i in 1:length(mu)
        slopes[i] = (linreg(uv[4:end], mod4[4:end,i])[2])^(-1)
        nor[j,i,:] = mod4[:,i]
      end
      sl[j,:] = slopes
      a,b[j] = linreg(mu, slopes)
      println(b[j])
    end
    str = "stats/nvmu"*string(n)*"_"*string(t)*"_"*string(tup)*".jld"
    filedict = Dict{String, Any}()
    push!(filedict, "mu" => mu)
    push!(filedict, "nc2" => [vecnorm(cormat) for cormat in cormats])
    push!(filedict, "n" => n)
    push!(filedict, "mVinv(sl)" => b)
    push!(filedict, "sl" => sl)
    push!(filedict, "norms" => nor)
    push!(filedict, "tuparr" => uv)
    save(str, filedict)
  end
end

function main()
  nvmu(2000000, 50000, [10, 12, 14, 16, 18, 20, 22, 24, 26], [25, 30, 35, 40])
end

main()
