#!/usr/bin/env julia

using Cumupdates
addprocs(6)
@everywhere using Cumupdates
import Cumupdates: cormatgen, cumulants, tcopulagmarg, gcopulatmarg, tdistdat, normdist

function c4normvsmu(t::Int = 200000, n::Int = 20, wsize::Int = 50000, mu::Vector{Int} = [10, 12, 14, 16, 18, 20, 22])
  cormat = cormatgen(n);
  k = div(t, wsize)
  norms = zeros(k+1, length(mu))
  tup = [wsize*i/t for i in 0:k]
  for j in 1:length(mu)
  x = normdist(cormat, t, mu[j])
  norms[1,:] .= cumnorms(x, 4, true, 2, 3)[4]
    for i in 1:k
      xup = tcopulagmarg(cormat, wsize, mu[j])
      x = vcat(x, xup)[(size(xup, 1)+1):end,:]
      norms[i+1,j] = cumupdatnorms(xup, true, 1)[4]
      println(tup[i])
      println(mu[j])
    end
  end
  norms
end

t = 500000
tup = 50000
mu = [10, 14, 18, 22, 26, 30]
n = 25
c4 = c4normvsmu(t, n, tup, mu)


slopes = zeros(Float64, mu)
for i in 1:length(mu)
 slopes[i] = (linreg([i*0.1 for i in 0:div(t, tup)], c4[:,i])[2])^(-1)
end

using PyPlot
plot(mu, slopes, "d")

linreg(mu, slopes)[2]
