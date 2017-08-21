#!/usr/bin/env julia

using Cumupdates
addprocs(6)
@everywhere using Cumupdates
import Cumupdates: cormatgen, cumulants, tcopulagmarg, gcopulatmarg, tdistdat, normdist

function c4normvsmu(cormat::Matrix{Float64}, t::Int = 200000, wsize::Int = 50000, mu::Vector{Int} = [10, 12, 14, 16, 18, 20, 22])
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

using PyPlot

n = 30
t = 2000000
tup = 50000
mu = [10, 12, 14, 16, 18, 20, 22, 24, 26]
nu = [10,12,14]
cormats = [cormatgen(n), 0.5*cormatgen(n)+0.5*eye(n), 0.25*cormatgen(n)+0.75*eye(n), cor(randn(2*n, n))]

nor = zeros(length(cormats), length(mu), div(t, tup)+1)
sl = zeros(length(cormats), length(mu))
b  = zeros(length(cormats))
for j in 1:length(cormats)
  mod4 = c4normvsmu(cormats[j], t, tup, mu)
  slopes = zeros(Float64, mu)
  nor[j,:,:] = mod4
  for i in 1:length(mu)
    uv = [i*(tup/t) for i in 0:div(t, tup)]
   slopes[i] = (linreg(uv[3:end], mod4[3:end,i])[2])^(-1)
  end
  sl[j,:] = slopes
  #plot([i*(tup/t) for i in 0:div(t, tup)], mod4[1:end,1], "--s")
  a,b[j] = linreg(mu, slopes)
  println(b[j])
end

b
sl


plot(mu, sl[1,:], "d")
plot(mu, [a+b*m for m in mu], "--")
