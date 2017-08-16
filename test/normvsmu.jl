#!/usr/bin/env julia

function c4normvsmu(t::Int = 200000, n::Int = 20, wsize::Int = 50000, mu::Vector{Int} = [10, 12, 14, 16, 18, 20, 22])
  cormat = cormatgen(n);
  k = div(t, wsize)
  norms = zeros(k+1, length(mu))
  tup = [wsize*i/t for i in 0:k]
  for j in 1:length(mu)
  x = transpose(rand(MvNormal(cormat),t))
  norms[1,:] .= cumnorms(x, 4, true, 1)[4]
    for i in 1:k
      xup = gendata(cormat, wsize, mu[j])
      x = vcat(x, xup)[(size(xup, 1)+1):end,:]
      norms[i+1,j] = cumupdatnorms(xup, true, 1)[4]
      println(tup[i])
      println(mu[j])
    end
  end
  norms
end

c4 = c4normvsmu(200000, 24, 10000, [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30])
using PyPlot

mu = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
slopes = zeros(Float64, mu)
for i in 1:length(mu)
 slopes[i] = (linreg([i*0.1 for i in 0:10], c4[:,i])[2])^(-1)
end
plot(mu, slopes)
