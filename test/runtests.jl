using FactCheck
using Cumulants
using SymmetricTensors
using Cumupdates
using JLD
import Cumupdates: rep, momentarray, moms2cums!, cums2moms, cnorms, normperelement

srand(43)

te = [-0.112639 0.124715 0.124715 0.268717 0.124715 0.268717 0.268717 0.046154]
st = convert(SymmetricTensor, (reshape(te, (2,2,2))))

facts("axiliary functions") do
  context("rep") do
    @fact rep((1,2,3)) --> 6
    @fact rep((1,2,2)) --> 3
    @fact rep((1,1,1)) --> 1
  end
  context("vecnorm") do
    @fact vecnorm(st) --> roughly(0.5273572868359742)
    @fact vecnorm(st, 1) --> roughly(1.339089)
    @fact vecnorm(st, 2.5) --> roughly(vecnorm(te, 2.5))
  end
end

X = randn(1000, 20)
Xup = rand(50, 20)
l = size(Xup, 1) + 1
Xprim = vcat(X, Xup)[l:end,:]

facts("moment updates") do
  context("simple test") do
    x = ones(6, 2)
    y = 2*ones(2,2)
    m = moment(x, 3)
    @fact convert(Array, moment(vcat(x,y)[3:end,:],3)) --> roughly(convert(Array, momentupdat(m, x, y)))
  end
  m4 = moment(X, 4)
  m5 = moment(X, 5)
  context("block size = 2") do
    mm4 = momentupdat(m4, X, Xup)
    mm5 = momentupdat(m5, X, Xup)
    @fact convert(Array, mm4) --> roughly(convert(Array, moment(Xprim, 4)))
    @fact convert(Array, mm5) --> roughly(convert(Array, moment(Xprim, 5)))
  end
  context("block size = 3") do
    m = moment(X, 4, 3)
    aa = momentupdat(m, X, Xup)
    @fact convert(Array, aa) --> roughly(convert(Array, moment(Xprim, 4)))
  end
end

facts("moment exceptions") do
  x = ones(10,4);
  y = 2*ones(5,3);
  m = moment(x, 3);
  @fact_throws DimensionMismatch momentupdat(m, x, y)
  y = 2*ones(5,4)
  @fact_throws DimensionMismatch momentupdat(m, x[:, 1:3], y)
  y = 2*ones(15,4)
  @fact_throws BoundsError momentupdat(m, x, y)
end

facts("moments to cumulants") do
  m1 = moment(X, 1)
  m2 = moment(X, 2)
  m3 = moment(X, 3)
  m4 = moment(X, 4)
  m5 = moment(X, 5)
  c = cumulants(X, 5)
  m2c = momentarray(X, 5)
  context("moment array") do
    @fact  convert(Array, m2c[1]) --> roughly(convert(Array, m1))
    @fact  convert(Array, m2c[2]) --> roughly(convert(Array, m2))
    @fact  convert(Array, m2c[3]) --> roughly(convert(Array, m3))
    @fact  convert(Array, m2c[4]) --> roughly(convert(Array, m4))
    @fact  convert(Array, m2c[5]) --> roughly(convert(Array, m5))
  end
  context("moms2cums!") do
    moms2cums!(m2c)
    @fact convert(Array, c[1]) --> roughly(convert(Array, m2c[1]))
    @fact convert(Array, c[2]) --> roughly(convert(Array, m2c[2]))
    @fact convert(Array, c[3]) --> roughly(convert(Array, m2c[3]))
    @fact convert(Array, c[4]) --> roughly(convert(Array, m2c[4]))
    @fact convert(Array, c[5]) --> roughly(convert(Array, m2c[5]))
  end
  context("cums2moms") do
    mm = cums2moms(c);
    @fact convert(Array, mm[1]) --> roughly(convert(Array, m1))
    @fact convert(Array, mm[2]) --> roughly(convert(Array, m2))
    @fact convert(Array, mm[3]) --> roughly(convert(Array, m3))
    @fact convert(Array, mm[4]) --> roughly(convert(Array, m4))
    @fact convert(Array, mm[5]) --> roughly(convert(Array, m5))
  end
end

facts("cumulants update") do
  context("simple test") do
    x = ones(6, 2)
    y = 2*ones(2,2)
    c = cumulants(x, 3)
    cup = cumulantsupdat(c, x, y)
    cprim = cumulants(vcat(x,y)[3:end,:],3)
    @fact convert(Array, cup[1]) --> roughly(convert(Array, cprim[1]))
    @fact convert(Array, cup[2]) --> roughly(convert(Array, cprim[2]))
    @fact convert(Array, cup[3]) --> roughly(convert(Array, cprim[3]))
  end
  c = cumulants(X, 5)
  cc = cumulants(Xprim, 5)
  cup = cumulantsupdat(c, X, Xup)
  @fact convert(Array, cc[1]) --> roughly(convert(Array, cup[1]))
  @fact convert(Array, cc[2]) --> roughly(convert(Array, cup[2]))
  @fact convert(Array, cc[3]) --> roughly(convert(Array, cup[3]))
  @fact convert(Array, cc[4]) --> roughly(convert(Array, cup[4]))
  @fact convert(Array, cc[5]) --> roughly(convert(Array, cup[5]))
end

facts("cumulants exceptions") do
  x = ones(10,4);
  y = 2*ones(5,3);
  c = cumulants(x, 4);
  @fact_throws DimensionMismatch cumulantsupdat(c, x, y)
  y = 2*ones(5,4)
  @fact_throws DimensionMismatch cumulantsupdat(c, x[:, 1:3], y)
  @fact_throws MethodError cumulantsupdat([c[1], c[3], c[4]], x, y)
  @fact_throws MethodError cumulantsupdat([c[1], c[2], c[4], c[4]], x, y)
  @fact_throws MethodError cumulantsupdat([c[1], c[2], c[4], c[3]], x, y)
  @fact_throws UndefRefError cumulantsupdat([c[2], c[3], c[4]], x, y)
  @fact_throws UndefRefError cumulantsupdat([c[2], c[4]], x, y)
  y = 2*ones(15,4)
  @fact_throws BoundsError cumulantsupdat(c, x, y)
end


facts("norms of arrys of cumulants") do
  c = cumulants(X, 4)
  v = [vecnorm(c[i])/(c[i].dats^i) for i in 1:length(c)]
  v1 = hcat([vecnorm(c[i]) for i in 1:2], [vecnorm(c[i])/(vecnorm(c[2])^(i/2)) for i in 3:length(c)])
  context("vector of norms") do
    @fact cnorms(c, false) --> roughly([vecnorm(cum) for cum in c])
    @fact cnorms(c, false, 1.5) --> roughly([vecnorm(cum, 1.5) for cum in c])
    @fact normperelement(c) --> roughly(v)
  end
  context("cumulants norms") do
    @fact cumnorms(X) --> roughly(v1)
    @fact load("/tmp/cumdata.jld", "x") --> roughly(X)
    cload = load("/tmp/cumdata.jld", "c")
    @fact convert(Array, cload[4]) -->  roughly(convert(Array, c[4]))
  end
  context("updated cumulants norms") do
    Xprim = vcat(X, Xup)[(size(Xup, 1)+1):end,:]
    cup = cumulants(Xprim, 4)
    @fact cumupdatnorms(Xup, false) --> roughly([vecnorm(cum) for cum in cup])
    @fact load("/tmp/cumdata.jld", "x") --> roughly(Xprim)
    cload = load("/tmp/cumdata.jld", "c")
    @fact convert(Array, cload[1]) -->  roughly(convert(Array, cup[1]))
    @fact convert(Array, cload[2]) -->  roughly(convert(Array, cup[2]))
    @fact convert(Array, cload[3]) -->  roughly(convert(Array, cup[3]))
    @fact convert(Array, cload[4]) -->  roughly(convert(Array, cup[4]))
  end
  context("load exception") do
    save("/tmp/cumdata.jld", "t", 0.0)
    @fact_throws ErrorException cumupdatnorms(Xup, false)
  end
end
