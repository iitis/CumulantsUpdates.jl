using FactCheck
using Cumulants
using SymmetricTensors
using Cumupdates
import Cumupdates: rep, vecnorm

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
    @fact vecnorm(st; k = 1) --> roughly(1.339089)
  end
end

X = randn(1000, 20)
Xup = rand(50, 20)
l = size(Xup, 1) + 1
Xprim = vcat(X, Xup)[l:end,:]

facts("moment updates") do
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
  c = cumulants(X, 5)
  cc = cumulants(Xprim, 5)
  cup = cumulantsupdat(c, X, Xup)
  @fact convert(Array, cc[1]) --> roughly(convert(Array, cup[1]))
  @fact convert(Array, cc[2]) --> roughly(convert(Array, cup[2]))
  @fact convert(Array, cc[3]) --> roughly(convert(Array, cup[3]))
  @fact convert(Array, cc[4]) --> roughly(convert(Array, cup[4]))
  @fact convert(Array, cc[5]) --> roughly(convert(Array, cup[5]))
end
