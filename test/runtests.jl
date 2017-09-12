using Base.Test
using Cumulants
using SymmetricTensors
using Cumupdates
using JLD
import Cumupdates: rep, momentarray, moms2cums!, cums2moms, cnorms, normperelement
import Cumupdates: invers_gen, clcopulagen, cormatgen, tcopulagen, gcopulagen, u2stnormal, u2tdist
import Cumupdates: gcopulatmarg, tdistdat, tcopulagmarg, gcopulatmarg, normdist

srand(43)

te = [-0.112639 0.124715 0.124715 0.268717 0.124715 0.268717 0.268717 0.046154]
st = convert(SymmetricTensor, (reshape(te, (2,2,2))))

@testset "axiliary functions" begin
  @testset "rep" begin
    @test rep((1,2,3)) == 6
    @test rep((1,2,2)) == 3
    @test rep((1,1,1)) == 1
  end
  @testset "vecnorm" begin
    @test vecnorm(st) ≈ 0.5273572868359742
    @test vecnorm(st, 1) ≈ 1.339089
    @test vecnorm(st, 2.5) ≈ vecnorm(te, 2.5)
  end
end

X = randn(1000, 20)
Xup = rand(50, 20)
l = size(Xup, 1) + 1
Xprim = vcat(X, Xup)[l:end,:]

@testset "moment updates" begin
  @testset "simple test" begin
    x = ones(6, 2)
    y = 2*ones(2,2)
    m = moment(x, 3)
    @test convert(Array, moment(vcat(x,y)[3:end,:],3)) ≈ convert(Array, momentupdat(m, x, y))
  end
  m4 = moment(X, 4)
  m5 = moment(X, 5)
  @testset "block size = 2" begin
    mm4 = momentupdat(m4, X, Xup)
    mm5 = momentupdat(m5, X, Xup)
    @test convert(Array, mm4) ≈ convert(Array, moment(Xprim, 4))
    @test convert(Array, mm5) ≈ convert(Array, moment(Xprim, 5))
  end
  @testset "block size = 3" begin
    m = moment(X, 4, 3)
    aa = momentupdat(m, X, Xup)
    @test convert(Array, aa) ≈ convert(Array, moment(Xprim, 4))
  end
end
#=
@testset "moment exceptions" begin
  x = ones(10,4);
  y = 2*ones(5,3);
  m = moment(x, 3);
  @test_throws DimensionMismatch momentupdat(m, x, y)
  y = 2*ones(5,4)
  @test_throws DimensionMismatch momentupdat(m, x[:, 1:3], y)
  y = 2*ones(15,4)
  @test_throws BoundsError momentupdat(m, x, y)
end
=#
@testset "moments to cumulants" begin
  m1 = moment(X, 1)
  m2 = moment(X, 2)
  m3 = moment(X, 3)
  m4 = moment(X, 4)
  m5 = moment(X, 5)
  c = cumulants(X, 5)
  m2c = momentarray(X, 5)
  @testset "moment array" begin
    @test convert(Array, m2c[1]) ≈ convert(Array, m1)
    @test convert(Array, m2c[2]) ≈ convert(Array, m2)
    @test convert(Array, m2c[3]) ≈ convert(Array, m3)
    @test convert(Array, m2c[4]) ≈ convert(Array, m4)
    @test convert(Array, m2c[5]) ≈ convert(Array, m5)
  end
  @testset "moms2cums!" begin
    moms2cums!(m2c)
    @test convert(Array, c[1]) ≈ convert(Array, m2c[1])
    @test convert(Array, c[2]) ≈ convert(Array, m2c[2])
    @test convert(Array, c[3]) ≈ convert(Array, m2c[3])
    @test convert(Array, c[4]) ≈ convert(Array, m2c[4])
    @test convert(Array, c[5]) ≈ convert(Array, m2c[5])
  end
  @testset "cums2moms" begin
    mm = cums2moms(c);
    @test convert(Array, mm[1]) ≈ convert(Array, m1)
    @test convert(Array, mm[2]) ≈ convert(Array, m2)
    @test convert(Array, mm[3]) ≈ convert(Array, m3)
    @test convert(Array, mm[4]) ≈ convert(Array, m4)
    @test convert(Array, mm[5]) ≈ convert(Array, m5)
  end
end

@testset "cumulants update" begin
  @testset "simple test" begin
    x = ones(6, 2)
    y = 2*ones(2,2)
    c = cumulants(x, 3)
    cup = cumulantsupdat(c, x, y)
    cprim = cumulants(vcat(x,y)[3:end,:],3)
    @test convert(Array, cup[1]) ≈ convert(Array, cprim[1])
    @test convert(Array, cup[2]) ≈ convert(Array, cprim[2])
    @test convert(Array, cup[3]) ≈ convert(Array, cprim[3])
  end
  @testset "test on wider data" begin
    c = cumulants(X, 5)
    cc = cumulants(Xprim, 5)
    cup = cumulantsupdat(c, X, Xup)
    @test convert(Array, cc[1]) ≈ convert(Array, cup[1])
    @test convert(Array, cc[2]) ≈ convert(Array, cup[2])
    @test convert(Array, cc[3]) ≈ convert(Array, cup[3])
    @test convert(Array, cc[4]) ≈ convert(Array, cup[4])
    @test convert(Array, cc[5]) ≈ convert(Array, cup[5])
  end
end
#=
@testset "cumulants exceptions" begin
  x = ones(10,4);
  y = 2*ones(5,3);
  c = cumulants(x, 4);
  @test_throws DimensionMismatch cumulantsupdat(c, x, y)
  y = 2*ones(5,4)
  @test_throws DimensionMismatch cumulantsupdat(c, x[:, 1:3], y)
  @test_throws MethodError cumulantsupdat([c[1], c[3], c[4]], x, y)
  @test_throws MethodError cumulantsupdat([c[1], c[2], c[4], c[4]], x, y)
  @test_throws MethodError cumulantsupdat([c[1], c[2], c[4], c[3]], x, y)
  @test_throws UndefRefError cumulantsupdat([c[2], c[3], c[4]], x, y)
  @test_throws UndefRefError cumulantsupdat([c[2], c[4]], x, y)
  y = 2*ones(15,4)
  @test_throws BoundsError cumulantsupdat(c, x, y)
end
=#
@testset "generate data" begin
  @testset "axiliary functions" begin
    @test invers_gen([1., 2., 3., 4., 5.], 3.2) ≈ [0.638608, 0.535014, 0.478181, 0.44034, 0.412558] atol=1.0e-5
    srand(43)
    @test cormatgen(3) ≈ [1.0 -0.152492 0.629667; -0.152492 1.0 0.273588; 0.629667 0.273588 1.0] atol=1.0e-5
  end
  @testset "generate data from copuls" begin
    srand(43)
    @test clcopulagen(2,2) ≈ [0.982589 0.207961; 2.95138 2.77698] atol=1.0e-5
    srand(43)
    @test tcopulagen([[1. 0.5];[0.5 1.]], 2) ≈ [0.581625 0.792144; 0.76935 0.968669] atol=1.0e-5
    srand(43)
    @test gcopulagen([[1. 0.5];[0.5 1.]], 2) ≈ [0.589188 0.815308; 0.708285 0.924962] atol=1.0e-5
  end
  @testset "transform marginals" begin
    @test u2stnormal([0.2 0.4; 0.4 0.6; 0.6 0.8]) ≈ [-0.841621 -0.253347; -0.253347 0.253347; 0.253347 0.841621] atol=1.0e-5
    @test u2tdist([0.2 0.4; 0.4 0.6; 0.6 0.8], 10) ≈ [-0.879058  -0.260185; -0.260185 0.260185; 0.260185 0.879058] atol=1.0e-5
  end
  @testset "generates data given copula nad marginal" begin
    srand(43)
    @test gcopulatmarg([1. 0.5;0.5 1.], 2, 10) ≈ [0.231461 0.939967; 0.566673 1.55891] atol=1.0e-5
    srand(43)
    @test tcopulagmarg([1. 0.5;0.5 1.], 2, 10) ≈ [0.200129 0.784859; 0.87 2.0739] atol=1.0e-5
    srand(43)
    @test tdistdat([1. 0.5;0.5 1.], 2, 10) ≈ [0.205402 0.81778; 0.909864 2.388] atol=1.0e-4
    srand(43)
    @test normdist([1. 0.5;0.5 1.], 2, 10) ≈ [0.225457 0.897627; 0.548381 1.43926] atol=1.0e-5
  end
end

@testset "norms of arrys of cumulants" begin
  c = cumulants(X, 4)
  v = [vecnorm(c[i])/(c[i].dats^i) for i in 1:length(c)]
  v1 = [vecnorm(c[1]), vecnorm(c[2]), vecnorm(c[3])/(vecnorm(c[2])^(3/2)), vecnorm(c[4])/(vecnorm(c[2])^2)]
  @testset "vector of norms" begin
    @test cnorms(c, false) ≈ [vecnorm(cum) for cum in c]
    @test cnorms(c, false, 1.5) ≈ [vecnorm(cum, 1.5) for cum in c]
    @test normperelement(c) ≈ v
  end
  @testset "cumulants norms" begin
    @test cumnorms(X) ≈ v1
    @test load("/tmp/cumdata.jld", "x") ≈ X
    cload = load("/tmp/cumdata.jld", "c")
    @test convert(Array, cload[4]) ≈ convert(Array, c[4])
  end
  @testset "updated cumulants norms" begin
    Xprim = vcat(X, Xup)[(size(Xup, 1)+1):end,:]
    cup = cumulants(Xprim, 4)
    @test cumupdatnorms(Xup, false) ≈ [vecnorm(cum) for cum in cup]
    @test load("/tmp/cumdata.jld", "x") ≈ Xprim 
    cload = load("/tmp/cumdata.jld", "c")
    @test convert(Array, cload[1]) ≈ convert(Array, cup[1])
    @test convert(Array, cload[2]) ≈ convert(Array, cup[2])
    @test convert(Array, cload[3]) ≈ convert(Array, cup[3])
    @test convert(Array, cload[4]) ≈ convert(Array, cup[4])
  end
  @testset "load exception" begin
    save("/tmp/cumdata.jld", "t", 0.0)
    @test_throws ErrorException cumupdatnorms(Xup, false)
  end
end
