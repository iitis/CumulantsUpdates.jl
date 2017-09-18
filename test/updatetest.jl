srand(43)

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
  c = cumulants(X, 5)
  cc = cumulants(Xprim, 5)
  cup = cumulantsupdat(c, X, Xup)
  @testset "simple test" begin
    x = ones(6, 2)
    y = 2*ones(2,2)
    c1 = cumulants(x, 3)
    cup1 = cumulantsupdat(c1, x, y)
    cprim = cumulants(vcat(x,y)[3:end,:],3)
    @test convert(Array, cup1[1]) ≈ convert(Array, cprim[1])
    @test convert(Array, cup1[2]) ≈ convert(Array, cprim[2])
    @test convert(Array, cup1[3]) ≈ convert(Array, cprim[3])
  end
  @testset "test on wider data" begin
    @test convert(Array, cc[1]) ≈ convert(Array, cup[1])
    @test convert(Array, cc[2]) ≈ convert(Array, cup[2])
    @test convert(Array, cc[3]) ≈ convert(Array, cup[3])
    @test convert(Array, cc[4]) ≈ convert(Array, cup[4])
    @test convert(Array, cc[5]) ≈ convert(Array, cup[5])
  end
  addprocs(2)
  eval(Expr(:toplevel, :(@everywhere using Cumupdates)))
  @testset "multiprocessing cumulants update" begin
    cupp = cumulantsupdat(c[1:4], X, Xup)
    @test convert(Array, cc[1]) ≈ convert(Array, cupp[1])
    @test convert(Array, cc[2]) ≈ convert(Array, cupp[2])
    @test convert(Array, cc[3]) ≈ convert(Array, cupp[3])
    @test convert(Array, cc[4]) ≈ convert(Array, cupp[4])
  end
end

@testset "cumulants exceptions" begin
  x = ones(10,4);
  y = 2*ones(5,3);
  c1 = cumulants(x, 4);
  @test_throws DimensionMismatch cumulantsupdat(c1, x, y)
  y = 2*ones(5,4)
  @test_throws DimensionMismatch cumulantsupdat(c1, x[:, 1:3], y)
  @test_throws MethodError cumulantsupdat([c1[1], c1[3], c1[4]], x, y)
  @test_throws MethodError cumulantsupdat([c1[1], c1[2], c1[4], c1[4]], x, y)
  @test_throws MethodError cumulantsupdat([c1[1], c1[2], c1[4], c1[3]], x, y)
  @test_throws MethodError cumulantsupdat([c1[2], c1[3], c1[4]], x, y)
  @test_throws MethodError cumulantsupdat([c1[2], c1[4]], x, y)
  y = 2*ones(15,4)
  @test_throws BoundsError cumulantsupdat(c1, x, y)
end
