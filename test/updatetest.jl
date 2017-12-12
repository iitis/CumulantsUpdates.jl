srand(43)

X = randn(100, 20)
Xup = rand(25, 20)
l = size(Xup, 1) + 1
Xprim = vcat(X, Xup)[l:end,:]
@testset "data updat" begin
  @test dataupdat(X, Xup) ≈ Xprim
end


@testset "moment updates" begin
  x = ones(6, 2)
  y = 2*ones(2,2)
  M3 = moment(x, 3)
  M4 = moment(x, 4)
  M3up = momentupdat(M3, x, y)
  @testset "simple test" begin
    Mup = moment(vcat(x,y)[3:end,:],3)
    @test convert(Array, Mup) ≈ convert(Array, M3up)
  end
  @testset "moment array" begin
    Ma = momentarray(x, 4)
    @test convert(Array, Ma[3]) ≈ convert(Array, M3)
    @test convert(Array, Ma[4]) ≈ convert(Array, M4)
    MM = momentupdat(Ma, x, y)
    @test convert(Array, M3up) ≈convert(Array, MM[3])
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
  m2c = [m1, m2, m3, m4, m5]
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

@testset "simple cumulants update" begin
  x = ones(6, 2)
  y = 2*ones(2,2)
  c1 = cumulantsupdat(x, 3, 2)
  cup1 = cumulantsupdat(x, y, 3)
  cprim = cumulants(vcat(x,y)[3:end,:],3)
  @test convert(Array, cup1[1]) ≈ convert(Array, cprim[1])
  @test convert(Array, cup1[2]) ≈ convert(Array, cprim[2])
  @test convert(Array, cup1[3]) ≈ convert(Array, cprim[3])
end


@testset "cumulants updates using cache" begin
    Xprim = vcat(X, Xup)[(size(Xup, 1)+1):end,:]
    cup = cumulants(Xprim, 4)
    c = cumulants(X)
    cc = cumulantsupdat(X)
    @test convert(Array, cc[3]) ≈ convert(Array, c[3])
    @test convert(Array, cc[4]) ≈ convert(Array, c[4])
    cc = cumulantsupdat(X, Xup, 4)
    @test convert(Array, cc[3]) ≈ convert(Array, cup[3])
    @test convert(Array, cc[4]) ≈ convert(Array, cup[4])
    Xprpr = dataupdat(Xprim, Xup)
    cc = cumulantsupdat(Xprim, Xup, 4)
    ccc = cumulants(Xprpr, 4)
    @test convert(Array, cc[3]) ≈ convert(Array, ccc[3])
    @test convert(Array, cc[4]) ≈ convert(Array, ccc[4])
    @test convert(Array, cumulantsupdat(X, Xup, 4)[4]) ≈ convert(Array, cup[4])
end


@testset "larger data cumulants update" begin
  cc = cumulants(Xprim, 5)
  c = cumulantsupdat(X, 5)
  cup = cumulantsupdat(X, Xup, 5)
  @testset "test on wider data" begin
    @test convert(Array, cc[1]) ≈ convert(Array, cup[1])
    @test convert(Array, cc[2]) ≈ convert(Array, cup[2])
    @test convert(Array, cc[3]) ≈ convert(Array, cup[3])
    @test convert(Array, cc[4]) ≈ convert(Array, cup[4])
    @test convert(Array, cc[5]) ≈ convert(Array, cup[5])
  end
  addprocs(2)
  eval(Expr(:toplevel, :(@everywhere using CumulantsUpdates)))
  @testset "multiprocessing cumulants update" begin
    c = cumulantsupdat(X, 4)
    cup = cumulantsupdat(X, Xup, 4)
    @test convert(Array, cc[1]) ≈ convert(Array, cup[1])
    @test convert(Array, cc[2]) ≈ convert(Array, cup[2])
    @test convert(Array, cc[3]) ≈ convert(Array, cup[3])
    @test convert(Array, cc[4]) ≈ convert(Array, cup[4])
  end
end

@testset "cumulants update exceptions" begin
  x = ones(10,4);
  y = 2*ones(5,3);
  c = cumulantsupdat(x, 4)
  @test_throws ArgumentError cumulantsupdat(x, y)
  y = 2*ones(5,4)
  @test_throws ArgumentError cumulantsupdat(x[:, 1:3], y)
  y = 2*ones(15,4)
  @test_throws BoundsError cumulantsupdat(x, y)
end
