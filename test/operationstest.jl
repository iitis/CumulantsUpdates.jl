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
