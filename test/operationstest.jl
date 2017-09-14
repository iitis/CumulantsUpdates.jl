te = [-0.112639 0.124715 0.124715 0.268717 0.124715 0.268717 0.268717 0.046154]
st = convert(SymmetricTensor, (reshape(te, (2,2,2))))
srand(43)

X = randn(1000, 20)
Xup = rand(50, 20)
l = size(Xup, 1) + 1
Xprim = vcat(X, Xup)[l:end,:]

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
  @testset "cumulants norms and updates" begin
    Xprim = vcat(X, Xup)[(size(Xup, 1)+1):end,:]
    cup = cumulants(Xprim, 4)
    @test cumnorms(X, 3, true, 2, 3, false) ≈ v1[1:3]
    nor, xx = cumupdatnorms(X, Xup, 4, false)
    @test nor ≈ cnorms(cup, false)
    @test xx == Xprim
    @test cumnorms(X) ≈ v1
    @test cumupdatnorms(X, Xup, 4, true)[1] ≈ cnorms(cup, true)
    @test cumupdatnorms(X, Xup, 4, true, 2.5)[1] ≈ cnorms(cup, true, 2.5)
    Xprpr = vcat(Xprim, Xup)[(size(Xup, 1)+1):end,:]
    @test cumupdatnorms(Xprim, Xup, 3, false, 2, 3, false)[1] ≈ cnorms(cumulants(Xprpr, 3), false)
    @test cumupdatnorms(Xprim, Xup, 3, true, 2.5)[1] ≈ cnorms(cumulants(Xprpr, 3), true, 2.5)
  end
end
