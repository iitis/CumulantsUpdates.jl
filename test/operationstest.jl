te = [-0.112639 0.124715 0.124715 0.268717 0.124715 0.268717 0.268717 0.046154]
st = convert(SymmetricTensornew, (reshape(te, (2,2,2))))
srand(43)
c = cumulants(randn(1000, 20), 4)

@testset "axiliary functions" begin
    @test rep((1,2,3)) == 6
    @test rep((1,2,2)) == 3
    @test rep((1,1,1)) == 1
end

@testset "vecnorm" begin
    @test vecnorm(st) ≈ 0.5273572868359742
    @test vecnorm(st, 1) ≈ 1.339089
    @test vecnorm(st, 2.5) ≈ vecnorm(te, 2.5)
    n = vecnorm(c[2])
    @test n ≈ vecnorm(convert(Array, c[2]))
    @test cnorms(c) ≈ [vecnorm(c[3])/(n^(3/2)), vecnorm(c[4])/(n^2)]
end
