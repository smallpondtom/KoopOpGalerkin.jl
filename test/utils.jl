using Test
using KoopOpGalerkin

const KOG = KoopOpGalerkin

@testset "number of polynomial basis" begin
    @test KOG.numOfBasis(0,0) == 1
    @test KOG.numOfBasis(4,2) == 15
    @test KOG.numOfBasis(2,4) == 15
end