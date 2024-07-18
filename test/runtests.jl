using FuseUtils
using Test

@testset "FuseUtils.jl" begin
    N = 1000
    x = range(0, 1, N + 1)
    @test isapprox(trapz(x, x .^ 2), 1.0 / 3.0; rtol=1e-5)
    @test isapprox(cumtrapz(x, x .^ 2)[(N ÷ 2) + 1], 0.125 / 3.0; rtol=1e-5)
end
