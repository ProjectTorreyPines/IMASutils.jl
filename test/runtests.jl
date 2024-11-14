using IMASutils
using Test

@testset "integration" begin
    N = 1000
    x = range(0, 1, N + 1)
    @test isapprox(trapz(x, x .^ 2), 1.0 / 3.0; rtol=1e-5)
    @test isapprox(cumtrapz(x, x .^ 2)[(N ÷ 2) + 1], 0.125 / 3.0; rtol=1e-5)

    g = rand(N + 1)
    @test trapz(x, g .* sin.(x)) ≈ trapz(x, (k, x) -> g[k] * sin(x))
    @test cumtrapz(x, g .* sin.(x))[(N ÷ 2) + 1] ≈ cumtrapz(x, (k, x) -> g[k] * sin(x))[(N ÷ 2) + 1]
end
