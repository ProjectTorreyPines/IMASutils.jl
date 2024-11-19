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

@testset "contour" begin
    xaxis, yaxis, vaxis = -0.25, 0.5, 1.0
    nx, ny = 10, 11
    x_coords = range(-2, 2, nx)
    y_coords = range(-1, 1, ny)
    values = [(x - xaxis)^2 + (y - yaxis)^2 for x in x_coords, y in y_coords] .+ vaxis

    x_cache, y_cache = contour_cache(values; aggression_level=1)
    @test length(x_cache)  == length(y_cache) == ceil(Int, π * sqrt(2.0 * (nx^2 + ny^2)))
    x_cache, y_cache = contour_cache(values; aggression_level=3)
    @test length(x_cache) == length(y_cache) == nx * ny
    x_cache, y_cache = contour_cache(values)
    @test length(x_cache) == length(y_cache) == 2.0 * (nx + ny)

    # test closed
    level = vaxis + 0.2
    x_contour, y_contour = contour_from_midplane!(x_cache, y_cache, values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test (x_contour[1] == x_contour[end]) && (y_contour[1] == y_contour[end])
    @test (length(x_contour) == length(y_contour) == 13)
    @test (x_contour isa SubArray) && (y_contour isa SubArray)

    xc2, yc2 = contour_from_midplane(values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test (xc2 isa Vector) && (yc2 isa Vector)
    @test all(x_contour .== xc2) && all(y_contour .== yc2)

    # test open
    level = vaxis + 1.0
    x_contour, y_contour = contour_from_midplane!(x_cache, y_cache, values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test (x_contour[1] != x_contour[end]) || (y_contour[1] != y_contour[end])
    @test (x_contour[1] in extrema(x_coords)) || (y_contour[1] in extrema(y_coords))
    @test (x_contour[end] in extrema(x_coords)) || (y_contour[end] in extrema(y_coords))
    @test (length(x_contour) == length(y_contour) == 22)

    # test empty
    level = vaxis + 6.0
    x_contour, y_contour = contour_from_midplane!(x_cache, y_cache, values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test isempty(x_contour) && isempty(y_contour)
    level = vaxis + 0.01
    x_contour, y_contour = contour_from_midplane!(x_cache, y_cache, values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test isempty(x_contour) && isempty(y_contour)


end