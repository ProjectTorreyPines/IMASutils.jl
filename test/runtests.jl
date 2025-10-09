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

function is_ccw(x, y)
    N = length(x)
    a = 1
    b = N ÷ 3
    c = 2b
    A = (x[a], y[a])
    B = (x[b], y[b])
    C = (x[c], y[c])
    return (C[2] - A[2]) * (B[1] - A[1]) >= (B[2] - A[2]) * (C[1] - A[1])
end

@testset "contour" begin
    xaxis, yaxis, vaxis = -0.25, 0.5, 1.0
    nx, ny = 10, 11
    x_coords = range(-2, 2, nx)
    y_coords = range(-1, 1, ny)
    values = [(x - xaxis)^2 + (y - yaxis)^2 for x in x_coords, y in y_coords] .+ vaxis

    x_cache, y_cache = contour_cache(x_coords, y_coords; aggression_level=1)
    @test length(x_cache)  == length(y_cache) == ceil(Int, π * sqrt(2.0 * (nx^2 + ny^2)))
    x_cache, y_cache = contour_cache(x_coords, y_coords; aggression_level=3)
    @test length(x_cache) == length(y_cache) == nx * ny
    x_cache, y_cache = contour_cache(x_coords, y_coords)
    @test length(x_cache) == length(y_cache) == 2.0 * (nx + ny)

    # test closed
    level = vaxis + 0.2
    x_contour, y_contour = contour_from_midplane!(x_cache, y_cache, values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test (x_contour[1] == x_contour[end]) && (y_contour[1] == y_contour[end])
    @test (length(x_contour) == length(y_contour) == 13)
    @test (x_contour isa SubArray) && (y_contour isa SubArray)
    @test is_ccw(x_contour, y_contour)

    xc2, yc2 = contour_from_midplane(values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test (xc2 isa Vector) && (yc2 isa Vector)
    @test all(x_contour .== xc2) && all(y_contour .== yc2)

     # test dynamic resizing
     x_dynamic, y_dynamic = Vector{Float64}(undef, 1), Vector{Float64}(undef, 1)
     level = vaxis + 0.2
     xc2, yc2 = contour_from_midplane!(x_dynamic, y_dynamic, values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
     @test all(x_contour .== xc2) && all(y_contour .== yc2)
     @test (length(x_dynamic) == length(y_dynamic) == 2^4)  # should double four times

    # test reversed
    level = -(vaxis + 0.2)
    xc2, yc2 = contour_from_midplane!(x_cache, y_cache, .-values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test all(x_contour .== xc2) && all(y_contour .== yc2)
    @test is_ccw(x_contour, y_contour)

    # test open
    level = vaxis + 1.0
    x_contour, y_contour = contour_from_midplane!(x_cache, y_cache, values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test (x_contour[1] != x_contour[end]) || (y_contour[1] != y_contour[end])
    @test (x_contour[1] in extrema(x_coords)) || (y_contour[1] in extrema(y_coords))
    @test (x_contour[end] in extrema(x_coords)) || (y_contour[end] in extrema(y_coords))
    @test (length(x_contour) == length(y_contour) == 22)
    @test is_ccw(x_contour, y_contour)

    # test empty
    level = vaxis + 6.0
    x_contour, y_contour = contour_from_midplane!(x_cache, y_cache, values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test isempty(x_contour) && isempty(y_contour)
    level = vaxis + 0.01
    x_contour, y_contour = contour_from_midplane!(x_cache, y_cache, values, x_coords, y_coords, level, xaxis, yaxis, vaxis)
    @test isempty(x_contour) && isempty(y_contour)

end

@testset "buffer" begin
    # Test preallocate_buffer
    @testset "preallocate_buffer" begin
        nbuffers = 4
        template = zeros(100)
        pool = preallocate_buffer(template, nbuffers)

        @test pool isa Channel
        @test pool.sz_max == nbuffers

        # Should be able to take all buffers
        buffers = [take!(pool) for _ in 1:nbuffers]
        @test length(buffers) == nbuffers
        @test all(b -> b isa Vector{Float64} && length(b) == 100, buffers)

        # All buffers should be independent copies
        buffers[1][1] = 999.0
        @test buffers[2][1] == 0.0

        # Return buffers
        for b in buffers
            put!(pool, b)
        end

        # Test default nbuffers (should be Threads.nthreads())
        pool2 = preallocate_buffer(zeros(10))
        @test pool2.sz_max == Threads.nthreads()

        # Test function initializer signature
        pool3 = preallocate_buffer(() -> ones(50), 3)
        @test pool3 isa Channel
        @test pool3.sz_max == 3

        # Buffers should be initialized with custom values
        buffers3 = [take!(pool3) for _ in 1:3]
        @test all(b -> all(b .== 1.0), buffers3)

        # Each buffer should be independent
        buffers3[1][1] = 42.0
        @test buffers3[2][1] == 1.0

        # Return buffers
        for b in buffers3
            put!(pool3, b)
        end
    end

    # Test with_buffer
    @testset "with_buffer" begin
        pool = preallocate_buffer(zeros(50), 1)  # Use 1 buffer to test reuse

        # Basic usage
        result = with_buffer(pool) do buf
            buf .= 1:50
            sum(buf)
        end
        @test result == sum(1:50)

        # Buffer should be returned to pool
        @test isready(pool)

        # Test that buffer is actually reused (with 1 buffer pool)
        result2 = with_buffer(pool) do buf
            # Buffer should still contain previous data
            buf[1]
        end
        @test result2 == 1.0  # From previous use

        # Test error handling - buffer should still be returned
        @test_throws DivideError begin
            with_buffer(pool) do buf
                buf .= 42.0
                1 ÷ 0  # Force error
            end
        end

        # Pool should still have the buffer available
        @test isready(pool)

        # And the buffer should contain the data from failed operation
        result3 = with_buffer(pool) do buf
            buf[1]
        end
        @test result3 == 42.0

        # Test with multiple buffers in pool
        pool2 = preallocate_buffer(zeros(10), 3)
        results = Float64[]
        for i in 1:3
            push!(results, with_buffer(pool2) do buf
                buf .= i
                sum(buf)
            end)
        end
        @test results == [10.0, 20.0, 30.0]
    end

    # Test with_buffers (multiple buffers)
    @testset "with_buffers" begin
        pool1 = preallocate_buffer(zeros(10), 2)
        pool2 = preallocate_buffer(zeros(5), 2)
        pool3 = preallocate_buffer(zeros(3), 2)

        # Test with 1 buffer (base case)
        result = with_buffers(pool1) do buf
            buf .= 1.0
            sum(buf)
        end
        @test result == 10.0

        # Test with 2 buffers
        result = with_buffers(pool1, pool2) do buf1, buf2
            buf1 .= 1:10
            buf2 .= 1:5
            sum(buf1) + sum(buf2)
        end
        @test result == sum(1:10) + sum(1:5)

        # Test with 3 buffers
        result = with_buffers(pool1, pool2, pool3) do buf1, buf2, buf3
            buf1 .= 10.0
            buf2 .= 5.0
            buf3 .= 1.0
            sum(buf1) + sum(buf2) + sum(buf3)
        end
        @test result == 100.0 + 25.0 + 3.0

        # Test error handling with multiple buffers
        @test_throws DivideError begin
            with_buffers(pool1, pool2, pool3) do buf1, buf2, buf3
                buf1 .= 999.0
                buf2 .= 888.0
                buf3 .= 777.0
                1 ÷ 0
            end
        end

        # All buffers should be returned
        @test isready(pool1) && isready(pool2) && isready(pool3)
    end

    # Test thread safety
    @testset "thread safety" begin
        if Threads.nthreads() > 1
            nthreads = Threads.nthreads()
            pool = preallocate_buffer(zeros(100), nthreads)

            # Concurrent access from multiple threads
            results = Vector{Float64}(undef, nthreads)
            Threads.@threads for i in 1:nthreads
                results[i] = with_buffer(pool) do buf
                    buf .= i
                    sleep(0.001)  # Small delay to test concurrency
                    sum(buf)
                end
            end

            # Each thread should get its own buffer and compute correctly
            @test all(i -> results[i] == 100.0 * i, 1:nthreads)

            # Pool should have all buffers back
            @test pool.n_avail_items == nthreads
        else
            @info "Skipping multi-threaded tests (single thread detected)"
        end
    end

    # Test with different data types
    @testset "different types" begin
        # Test with Matrix
        pool_matrix = preallocate_buffer(zeros(Int, 5, 5), 2)
        result = with_buffer(pool_matrix) do buf
            buf .= 1
            sum(buf)
        end
        @test result == 25

        # Test with custom struct
        mutable struct TestStruct
            x::Float64
            y::Int
        end

        pool_struct = preallocate_buffer(TestStruct(0.0, 0), 2)
        result = with_buffer(pool_struct) do buf
            buf.x = 3.5
            buf.y = 7
            buf.x + buf.y
        end
        @test result == 10.5
    end

    # Test realistic use case
    @testset "realistic computation" begin
        # Simulate a common pattern: computing with temporary arrays
        pool = preallocate_buffer(zeros(1000), 4)

        function compute_with_temp(x, pool)
            with_buffer(pool) do temp
                temp .= sin.(x)
                temp .^= 2
                sum(temp)
            end
        end

        x = range(0, 2π, 1000)
        result1 = compute_with_temp(x, pool)
        result2 = compute_with_temp(x, pool)

        @test result1 ≈ result2 ≈ sum(sin.(x).^2)
    end
end