"""
    preallocate_buffer(A, nbuffers=Threads.nthreads())

Create a channel of preallocated buffer copies for thread-safe reuse.

# Arguments
- `A`: Template object to copy (typically an array). Each buffer slot is initialized
       as `deepcopy(A)` to ensure independence.
- `nbuffers`: Number of buffer copies (defaults to thread count)

# Example
buffer_pool = preallocate_buffer(zeros(100))
# Use with with_buffer() for thread-safe operations:
result = with_buffer(buffer_pool) do buf
    @. buf = sin(1:100)
    sum(buf)
end

"""
preallocate_buffer(A, nbuffers::Int=Threads.nthreads()) = preallocate_buffer(() -> deepcopy(A), nbuffers)

"""
    preallocate_buffer(init::Function, nbuffers=Threads.nthreads())

Create a channel of preallocated buffers using a custom initializer function.

Use this when buffers cannot/should not be deep-copied (e.g., FFT plans, objects
with external pointers) or when you want custom initialization.

# Arguments
- `init`: A zero-argument function that constructs and returns one fresh buffer.
- `nbuffers`: Number of buffer instances to create (defaults to thread count)

# Examples
# 1) Arrays with custom contents
pool = preallocate_buffer(() -> ones(256))

# 2) FFTW plans (do NOT deepcopy plans)
using FFTW
Ppool = preallocate_buffer(() -> plan_fft(zeros(ComplexF64, 1024)))

# Using with with_buffer()
s = with_buffer(Ppool) do P
    x = rand(ComplexF64, 1024)
    y = P * x
    sum(abs2, y)
end

"""
function preallocate_buffer(init::F, nbuffers::Int=Threads.nthreads()) where {F<:Function}
    @assert nbuffers >= 1
    first = init()
    ch = Channel{typeof(first)}(nbuffers)
    put!(ch, first)
    for _ in 2:nbuffers
        put!(ch, init())
    end
    return ch
end

"""
    with_buffer(f::Function, ch::Channel)

Borrow a buffer from the channel, execute function `f` with it, and return it.
Ensures buffer is returned even if an error occurs.

# Example
```julia
buffer_pool = preallocate_buffer(zeros(100))
result = with_buffer(buffer_pool) do buf
    buf .= calculate_something()
    sum(buf)
end
```
"""
@inline function with_buffer(f::F, ch::Channel{T}) where {T,F<:Function}
    buf = take!(ch)
    try
        return f(buf)
    finally
        put!(ch, buf)
    end
end

"""
    with_buffers(f::Function, channels::Channel...)

Borrow multiple buffers simultaneously and execute function `f` with them.
Ensures all buffers are returned even if an error occurs.

# Example
```julia
pool1 = preallocate_buffer(zeros(100))
pool2 = preallocate_buffer(zeros(50))
result = with_buffers(pool1, pool2) do buf1, buf2
    buf1 .= calculate_x()
    buf2 .= calculate_y()
    dot(buf1, buf2[1:100])
end
```
"""
# Base case: 1 channel
@inline with_buffers(f::Function, ch::Channel) = with_buffer(f, ch)

# Recursive case: 2+ channels
@inline function with_buffers(f::Function, ch1::Channel, chs::Channel...)
    with_buffer(ch1) do b1
        # capture b1 and keep borrowing the rest
        with_buffers((args...) -> f(b1, args...), chs...)
    end
end