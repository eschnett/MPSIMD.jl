using Base.Iterators
using Test
using Arbitrary

using MPSIMD



# Define type
const IntN128 = IntN{Float64, Int32, 4}
const maxval = big(2)^120       # don't go near overflow yet


# Generate arbitrary values
isinbounds(x) = abs(x) <= maxval
xs = collect(take(Base.Iterators.filter(isinbounds, arbitrary(BigInt)), 100))
ys = collect(take(Base.Iterators.filter(isinbounds, arbitrary(BigInt)), 100))
zs = collect(take(Base.Iterators.filter(isinbounds, arbitrary(BigInt)), 100))



@testset "Internal representation" begin
    z = IntN128()
    @test Tuple(z.digits) == (0, 0, 0, 0)
    i0 = IntN128(Int32(0))
    @test Tuple(i0.digits) == (0, 0, 0, 0)
    i1 = IntN128(Int32(1))
    @test Tuple(i1.digits) == (1, 0, 0, 0)
    @test Tuple(normalize(i1).digits) == (1, 0, 0, 0)
    im1 = IntN128(Int32(-1))
    @test Tuple(im1.digits) == (-1, 0, 0, 0)
    @test Tuple(normalize(im1).digits) == (-1, 0, 0, 0)
end



@testset "Conversions" begin
    for x in xs
        @test BigInt(IntN128(x)) == x
    end
end



@testset "Addition and subtraction" begin
    for x in xs
        @test BigInt(+IntN128(x)) == +x
        @test BigInt(-IntN128(x)) == -x
    end
    for (x,y) in zip(xs, ys)
        @test BigInt(IntN128(x) + IntN128(y)) == x + y
        @test BigInt(IntN128(x) - IntN128(y)) == x - y
    end
    s = IntN128(0)
    for x in xs
        s += IntN128(x)
    end
    @test BigInt(s) == sum(xs)
end



@testset "Scaling" begin
    for (x,y) in zip(xs, ys)
        y1 = y % Int32
        @test BigInt(IntN128(x) * y1) == x * y1
    end
end



# +z
# -z
# 
# i0 = zero(IntN128)
# i1 = IntN128(1)
# 
# @test i0 === normalize(i0)
# @test i1 === normalize(i1)
# 
# @show normalize(i0)
# @show normalize(i1)
# @show normalize(-i1)
# 
# @test i1 - i1 === i0
# @test i1 + -i1 === i0
# 
# using InteractiveUtils
# # @code_native MPSIMD.shift_up(i1.digits)
# # @code_native normalize1(i1)
# @code_native i0 + i0
