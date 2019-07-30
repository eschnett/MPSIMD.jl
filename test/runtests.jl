using Base.Iterators
using Test
using Arbitrary

using MPSIMD



# Define type
const IntN128 = IntN{Float64, Int32, 4}
const maxval = big(2)^120       # don't go near overflow yet
const maxval2 = typeof(maxval)(sqrt(maxval))



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
    @test BigInt(zero(IntN128)) == 0
    @test BigInt(one(IntN128)) == 1
    for x in xs
        @test BigInt(IntN128(x)) == x
    end
end



@testset "Comparisons" begin
    for x in xs
        @test signbit(IntN128(x)) === signbit(x)
        @test BigInt(abs(IntN128(x))) == abs(x)
        @test cmp(IntN128(x), IntN128(x)) === 0
        @test (IntN128(x) == IntN128(x)) === true
        @test (IntN128(x) != IntN128(x)) === false
        @test (IntN128(x) < IntN128(x)) === false
        @test (IntN128(x) <= IntN128(x)) === true
        @test (IntN128(x) > IntN128(x)) === false
        @test (IntN128(x) >= IntN128(x)) === true
        @test isequal(IntN128(x), IntN128(x)) === true
        @test isless(IntN128(x), IntN128(x)) === false
    end
    for (x,y) in zip(xs, ys)
        @test cmp(IntN128(x), IntN128(y)) === cmp(x, y)
        @test (IntN128(x) == IntN128(y)) === (x == y)
        @test (IntN128(x) != IntN128(y)) === (x != y)
        @test (IntN128(x) < IntN128(y)) === (x < y)
        @test (IntN128(x) <= IntN128(y)) === (x <= y)
        @test (IntN128(x) > IntN128(y)) === (x > y)
        @test (IntN128(x) >= IntN128(y)) === (x >= y)
        @test isequal(IntN128(x), IntN128(y)) == isequal(x, y)
        @test isless(IntN128(x), IntN128(y)) == isless(x, y)
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



@testset "Multiplication" begin
    for (x,y) in zip(xs, ys)
        if abs(x * y) <= maxval
            @test BigInt(IntN128(x) * IntN128(y)) == x * y
        end
        x1 = x % maxval2
        y1 = y % maxval2
        @test BigInt(IntN128(x1) * IntN128(y1)) == x1 * y1
        if x != 0
            x2 = x
            y2 = y ÷ x
            @test BigInt(IntN128(x2) * IntN128(y2)) == x2 * y2
        end
        if x1 != 0
            x3 = x1
            y3 = y ÷ x1
            @test BigInt(IntN128(x3) * IntN128(y3)) == x3 * y3
        end
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



using BenchmarkTools



@benchmark x * y setup=((x,y)=rand(Int16,2)*big(2)^(128÷2-32))
@benchmark x * y setup=((x,y)=rand(Int16,2)*big(2)^(256÷2-32))
@benchmark x * y setup=((x,y)=rand(Int16,2)*big(2)^(512÷2-32))
@benchmark x * y setup=((x,y)=rand(Int16,2)*big(2)^(1024÷2-32))



using BitIntegers

@benchmark Int128(x) * Int128(y) setup=((x,y)=rand(Int,2))
@benchmark Int256(x) * Int256(y) setup=((x,y)=rand(Int,2))
@benchmark Int512(x) * Int512(y) setup=((x,y)=rand(Int,2))
@benchmark Int1024(x) * Int1024(y) setup=((x,y)=rand(Int,2))



const IntN256 = IntN{Float64, Int32, 8}
const IntN512 = IntN{Float64, Int32, 16}
const IntN1024 = IntN{Float64, Int32, 32}

@benchmark IntN128(x) * IntN128(y) setup=((x,y)=rand(Int,2))
# @benchmark IntN256(x) * IntN256(y) setup=((x,y)=rand(Int,2))
# @benchmark IntN512(x) * IntN512(y) setup=((x,y)=rand(Int,2))
# @benchmark IntN1024(x) * IntN1024(y) setup=((x,y)=rand(Int,2))
