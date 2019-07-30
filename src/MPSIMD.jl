module MPSIMD

using SIMD

struct Internal{F,I,N}
end

export IntN
struct IntN{F,I,N} #<: Signed
    # TODO: add a flag memoizing whether the representation is normalized

    digits::Vec{N,F}

    function IntN(::Internal{F1, I1, N1},
                  digits::Vec{N1, F1}) where {F1, I1, N1}
        @assert F1 <: AbstractFloat
        @assert I1 <: Signed
        @assert N1 isa Integer
        @assert N1 > 0
        new{F1, I1, N1}(digits)
    end
    function IntN(::Internal{F1, I1, N1}) where {F1, I1, N1}
        @assert F1 <: AbstractFloat
        @assert I1 <: Signed
        @assert N1 isa Integer
        @assert N1 > 0
        new{F1, I1, N1}(Vec{N1, F1}(0))
    end
end

function IntN{F,I,N}() where {F,I,N}
    IntN(Internal{F,I,N}())
end



function Base.show(io::IO, x::IntN{F,I,N}) where {F,I,N}
    print(io, "($F,$I,$N)[")
    join(io, [string(Integer(x.digits[n]), base=16) for n in N:-1:1], ", ")
    print(io, "]")
end



export normalize
function normalize(x::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    nbits = 8 * sizeof(I)
    carry = F(0)
    xs = x.digits
    rs = zero(typeof(xs))
    for n in 1:N
        rn = xs[n] + carry
        carry = round(rn * ldexp(F(1), -nbits))
        rn = rn - carry * ldexp(F(1), nbits)
        rs = setindex(rs, rn, n)
    end
    @assert carry == 0
    # @assert carry == 0 || carry == -1
    # rs = setindex(rs, rs[N] + carry * ldexp(F(1), nbits), N)
    IntN(Internal{F,I,N}(), rs)
end



function IntN{F,I,N}(x::IntN{F,I,N}) where {F,I,N}
    x
end

function IntN{F,I,N}(i::J) where {J<:Integer, F,I,N}
    nbits = 8 * sizeof(I)
    if Base.hastypemax(J) &&
        typemin(I) <= typemin(J) && typemax(J) <= typemax(I)
        ds = setindex(IntN{F,I,N}().digits, F(i), 1)
        return IntN(Internal{F,I,N}(), ds)
    end
    ds = IntN{F,I,N}().digits
    for n in 1:N
        i0 = i
        i = (i >> nbits) + (i < 0)
        d = F(i0 - (i << nbits))
        ds = setindex(ds, d, n)
    end
    @assert i == 0 || i == -1
    ds = setindex(ds, ds[N] + (i << nbits), N)
    IntN(Internal{F,I,N}(), ds)
end

function (::Type{J})(x::IntN{F,I,N}) where {J<:Integer, F,I,N}
    nbits = 8 * sizeof(I)
    r = J(0)
    for n in N:-1:1
        r <<= nbits
        r += J(x.digits[n])
    end
    r
end



function Base.cmp(x::IntN{F,I,N}, y::IntN{F,I,N})::Int where {F,I,N}
    J = SIMD.int_type(F)
    x1 = normalize(x)
    y1 = normalize(y)
    # r = cmp(x1.digits, y1.digits)
    r = vifelse(x1.digits == y1.digits, Vec{N, J}(0),
                vifelse(x1.digits < y1.digits, Vec{N, J}(-1), Vec{N, J}(1)))
    for n in N:-1:1
        r[n] != 0 && return r[n]
    end
    0
end

function Base. ==(x::IntN{F,I,N}, y::IntN{F,I,N})::Bool where {F,I,N}
    cmp(x, y) == 0
end

function Base.isless(x::IntN{F,I,N}, y::IntN{F,I,N})::Bool where {F,I,N}
    cmp(x, y) < 0
end

function Base.signbit(x::IntN{F,I,N})::Bool where {F,I,N}
    cmp(x, zero(x)) < 0
end

function Base.flipsign(x::IntN{F,I,N},
                       y::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    ifelse(signbit(y), -x, x)
end

function Base.abs(x::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    flipsign(x, x)
end

function Base.copysign(x::IntN{F,I,N},
                       y::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    flipsign(abs(x), y)
end



function Base.zero(::Type{IntN{F,I,N}})::IntN{F,I,N} where {F,I,N}
    IntN{F,I,N}()
end

function Base.one(::Type{IntN{F,I,N}})::IntN{F,I,N} where {F,I,N}
    IntN{F,I,N}(1)
end

function Base.zero(::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    zero(IntN{F,I,N})
end

function Base.one(::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    one(IntN{F,I,N})
end



function shift_up(x::Vec{N,F})::Vec{N,F} where {N,F}
    # @assert x[N] == 0
    shufflevector(x, zero(typeof(x)), Val{ntuple(n -> n==1 ? 2*N-1 : n-2, N)})
end

function shift_up(x::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    IntN(Internal{F,I,N}(), shift_up(x.digits))
end

export normalize1
function normalize1(::Type{I}, x::Vec{N,F}) where {I,N,F}
    nbits = 8 * sizeof(I)
    rhi = round(x * ldexp(F(1), -nbits))
    rlo = x - rhi * ldexp(F(1), nbits)
    rlo + shift_up(rhi)
end

function normalize1(x::IntN{F,I,N}) where {F,I,N}
    IntN(Internal{F,I,N}(), normalize1(I, x.digits))
end



function Base. +(x::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    x
end

function Base. -(x::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    IntN(Internal{F,I,N}(), -x.digits)
end



export fastadd
@inline function fastadd(x::IntN{F,I,N},
                  y::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    IntN(Internal{F,I,N}(), x.digits + y.digits)
end

export fastsub
@inline function fastsub(x::IntN{F,I,N},
                  y::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    IntN(Internal{F,I,N}(), x.digits - y.digits)
end

function Base. +(x::IntN{F,I,N},
                 y::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    normalize1(fastadd(x, y))
end

function Base. -(x::IntN{F,I,N},
                 y::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    normalize1(fastsub(x, y))
end



export fastmul
@inline function fastmul(x::IntN{F,I,N}, y::I)::IntN{F,I,N} where {F,I,N}
    fastmul(x, F(y))
end

@inline function fastmul(x::IntN{F,I,N}, y::F)::IntN{F,I,N} where {F,I,N}
    nbits = 8 * sizeof(I)
    nbits2 = nbits ÷ 2
    yhi = round(y * ldexp(F(1), -nbits2)) * ldexp(F(1), nbits2)
    ylo = y - yhi
    rlo = x.digits * ylo
    rhi = normalize1(I, x.digits * yhi)
    r = rlo + rhi
    IntN(Internal{F,I,N}(), r)
end

@inline function fastmul(x::IntN{F,I,N}, y::IntN{F,I,N}) where {F,I,N}
    r = zero(x)
    for n in N:-1:1
        r = fastadd(shift_up(r), fastmul(x, y.digits[n]))
    end
    r
end

function Base. *(x::IntN{F,I,N}, y::I)::IntN{F,I,N} where {F,I,N}
    normalize1(fastmul(x, y))
end

function Base. *(x::IntN{F,I,N}, y::IntN{F,I,N})::IntN{F,I,N} where {F,I,N}
    normalize1(fastmul(x, y))
end

end
