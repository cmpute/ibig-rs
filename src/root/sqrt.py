import math, random

def fixpoint(x, f):
    xn = f(x)
    while x < xn:
        x = xn
        xn = f(x)
    while x > xn:
        x = xn
        xn = f(x)
    return x

# newton's method on sqrt(n)
def sqrt_word(n):
    guess = 1 << (n.bit_length() // 2)
    # guess = int(math.sqrt(n))
    eq = lambda x: (n // x + x) // 2
    return fixpoint(guess, eq)

# rsqrt_tab[i] / 256 estimates (sqrt(32) / sqrt(32 + i)) (i.e. 8bits precision)
rsqrt_tab = [
    0xfe, 0xfa, 0xf7, 0xf3, 0xf0, 0xec, 0xe9, 0xe6, 0xe4, 0xe1, 0xde, 0xdc,
    0xd9, 0xd7, 0xd4, 0xd2, 0xd0, 0xce, 0xcc, 0xca, 0xc8, 0xc6, 0xc4, 0xc2,
    0xc1, 0xbf, 0xbd, 0xbc, 0xba, 0xb9, 0xb7, 0xb6, 0xb4, 0xb3, 0xb2, 0xb0,
    0xaf, 0xae, 0xac, 0xab, 0xaa, 0xa9, 0xa8, 0xa7, 0xa6, 0xa5, 0xa3, 0xa2,
    0xa1, 0xa0, 0x9f, 0x9e, 0x9e, 0x9d, 0x9c, 0x9b, 0x9a, 0x99, 0x98, 0x97,
    0x97, 0x96, 0x95, 0x94, 0x93, 0x93, 0x92, 0x91, 0x90, 0x90, 0x8f, 0x8e,
    0x8e, 0x8d, 0x8c, 0x8c, 0x8b, 0x8a, 0x8a, 0x89, 0x89, 0x88, 0x87, 0x87,
    0x86, 0x86, 0x85, 0x84, 0x84, 0x83, 0x83, 0x82, 0x82, 0x81, 0x81, 0x80,
]

# newton's method on 1/sqrt(n)
# x_{i+1} = x_i * (3 - n*x_i^2) / 2
# ref: https://stackoverflow.com/a/70428846/5960776
def sqrt_word_v2(n):
    # util: get the highest 32bits of a 64bit int
    hi32 = lambda a: a >> 32

    # step1: normalize input such that highest or second highest bit are set
    shift = (64 - n.bit_length()) & ((1 << 64) - 2)
    n <<= shift

    # step2: lookup initial estimation of 1/√x. The lookup table uses the highest 7 bits,
    # since the input is normalized, the lookup index must be larger than 2**(7-2) = 32
    n32 = hi32(n)
    r = rsqrt_tab[(n32 >> 25) - 32] # 8bit int
    print("r0:", r)
    
    # step3: first Newton iteration (without dividing by 2)
    # r should be an estimation of 1/√x with 16 bits precision
    r = ((3 * r) << 22) - hi32(n32 * ((r * r * r) << 8)) # 31 bits 3r - n*r^3 (highest bit will always be zero due to subtraction)
    print("r1:", r.bit_length(), r)

    # step4: second Newton iteration (without dividing by 2)
    # r should be an estimation of 1/√x with 32 bits precision now
    t = (3 << 28) - hi32(r * hi32(r * n32)) # 30bit 3 - r*r*n
    print("t:", t.bit_length(), t)
    r = hi32(r * t) # 28 bits (remove one additional bit, with the same reason as above)
    print("r2:", r.bit_length(), r)

    # step5: √x = x * 1/√x
    r <<= 4 # pad to 32 bits
    s = hi32(r * n32)
    s *= 2 # why?
    s -= 10 # make sure it's an underestimate
    print("s1:", s.bit_length(), s)

    # step6: third Newton iteration on √x
    e = n - s * s
    s = s + hi32((e >> 32) * r)
    print("s2:", s)

    # step6: fix the estimation error, at most 2 steps are needed
    # if we use more bits to estimate the initial guess, less steps can be required
    e = n - s * s
    if e >= s * 4 + 4:
        s += 2
    elif e >= s * 2 + 1:
        s += 1

    # step 7: restore the normalization
    return s >> (shift // 2)

# rsqrt_tab2[i] + 0x100 estimates (sqrt(32) / sqrt(32 + i)) (9 bits precision)
rsqrt_tab2 = [
    0xfc, 0xf4, 0xed, 0xe6, 0xdf, 0xd9, 0xd3, 0xcd, 0xc7, 0xc2, 0xbc, 0xb7,
    0xb2, 0xad, 0xa9, 0xa4, 0xa0, 0x9c, 0x98, 0x94, 0x90, 0x8c, 0x88, 0x85,
    0x81, 0x7e, 0x7b, 0x77, 0x74, 0x71, 0x6e, 0x6b, 0x69, 0x66, 0x63, 0x61,
    0x5e, 0x5b, 0x59, 0x57, 0x54, 0x52, 0x50, 0x4d, 0x4b, 0x49, 0x47, 0x45,
    0x43, 0x41, 0x3f, 0x3d, 0x3b, 0x39, 0x37, 0x36, 0x34, 0x32, 0x30, 0x2f,
    0x2d, 0x2c, 0x2a, 0x28, 0x27, 0x25, 0x24, 0x22, 0x21, 0x1f, 0x1e, 0x1d,
    0x1b, 0x1a, 0x19, 0x17, 0x16, 0x15, 0x14, 0x12, 0x11, 0x10, 0x0f, 0x0d,
    0x0c, 0x0b, 0x0a, 0x09, 0x08, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01,
]

def sqrt_word_v2b(n):
    hi32 = lambda a: a >> 32

    # step1: normalize input such that highest or second highest bit are set
    shift = (64 - n.bit_length()) & ((1 << 64) - 2)
    n <<= shift

    # step2: lookup initial estimation of 32/√x. The lookup table uses the highest 7 bits
    n32 = hi32(n)
    r = 0x100 | rsqrt_tab2[(n32 >> 25) - 32] # 9bit int
    print("r0:", r)
    
    # step3: first Newton iteration (without dividing by 2)
    # r should be a 31 bit estimation of 1/√x (*2) now
    r = ((3 * r) << 21) - hi32(n32 * ((r * r * r) << 5)) # 32bit 3r - n*r^3
    print("r1:", r.bit_length(), r)

    # step4: second Newton iteration (without dividing by 2)
    # r should be a 28 bit estimation of 1/√x (*16) now
    t = (3 << 28) - hi32(r * hi32(r * n32)) # 32bit 3 - r*r*n
    r = hi32(r * t)
    print("r2:", r.bit_length(), r)

    # step5: √x = x * 1/√x
    r <<= 4
    s = hi32(r * n32)
    s *= 2
    s -= 10 # make sure it's an underestimate
    print("s1:", s)

    # step6: third Newton iteration on √x
    e = n - s * s
    s = s + hi32((e >> 32) * r)
    print("s2:", s)

    # step6: fix the estimation error, at most 2 steps are needed
    # if we use more bits to estimate the initial guess, less steps can be required
    e = n - s * s
    print("e", e, "2s", s * 2)
    if e >= s * 4 + 4:
        s += 2
    elif e >= s * 2 + 1:
        s += 1

    # step 7: restore the normalization
    return s >> (shift // 2)

# binary iterative method
# XXX: this should be generalizable to cubic root?
# ref: https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Binary_numeral_system_(base_2)
# ref2: Henry S. Warren, Jr., "Hacker's Delight, 2nd e.d", p. 286
def sqrt_word_v3(n):
    x, c = n, 0

    # dₙ which starts at the highest power of four <= n
    d = 1 << 62 # The second-to-top bit is set.
                # Same as ((unsigned) INT32_MAX + 1) / 2.
    while d > n:
        d >>= 2

    # for dₙ … d₀
    while (d != 0):
        b = c | d
        c >>= 1
        if (x >= b):
            x -= b
            c = c | d
        d >>= 2
    return c

# newton's method on cubic root
def cbrt_word(n):
    guess = 1 << (n.bit_length() // 3)
    # guess = int(n ** 0.333)
    eq = lambda x: (n // (x**2) + x * 2) // 3
    return fixpoint(guess, eq)

# binary iterative method
def cbrt_word2(n):
    x = n
    y2 = 0
    y = 0
    smax = n.bit_length() // 3 + 1
    for s in reversed(range(smax)):
        s *= 3
        y2 *= 4
        y *= 2
        b = 3 * (y2 + y) + 1
        if (x >> s) >= b:
            x -= b << s
            y2 += 2 * y + 1
            y += 1

    return y

n = random.randrange(2, 2**64)
# n = 12345678901234567891
print(n, "GT sqrt", math.isqrt(n))
assert(sqrt_word(n) == math.isqrt(n))
assert(sqrt_word_v2(n) == math.isqrt(n))
assert(sqrt_word_v2b(n) == math.isqrt(n))
assert(sqrt_word_v3(n) == math.isqrt(n))

# n //= 2
r = cbrt_word(n)
assert(r**3 <= n and (r+1)**3 >= n)
assert(cbrt_word2(n) == r)

def generate(i, min_bits=5, rep_bits=8):
    '''
    Generate the estimation of 1/sqrt(i + 1/2**min_bits) * 2**rep_bits

    The result table can be used in the following way:
    the input number is first normalized so that the highest two bits are not both zero
    then select highest `(min_bits + 2)` bits to query the table, the result will be a fixed point
    estimation with `rep_bits` effective bits.
    '''
    x1 = 1.0 + i * 1.0 / (1 << min_bits)
    x2 = 1.0 + (i + 1) * 1.0 / (1 << min_bits)
    y = (1.0 / math.sqrt (x1) + 1.0 / math.sqrt (x2)) / 2
    val = int(y * (1 << rep_bits) + 0.5)
    return val

# for i in range(128 - 32):
#     print(hex(generate(i, rep_bits=9) - 256))

def generate_cbrt(i, min_bits=4, rep_bits=8):
    '''
    Generate the estimation of 1/cbrt(i + 1/2**min_bits) * 2**rep_bits

    The result table can be used in the following way:
    the input number is first normalized so that the highest three bits are not all zero
    then select highest `(min_bits + 3)` bits to query the table, the result will be a fixed point
    estimation with `rep_bits` effective bits.
    '''
    r3 = 1/3
    x1 = 1.0 + i * 1.0 / (1 << min_bits)
    x2 = 1.0 + (i + 1) * 1.0 / (1 << min_bits)
    y = (1.0 / x1**r3 + 1.0 / x2**r3) / 2
    val = int(y * (1 << rep_bits) + 0.5)
    return val

# for i in range(64 - 8):
#     print(hex(generate_cbrt(i, min_bits=3, rep_bits=9)))

# (rcbrt_tab[i] + 0x100) / 0x1ff ≈ (cbrt(8) / cbrt(8 + i))
rcbrt_tab = [
    0xf6, 0xe4, 0xd4, 0xc6, 0xb9, 0xae, 0xa4, 0x9b,
    0x92, 0x8a, 0x83, 0x7c, 0x76, 0x70, 0x6b, 0x66,
    0x61, 0x5c, 0x57, 0x53, 0x4f, 0x4b, 0x48, 0x44,
    0x41, 0x3e, 0x3b, 0x38, 0x35, 0x32, 0x2f, 0x2d,
    0x2a, 0x28, 0x25, 0x23, 0x21, 0x1f, 0x1d, 0x1b,
    0x19, 0x17, 0x15, 0x13, 0x11, 0x10, 0x0e, 0x0c,
    0x0b, 0x09, 0x08, 0x06, 0x05, 0x03, 0x02, 0x01
]


# iter for 1/cbrt(n): x(4 - ax^3)/3
def cbrt_word_v3(n):
    # use ∛2 to fix 1 bit error
    if n.bit_length() > 63:
        cbrt2 = 1385297844439 # 40bit estimation of ∛2, shoud be more than enough
        s, rem = cbrt_word_v3(n >> 1)
        s = s * cbrt2 >> 40

        e = n - s * s * s
        elim = 3 * s * (s + 1) + 1
        if e >= elim:
            print("outer fix")
            s += 1
            e -= elim

        return s, e

    print("cbrt", n)
    hi32 = lambda a: a >> 32

    # step1: normalize input such that at least one of the highest three bits are set
    shift = 63 - n.bit_length()
    shift -= shift % 3
    n <<= shift

    # step2: lookup initial estimation of 1/∛x. The lookup table uses the highest 6 bits
    # retrieved r ≈ ∛8 / ∛(n >> 57) * 0x200 = 1 / ∛(n >> 60) * 2^9 = 2^29 / ∛n.
    n32 = hi32(n)
    print("idx", (n32 >> 25) - 8)
    r = 0x100 | rcbrt_tab[(n32 >> 25) - 8] # 9bit int
    print("r0:", r)
    
    # step3: first Newton iteration
    # r = 2^52 / ∛n
    t = ((4 << 23) - hi32(n32 * (r * r * r)))
    r = r * t // 3 # 32bit
    print("r1:", r.bit_length(), r, 2**52 / n**(1/3))

    # step4: second Newton iteration
    # r = 2^48 / ∛n
    t = (4 << 28) - hi32(r * hi32(r * hi32(r * n32)))
    r = hi32(r * t) // 3 # 28bit
    print("r2:", r.bit_length(), r, 2**48 / n**(1/3))

    # step5: ∛x = x * (1/∛x)^2
    s = hi32(r * hi32(r * n32))
    print("s0:", s, n**(1/3), s < n**(1/3))

    # step6: fix the estimation error, at most 1 steps are needed
    # if we use more bits to estimate the initial guess, less steps can be required
    e = n - s * s * s
    elim = 3 * s * (s + 1) + 1
    while e >= elim:
        print("inner fix")
        s += 1
        e -= elim
        elim = 3 * s * (s + 1) + 1

    # step 7: restore the normalization
    return (s >> (shift // 3), e)

assert(cbrt_word_v3(12345678901234567891)[0] == 2311204)
assert(cbrt_word_v3(n)[0] == r, "error on n={}, cbrt={}".format(n, cbrt_word_v3(n)[0]))
print(cbrt_word_v3(6172839450617283945))
