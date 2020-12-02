// void modmul768(const uint64_t a[12], const uint64_t b[12], const uint64_t m[13], uint64_t c[12])

// m[12] contains the least significant word of the negated inverse of the modulus mod 2^768

// Function parameters

#define a   x0
#define b   x1
#define m   x2
#define c   x3

// Internal register roles

// - multiplication

#define bi  x2

#define a0  x4
#define a1  x5
#define a2  x6
#define a3  x7
#define a4  x8
#define a5  x9
#define a6  x10
#define a7  x11
#define a8  x12
#define a9  x13
#define a10 x14
#define a11 x15

// - reduction

#define r   x1
#define inv x2

#define m0  a0
#define m1  a1
#define m2  a2
#define m3  a3
#define m4  a4
#define m5  a5
#define m6  a6
#define m7  a7
#define m8  a8
#define m9  a9
#define m10 a10
#define m11 a11

// - both

#define pp  x0

#define t0  x16
#define t1  x17
#define t2  x19
// skip PR (x18)
#define t3  x20
#define t4  x21
#define t5  x22
#define t6  x23
#define t7  x24
#define t8  x25
#define t9  x26
#define t10 x27
#define t11 x28
#define t12 t0
#define t13 t1
#define t14 t2
#define t15 t3
#define t16 t4
#define t17 t5
#define t18 t6
#define t19 t7
#define t20 t8
#define t21 t9
#define t22 t10
#define t23 t11
#define t24 x0

#ifdef __APPLE__
#define modmul768 _modmul768
#endif

.global modmul768

.text
.p2align 6,,63

modmul768:
	// Allocate space for storing 26 registers
	// - 10 for callee-saves
	// -  1 for m
	// -  1 for 16-byte alignment of SP
	// - 14 for intermediate values

	sub	SP, SP, #26*8

	// Load a

	ldp	 a0,  a1, [a, # 0*8]
	ldp	 a2,  a3, [a, # 2*8]
	ldp	 a4,  a5, [a, # 4*8]
	ldp	 a6,  a7, [a, # 6*8]
	ldp	 a8,  a9, [a, # 8*8]
	ldp	a10, a11, [a, #10*8]

	// Preserve m, r19-r28    

	stp	x19, x20, [SP, # 0*8]
	stp	x21, x22, [SP, # 2*8]
	stp	x23, x24, [SP, # 4*8]
	stp	x25, x26, [SP, # 6*8]
	stp	x27, x28, [SP, # 8*8]
	str	  m,      [SP, #10*8]

	// a * b0

	ldr	bi, [b, #0*8]

	mul	 t0,  a0, bi;	str	t0, [SP, #12*8]
	mul	 t1,  a1, bi
	mul	 t2,  a2, bi
	mul	 t3,  a3, bi
	mul	 t4,  a4, bi
	mul	 t5,  a5, bi
	mul	 t6,  a6, bi
	mul	 t7,  a7, bi
	mul	 t8,  a8, bi
	mul	 t9,  a9, bi
	mul	t10, a10, bi
	mul	t11, a11, bi

	umulh	pp,  a0, bi;	adds	 t1,  t1, pp
	umulh	pp,  a1, bi;	adcs	 t2,  t2, pp
	umulh	pp,  a2, bi;	adcs	 t3,  t3, pp
	umulh	pp,  a3, bi;	adcs	 t4,  t4, pp
	umulh	pp,  a4, bi;	adcs	 t5,  t5, pp
	umulh	pp,  a5, bi;	adcs	 t6,  t6, pp
	umulh	pp,  a6, bi;	adcs	 t7,  t7, pp
	umulh	pp,  a7, bi;	adcs	 t8,  t8, pp
	umulh	pp,  a8, bi;	adcs	 t9,  t9, pp
	umulh	pp,  a9, bi;	adcs	t10, t10, pp
	umulh	pp, a10, bi;	adcs	t11, t11, pp
	umulh	pp, a11, bi;	adc	t12, XZR, pp

	// a * b1

	ldr	bi, [b, #1*8]

	mul	pp,  a0, bi;	adds	 t1,  t1, pp;	str	t1, [SP, #13*8]
	mul	pp,  a1, bi;	adcs	 t2,  t2, pp
	mul	pp,  a2, bi;	adcs	 t3,  t3, pp
	mul	pp,  a3, bi;	adcs	 t4,  t4, pp
	mul	pp,  a4, bi;	adcs	 t5,  t5, pp
	mul	pp,  a5, bi;	adcs	 t6,  t6, pp
	mul	pp,  a6, bi;	adcs	 t7,  t7, pp
	mul	pp,  a7, bi;	adcs	 t8,  t8, pp
	mul	pp,  a8, bi;	adcs	 t9,  t9, pp
	mul	pp,  a9, bi;	adcs	t10, t10, pp
	mul	pp, a10, bi;	adcs	t11, t11, pp
	mul	pp, a11, bi;	adcs	t12, t12, pp
				adc	t13, XZR, XZR

	umulh	pp,  a0, bi;	adds	 t2,  t2, pp
	umulh	pp,  a1, bi;	adcs	 t3,  t3, pp
	umulh	pp,  a2, bi;	adcs	 t4,  t4, pp
	umulh	pp,  a3, bi;	adcs	 t5,  t5, pp
	umulh	pp,  a4, bi;	adcs	 t6,  t6, pp
	umulh	pp,  a5, bi;	adcs	 t7,  t7, pp
	umulh	pp,  a6, bi;	adcs	 t8,  t8, pp
	umulh	pp,  a7, bi;	adcs	 t9,  t9, pp
	umulh	pp,  a8, bi;	adcs	t10, t10, pp
	umulh	pp,  a9, bi;	adcs	t11, t11, pp
	umulh	pp, a10, bi;	adcs	t12, t12, pp
	umulh	pp, a11, bi;	adc	t13, t13, pp

	// a * b2

	ldr	bi, [b, #2*8]

	mul	pp,  a0, bi;	adds	 t2,  t2, pp;	str	t2, [SP, #14*8]
	mul	pp,  a1, bi;	adcs	 t3,  t3, pp
	mul	pp,  a2, bi;	adcs	 t4,  t4, pp
	mul	pp,  a3, bi;	adcs	 t5,  t5, pp
	mul	pp,  a4, bi;	adcs	 t6,  t6, pp
	mul	pp,  a5, bi;	adcs	 t7,  t7, pp
	mul	pp,  a6, bi;	adcs	 t8,  t8, pp
	mul	pp,  a7, bi;	adcs	 t9,  t9, pp
	mul	pp,  a8, bi;	adcs	t10, t10, pp
	mul	pp,  a9, bi;	adcs	t11, t11, pp
	mul	pp, a10, bi;	adcs	t12, t12, pp
	mul	pp, a11, bi;	adcs	t13, t13, pp
				adc	t14, XZR, XZR

	umulh	pp,  a0, bi;	adds	 t3,  t3, pp
	umulh	pp,  a1, bi;	adcs	 t4,  t4, pp
	umulh	pp,  a2, bi;	adcs	 t5,  t5, pp
	umulh	pp,  a3, bi;	adcs	 t6,  t6, pp
	umulh	pp,  a4, bi;	adcs	 t7,  t7, pp
	umulh	pp,  a5, bi;	adcs	 t8,  t8, pp
	umulh	pp,  a6, bi;	adcs	 t9,  t9, pp
	umulh	pp,  a7, bi;	adcs	t10, t10, pp
	umulh	pp,  a8, bi;	adcs	t11, t11, pp
	umulh	pp,  a9, bi;	adcs	t12, t12, pp
	umulh	pp, a10, bi;	adcs	t13, t13, pp
	umulh	pp, a11, bi;	adc	t14, t14, pp

	// a * b3

	ldr	bi, [b, #3*8]

	mul	pp,  a0, bi;	adds	 t3,  t3, pp;	str	t3, [SP, #15*8]
	mul	pp,  a1, bi;	adcs	 t4,  t4, pp
	mul	pp,  a2, bi;	adcs	 t5,  t5, pp
	mul	pp,  a3, bi;	adcs	 t6,  t6, pp
	mul	pp,  a4, bi;	adcs	 t7,  t7, pp
	mul	pp,  a5, bi;	adcs	 t8,  t8, pp
	mul	pp,  a6, bi;	adcs	 t9,  t9, pp
	mul	pp,  a7, bi;	adcs	t10, t10, pp
	mul	pp,  a8, bi;	adcs	t11, t11, pp
	mul	pp,  a9, bi;	adcs	t12, t12, pp
	mul	pp, a10, bi;	adcs	t13, t13, pp
	mul	pp, a11, bi;	adcs	t14, t14, pp
				adc	t15, XZR, XZR

	umulh	pp,  a0, bi;	adds	 t4,  t4, pp
	umulh	pp,  a1, bi;	adcs	 t5,  t5, pp
	umulh	pp,  a2, bi;	adcs	 t6,  t6, pp
	umulh	pp,  a3, bi;	adcs	 t7,  t7, pp
	umulh	pp,  a4, bi;	adcs	 t8,  t8, pp
	umulh	pp,  a5, bi;	adcs	 t9,  t9, pp
	umulh	pp,  a6, bi;	adcs	t10, t10, pp
	umulh	pp,  a7, bi;	adcs	t11, t11, pp
	umulh	pp,  a8, bi;	adcs	t12, t12, pp
	umulh	pp,  a9, bi;	adcs	t13, t13, pp
	umulh	pp, a10, bi;	adcs	t14, t14, pp
	umulh	pp, a11, bi;	adc	t15, t15, pp

	// a * b4

	ldr	bi, [b, #4*8]

	mul	pp,  a0, bi;	adds	 t4,  t4, pp;	str	t4, [SP, #16*8]
	mul	pp,  a1, bi;	adcs	 t5,  t5, pp
	mul	pp,  a2, bi;	adcs	 t6,  t6, pp
	mul	pp,  a3, bi;	adcs	 t7,  t7, pp
	mul	pp,  a4, bi;	adcs	 t8,  t8, pp
	mul	pp,  a5, bi;	adcs	 t9,  t9, pp
	mul	pp,  a6, bi;	adcs	t10, t10, pp
	mul	pp,  a7, bi;	adcs	t11, t11, pp
	mul	pp,  a8, bi;	adcs	t12, t12, pp
	mul	pp,  a9, bi;	adcs	t13, t13, pp
	mul	pp, a10, bi;	adcs	t14, t14, pp
	mul	pp, a11, bi;	adcs	t15, t15, pp
				adc	t16, XZR, XZR

	umulh	pp,  a0, bi;	adds	 t5,  t5, pp
	umulh	pp,  a1, bi;	adcs	 t6,  t6, pp
	umulh	pp,  a2, bi;	adcs	 t7,  t7, pp
	umulh	pp,  a3, bi;	adcs	 t8,  t8, pp
	umulh	pp,  a4, bi;	adcs	 t9,  t9, pp
	umulh	pp,  a5, bi;	adcs	t10, t10, pp
	umulh	pp,  a6, bi;	adcs	t11, t11, pp
	umulh	pp,  a7, bi;	adcs	t12, t12, pp
	umulh	pp,  a8, bi;	adcs	t13, t13, pp
	umulh	pp,  a9, bi;	adcs	t14, t14, pp
	umulh	pp, a10, bi;	adcs	t15, t15, pp
	umulh	pp, a11, bi;	adc	t16, t16, pp

	// a * b5

	ldr	bi, [b, #5*8]

	mul	pp,  a0, bi;	adds	 t5,  t5, pp;	str	t5, [SP, #17*8]
	mul	pp,  a1, bi;	adcs	 t6,  t6, pp
	mul	pp,  a2, bi;	adcs	 t7,  t7, pp
	mul	pp,  a3, bi;	adcs	 t8,  t8, pp
	mul	pp,  a4, bi;	adcs	 t9,  t9, pp
	mul	pp,  a5, bi;	adcs	t10, t10, pp
	mul	pp,  a6, bi;	adcs	t11, t11, pp
	mul	pp,  a7, bi;	adcs	t12, t12, pp
	mul	pp,  a8, bi;	adcs	t13, t13, pp
	mul	pp,  a9, bi;	adcs	t14, t14, pp
	mul	pp, a10, bi;	adcs	t15, t15, pp
	mul	pp, a11, bi;	adcs	t16, t16, pp
				adc	t17, XZR, XZR

	umulh	pp,  a0, bi;	adds	 t6,  t6, pp
	umulh	pp,  a1, bi;	adcs	 t7,  t7, pp
	umulh	pp,  a2, bi;	adcs	 t8,  t8, pp
	umulh	pp,  a3, bi;	adcs	 t9,  t9, pp
	umulh	pp,  a4, bi;	adcs	t10, t10, pp
	umulh	pp,  a5, bi;	adcs	t11, t11, pp
	umulh	pp,  a6, bi;	adcs	t12, t12, pp
	umulh	pp,  a7, bi;	adcs	t13, t13, pp
	umulh	pp,  a8, bi;	adcs	t14, t14, pp
	umulh	pp,  a9, bi;	adcs	t15, t15, pp
	umulh	pp, a10, bi;	adcs	t16, t16, pp
	umulh	pp, a11, bi;	adc	t17, t17, pp

	// a * b6

	ldr	bi, [b, #6*8]

	mul	pp,  a0, bi;	adds	 t6,  t6, pp;	str	t6, [SP, #18*8]
	mul	pp,  a1, bi;	adcs	 t7,  t7, pp
	mul	pp,  a2, bi;	adcs	 t8,  t8, pp
	mul	pp,  a3, bi;	adcs	 t9,  t9, pp
	mul	pp,  a4, bi;	adcs	t10, t10, pp
	mul	pp,  a5, bi;	adcs	t11, t11, pp
	mul	pp,  a6, bi;	adcs	t12, t12, pp
	mul	pp,  a7, bi;	adcs	t13, t13, pp
	mul	pp,  a8, bi;	adcs	t14, t14, pp
	mul	pp,  a9, bi;	adcs	t15, t15, pp
	mul	pp, a10, bi;	adcs	t16, t16, pp
	mul	pp, a11, bi;	adcs	t17, t17, pp
				adc	t18, XZR, XZR

	umulh	pp,  a0, bi;	adds	 t7,  t7, pp
	umulh	pp,  a1, bi;	adcs	 t8,  t8, pp
	umulh	pp,  a2, bi;	adcs	 t9,  t9, pp
	umulh	pp,  a3, bi;	adcs	t10, t10, pp
	umulh	pp,  a4, bi;	adcs	t11, t11, pp
	umulh	pp,  a5, bi;	adcs	t12, t12, pp
	umulh	pp,  a6, bi;	adcs	t13, t13, pp
	umulh	pp,  a7, bi;	adcs	t14, t14, pp
	umulh	pp,  a8, bi;	adcs	t15, t15, pp
	umulh	pp,  a9, bi;	adcs	t16, t16, pp
	umulh	pp, a10, bi;	adcs	t17, t17, pp
	umulh	pp, a11, bi;	adc	t18, t18, pp

	// a * b7

	ldr	bi, [b, #7*8]

	mul	pp,  a0, bi;	adds	 t7,  t7, pp;	str	t7, [SP, #19*8]
	mul	pp,  a1, bi;	adcs	 t8,  t8, pp
	mul	pp,  a2, bi;	adcs	 t9,  t9, pp
	mul	pp,  a3, bi;	adcs	t10, t10, pp
	mul	pp,  a4, bi;	adcs	t11, t11, pp
	mul	pp,  a5, bi;	adcs	t12, t12, pp
	mul	pp,  a6, bi;	adcs	t13, t13, pp
	mul	pp,  a7, bi;	adcs	t14, t14, pp
	mul	pp,  a8, bi;	adcs	t15, t15, pp
	mul	pp,  a9, bi;	adcs	t16, t16, pp
	mul	pp, a10, bi;	adcs	t17, t17, pp
	mul	pp, a11, bi;	adcs	t18, t18, pp
				adc	t19, XZR, XZR

	umulh	pp,  a0, bi;	adds	 t8,  t8, pp
	umulh	pp,  a1, bi;	adcs	 t9,  t9, pp
	umulh	pp,  a2, bi;	adcs	t10, t10, pp
	umulh	pp,  a3, bi;	adcs	t11, t11, pp
	umulh	pp,  a4, bi;	adcs	t12, t12, pp
	umulh	pp,  a5, bi;	adcs	t13, t13, pp
	umulh	pp,  a6, bi;	adcs	t14, t14, pp
	umulh	pp,  a7, bi;	adcs	t15, t15, pp
	umulh	pp,  a8, bi;	adcs	t16, t16, pp
	umulh	pp,  a9, bi;	adcs	t17, t17, pp
	umulh	pp, a10, bi;	adcs	t18, t18, pp
	umulh	pp, a11, bi;	adc	t19, t19, pp

	// a * b8

	ldr	bi, [b, #8*8]

	mul	pp,  a0, bi;	adds	 t8,  t8, pp;	str	t8, [SP, #20*8]
	mul	pp,  a1, bi;	adcs	 t9,  t9, pp
	mul	pp,  a2, bi;	adcs	t10, t10, pp
	mul	pp,  a3, bi;	adcs	t11, t11, pp
	mul	pp,  a4, bi;	adcs	t12, t12, pp
	mul	pp,  a5, bi;	adcs	t13, t13, pp
	mul	pp,  a6, bi;	adcs	t14, t14, pp
	mul	pp,  a7, bi;	adcs	t15, t15, pp
	mul	pp,  a8, bi;	adcs	t16, t16, pp
	mul	pp,  a9, bi;	adcs	t17, t17, pp
	mul	pp, a10, bi;	adcs	t18, t18, pp
	mul	pp, a11, bi;	adcs	t19, t19, pp
				adc	t20, XZR, XZR

	umulh	pp,  a0, bi;	adds	 t9,  t9, pp
	umulh	pp,  a1, bi;	adcs	t10, t10, pp
	umulh	pp,  a2, bi;	adcs	t11, t11, pp
	umulh	pp,  a3, bi;	adcs	t12, t12, pp
	umulh	pp,  a4, bi;	adcs	t13, t13, pp
	umulh	pp,  a5, bi;	adcs	t14, t14, pp
	umulh	pp,  a6, bi;	adcs	t15, t15, pp
	umulh	pp,  a7, bi;	adcs	t16, t16, pp
	umulh	pp,  a8, bi;	adcs	t17, t17, pp
	umulh	pp,  a9, bi;	adcs	t18, t18, pp
	umulh	pp, a10, bi;	adcs	t19, t19, pp
	umulh	pp, a11, bi;	adc	t20, t20, pp

	// a * b9

	ldr	bi, [b, #9*8]

	mul	pp,  a0, bi;	adds	 t9,  t9, pp;	str	t9, [SP, #21*8]
	mul	pp,  a1, bi;	adcs	t10, t10, pp
	mul	pp,  a2, bi;	adcs	t11, t11, pp
	mul	pp,  a3, bi;	adcs	t12, t12, pp
	mul	pp,  a4, bi;	adcs	t13, t13, pp
	mul	pp,  a5, bi;	adcs	t14, t14, pp
	mul	pp,  a6, bi;	adcs	t15, t15, pp
	mul	pp,  a7, bi;	adcs	t16, t16, pp
	mul	pp,  a8, bi;	adcs	t17, t17, pp
	mul	pp,  a9, bi;	adcs	t18, t18, pp
	mul	pp, a10, bi;	adcs	t19, t19, pp
	mul	pp, a11, bi;	adcs	t20, t20, pp
				adc	t21, XZR, XZR

	umulh	pp,  a0, bi;	adds	t10, t10, pp
	umulh	pp,  a1, bi;	adcs	t11, t11, pp
	umulh	pp,  a2, bi;	adcs	t12, t12, pp
	umulh	pp,  a3, bi;	adcs	t13, t13, pp
	umulh	pp,  a4, bi;	adcs	t14, t14, pp
	umulh	pp,  a5, bi;	adcs	t15, t15, pp
	umulh	pp,  a6, bi;	adcs	t16, t16, pp
	umulh	pp,  a7, bi;	adcs	t17, t17, pp
	umulh	pp,  a8, bi;	adcs	t18, t18, pp
	umulh	pp,  a9, bi;	adcs	t19, t19, pp
	umulh	pp, a10, bi;	adcs	t20, t20, pp
	umulh	pp, a11, bi;	adc	t21, t21, pp

	// a * b10

	ldr	bi, [b, #10*8]

	mul	pp,  a0, bi;	adds	t10, t10, pp;	str	t10, [SP, #22*8]
	mul	pp,  a1, bi;	adcs	t11, t11, pp
	mul	pp,  a2, bi;	adcs	t12, t12, pp
	mul	pp,  a3, bi;	adcs	t13, t13, pp
	mul	pp,  a4, bi;	adcs	t14, t14, pp
	mul	pp,  a5, bi;	adcs	t15, t15, pp
	mul	pp,  a6, bi;	adcs	t16, t16, pp
	mul	pp,  a7, bi;	adcs	t17, t17, pp
	mul	pp,  a8, bi;	adcs	t18, t18, pp
	mul	pp,  a9, bi;	adcs	t19, t19, pp
	mul	pp, a10, bi;	adcs	t20, t20, pp
	mul	pp, a11, bi;	adcs	t21, t21, pp
				adc	t22, XZR, XZR

	umulh	pp,  a0, bi;	adds	t11, t11, pp
	umulh	pp,  a1, bi;	adcs	t12, t12, pp
	umulh	pp,  a2, bi;	adcs	t13, t13, pp
	umulh	pp,  a3, bi;	adcs	t14, t14, pp
	umulh	pp,  a4, bi;	adcs	t15, t15, pp
	umulh	pp,  a5, bi;	adcs	t16, t16, pp
	umulh	pp,  a6, bi;	adcs	t17, t17, pp
	umulh	pp,  a7, bi;	adcs	t18, t18, pp
	umulh	pp,  a8, bi;	adcs	t19, t19, pp
	umulh	pp,  a9, bi;	adcs	t20, t20, pp
	umulh	pp, a10, bi;	adcs	t21, t21, pp
	umulh	pp, a11, bi;	adc	t22, t22, pp

	// a * b11, replace a with m

	ldr	bi, [b, #11*8];				ldr	b, [SP, #10*8]

	mul	pp,  a0, bi;	adds	t11, t11, pp;	str	t11, [SP, #23*8]
	mul	pp,  a1, bi;	adcs	t12, t12, pp
	mul	pp,  a2, bi;	adcs	t13, t13, pp
	mul	pp,  a3, bi;	adcs	t14, t14, pp
	mul	pp,  a4, bi;	adcs	t15, t15, pp
	mul	pp,  a5, bi;	adcs	t16, t16, pp
	mul	pp,  a6, bi;	adcs	t17, t17, pp
	mul	pp,  a7, bi;	adcs	t18, t18, pp
	mul	pp,  a8, bi;	adcs	t19, t19, pp
	mul	pp,  a9, bi;	adcs	t20, t20, pp
	mul	pp, a10, bi;	adcs	t21, t21, pp
	mul	pp, a11, bi;	adcs	t22, t22, pp
				adc	t23, XZR, XZR

	umulh	pp,  a0, bi;	adds	t12, t12, pp
	umulh	pp,  a1, bi;	adcs	t13, t13, pp;	ldp	 m0,  m1, [b, # 0*8]
	umulh	pp,  a2, bi;	adcs	t14, t14, pp
	umulh	pp,  a3, bi;	adcs	t15, t15, pp;	ldp	 m2,  m3, [b, # 2*8]
	umulh	pp,  a4, bi;	adcs	t16, t16, pp
	umulh	pp,  a5, bi;	adcs	t17, t17, pp;	ldp	 m4,  m5, [b, # 4*8]
	umulh	pp,  a6, bi;	adcs	t18, t18, pp
	umulh	pp,  a7, bi;	adcs	t19, t19, pp;	ldp	 m6,  m7, [b, # 6*8]
	umulh	pp,  a8, bi;	adcs	t20, t20, pp
	umulh	pp,  a9, bi;	adcs	t21, t21, pp;	ldp	 m8,  m9, [b, # 8*8]
	umulh	pp, a10, bi;	adcs	t22, t22, pp
	umulh	pp, a11, bi;	adc	t23, t23, pp;	ldp	m10, m11, [b, #10*8]
							ldr	inv, [b, #12*8]
	//// Multiplication done ////

	// Store t12-t23		// Load t0-t11

	stp	t12, t13, [SP, #24*8];	ldp	 t0,  t1, [SP, #12*8]
	stp	t14, t15, [SP, #12*8];	ldp	 t2,  t3, [SP, #14*8]
	stp	t16, t17, [SP, #14*8];	ldp	 t4,  t5, [SP, #16*8]
	stp	t18, t19, [SP, #16*8];	ldp	 t6,  t7, [SP, #18*8]
	stp	t20, t21, [SP, #18*8];	ldp	 t8,  t9, [SP, #20*8]
	stp	t22, t23, [SP, #20*8];	ldp	t10, t11, [SP, #22*8]

	// Reduction

	mul	r, inv, t0

	mul	pp,  m0, r;	adds	 t0,  t0, pp
	mul	pp,  m1, r;	adcs	 t1,  t1, pp
	mul	pp,  m2, r;	adcs	 t2,  t2, pp
	mul	pp,  m3, r;	adcs	 t3,  t3, pp
	mul	pp,  m4, r;	adcs	 t4,  t4, pp
	mul	pp,  m5, r;	adcs	 t5,  t5, pp
	mul	pp,  m6, r;	adcs	 t6,  t6, pp
	mul	pp,  m7, r;	adcs	 t7,  t7, pp
	mul	pp,  m8, r;	adcs	 t8,  t8, pp
	mul	pp,  m9, r;	adcs	 t9,  t9, pp
	mul	pp, m10, r;	adcs	t10, t10, pp
	mul	pp, m11, r;	adcs	t11, t11, pp
				adc	t12, XZR, XZR

	umulh	pp,  m0, r;	adds	 t1,  t1, pp
	umulh	pp,  m1, r;	adcs	 t2,  t2, pp
	umulh	pp,  m2, r;	adcs	 t3,  t3, pp
	umulh	pp,  m3, r;	adcs	 t4,  t4, pp
	umulh	pp,  m4, r;	adcs	 t5,  t5, pp
	umulh	pp,  m5, r;	adcs	 t6,  t6, pp
	umulh	pp,  m6, r;	adcs	 t7,  t7, pp
	umulh	pp,  m7, r;	adcs	 t8,  t8, pp
	umulh	pp,  m8, r;	adcs	 t9,  t9, pp
	umulh	pp,  m9, r;	adcs	t10, t10, pp
	umulh	pp, m10, r;	adcs	t11, t11, pp
	umulh	pp, m11, r;	adc	t12, t12, pp

	mul	r, inv, t1

	mul	pp,  m0, r;	adds	 t1,  t1, pp
	mul	pp,  m1, r;	adcs	 t2,  t2, pp
	mul	pp,  m2, r;	adcs	 t3,  t3, pp
	mul	pp,  m3, r;	adcs	 t4,  t4, pp
	mul	pp,  m4, r;	adcs	 t5,  t5, pp
	mul	pp,  m5, r;	adcs	 t6,  t6, pp
	mul	pp,  m6, r;	adcs	 t7,  t7, pp
	mul	pp,  m7, r;	adcs	 t8,  t8, pp
	mul	pp,  m8, r;	adcs	 t9,  t9, pp
	mul	pp,  m9, r;	adcs	t10, t10, pp
	mul	pp, m10, r;	adcs	t11, t11, pp
	mul	pp, m11, r;	adcs	t12, t12, pp
				adc	t13, XZR, XZR

	umulh	pp,  m0, r;	adds	 t2,  t2, pp
	umulh	pp,  m1, r;	adcs	 t3,  t3, pp
	umulh	pp,  m2, r;	adcs	 t4,  t4, pp
	umulh	pp,  m3, r;	adcs	 t5,  t5, pp
	umulh	pp,  m4, r;	adcs	 t6,  t6, pp
	umulh	pp,  m5, r;	adcs	 t7,  t7, pp
	umulh	pp,  m6, r;	adcs	 t8,  t8, pp
	umulh	pp,  m7, r;	adcs	 t9,  t9, pp
	umulh	pp,  m8, r;	adcs	t10, t10, pp
	umulh	pp,  m9, r;	adcs	t11, t11, pp
	umulh	pp, m10, r;	adcs	t12, t12, pp
	umulh	pp, m11, r;	adc	t13, t13, pp

	mul	r, inv, t2

	mul	pp,  m0, r;	adds	 t2,  t2, pp
	mul	pp,  m1, r;	adcs	 t3,  t3, pp
	mul	pp,  m2, r;	adcs	 t4,  t4, pp
	mul	pp,  m3, r;	adcs	 t5,  t5, pp
	mul	pp,  m4, r;	adcs	 t6,  t6, pp
	mul	pp,  m5, r;	adcs	 t7,  t7, pp
	mul	pp,  m6, r;	adcs	 t8,  t8, pp
	mul	pp,  m7, r;	adcs	 t9,  t9, pp
	mul	pp,  m8, r;	adcs	t10, t10, pp
	mul	pp,  m9, r;	adcs	t11, t11, pp
	mul	pp, m10, r;	adcs	t12, t12, pp
	mul	pp, m11, r;	adcs	t13, t13, pp
				adc	t14, XZR, XZR

	umulh	pp,  m0, r;	adds	 t3,  t3, pp
	umulh	pp,  m1, r;	adcs	 t4,  t4, pp
	umulh	pp,  m2, r;	adcs	 t5,  t5, pp
	umulh	pp,  m3, r;	adcs	 t6,  t6, pp
	umulh	pp,  m4, r;	adcs	 t7,  t7, pp
	umulh	pp,  m5, r;	adcs	 t8,  t8, pp
	umulh	pp,  m6, r;	adcs	 t9,  t9, pp
	umulh	pp,  m7, r;	adcs	t10, t10, pp
	umulh	pp,  m8, r;	adcs	t11, t11, pp
	umulh	pp,  m9, r;	adcs	t12, t12, pp
	umulh	pp, m10, r;	adcs	t13, t13, pp
	umulh	pp, m11, r;	adc	t14, t14, pp

	mul	r, inv, t3

	mul	pp,  m0, r;	adds	 t3,  t3, pp
	mul	pp,  m1, r;	adcs	 t4,  t4, pp
	mul	pp,  m2, r;	adcs	 t5,  t5, pp
	mul	pp,  m3, r;	adcs	 t6,  t6, pp
	mul	pp,  m4, r;	adcs	 t7,  t7, pp
	mul	pp,  m5, r;	adcs	 t8,  t8, pp
	mul	pp,  m6, r;	adcs	 t9,  t9, pp
	mul	pp,  m7, r;	adcs	t10, t10, pp
	mul	pp,  m8, r;	adcs	t11, t11, pp
	mul	pp,  m9, r;	adcs	t12, t12, pp
	mul	pp, m10, r;	adcs	t13, t13, pp
	mul	pp, m11, r;	adcs	t14, t14, pp
				adc	t15, XZR, XZR

	umulh	pp,  m0, r;	adds	 t4,  t4, pp
	umulh	pp,  m1, r;	adcs	 t5,  t5, pp
	umulh	pp,  m2, r;	adcs	 t6,  t6, pp
	umulh	pp,  m3, r;	adcs	 t7,  t7, pp
	umulh	pp,  m4, r;	adcs	 t8,  t8, pp
	umulh	pp,  m5, r;	adcs	 t9,  t9, pp
	umulh	pp,  m6, r;	adcs	t10, t10, pp
	umulh	pp,  m7, r;	adcs	t11, t11, pp
	umulh	pp,  m8, r;	adcs	t12, t12, pp
	umulh	pp,  m9, r;	adcs	t13, t13, pp
	umulh	pp, m10, r;	adcs	t14, t14, pp
	umulh	pp, m11, r;	adc	t15, t15, pp

	mul	r, inv, t4

	mul	pp,  m0, r;	adds	 t4,  t4, pp
	mul	pp,  m1, r;	adcs	 t5,  t5, pp
	mul	pp,  m2, r;	adcs	 t6,  t6, pp
	mul	pp,  m3, r;	adcs	 t7,  t7, pp
	mul	pp,  m4, r;	adcs	 t8,  t8, pp
	mul	pp,  m5, r;	adcs	 t9,  t9, pp
	mul	pp,  m6, r;	adcs	t10, t10, pp
	mul	pp,  m7, r;	adcs	t11, t11, pp
	mul	pp,  m8, r;	adcs	t12, t12, pp
	mul	pp,  m9, r;	adcs	t13, t13, pp
	mul	pp, m10, r;	adcs	t14, t14, pp
	mul	pp, m11, r;	adcs	t15, t15, pp
				adc	t16, XZR, XZR

	umulh	pp,  m0, r;	adds	 t5,  t5, pp
	umulh	pp,  m1, r;	adcs	 t6,  t6, pp
	umulh	pp,  m2, r;	adcs	 t7,  t7, pp
	umulh	pp,  m3, r;	adcs	 t8,  t8, pp
	umulh	pp,  m4, r;	adcs	 t9,  t9, pp
	umulh	pp,  m5, r;	adcs	t10, t10, pp
	umulh	pp,  m6, r;	adcs	t11, t11, pp
	umulh	pp,  m7, r;	adcs	t12, t12, pp
	umulh	pp,  m8, r;	adcs	t13, t13, pp
	umulh	pp,  m9, r;	adcs	t14, t14, pp
	umulh	pp, m10, r;	adcs	t15, t15, pp
	umulh	pp, m11, r;	adc	t16, t16, pp

	mul	r, inv, t5

	mul	pp,  m0, r;	adds	 t5,  t5, pp
	mul	pp,  m1, r;	adcs	 t6,  t6, pp
	mul	pp,  m2, r;	adcs	 t7,  t7, pp
	mul	pp,  m3, r;	adcs	 t8,  t8, pp
	mul	pp,  m4, r;	adcs	 t9,  t9, pp
	mul	pp,  m5, r;	adcs	t10, t10, pp
	mul	pp,  m6, r;	adcs	t11, t11, pp
	mul	pp,  m7, r;	adcs	t12, t12, pp
	mul	pp,  m8, r;	adcs	t13, t13, pp
	mul	pp,  m9, r;	adcs	t14, t14, pp
	mul	pp, m10, r;	adcs	t15, t15, pp
	mul	pp, m11, r;	adcs	t16, t16, pp
				adc	t17, XZR, XZR

	umulh	pp,  m0, r;	adds	 t6,  t6, pp
	umulh	pp,  m1, r;	adcs	 t7,  t7, pp
	umulh	pp,  m2, r;	adcs	 t8,  t8, pp
	umulh	pp,  m3, r;	adcs	 t9,  t9, pp
	umulh	pp,  m4, r;	adcs	t10, t10, pp
	umulh	pp,  m5, r;	adcs	t11, t11, pp
	umulh	pp,  m6, r;	adcs	t12, t12, pp
	umulh	pp,  m7, r;	adcs	t13, t13, pp
	umulh	pp,  m8, r;	adcs	t14, t14, pp
	umulh	pp,  m9, r;	adcs	t15, t15, pp
	umulh	pp, m10, r;	adcs	t16, t16, pp
	umulh	pp, m11, r;	adc	t17, t17, pp

	mul	r, inv, t6

	mul	pp,  m0, r;	adds	 t6,  t6, pp
	mul	pp,  m1, r;	adcs	 t7,  t7, pp
	mul	pp,  m2, r;	adcs	 t8,  t8, pp
	mul	pp,  m3, r;	adcs	 t9,  t9, pp
	mul	pp,  m4, r;	adcs	t10, t10, pp
	mul	pp,  m5, r;	adcs	t11, t11, pp
	mul	pp,  m6, r;	adcs	t12, t12, pp
	mul	pp,  m7, r;	adcs	t13, t13, pp
	mul	pp,  m8, r;	adcs	t14, t14, pp
	mul	pp,  m9, r;	adcs	t15, t15, pp
	mul	pp, m10, r;	adcs	t16, t16, pp
	mul	pp, m11, r;	adcs	t17, t17, pp
				adc	t18, XZR, XZR

	umulh	pp,  m0, r;	adds	 t7,  t7, pp
	umulh	pp,  m1, r;	adcs	 t8,  t8, pp
	umulh	pp,  m2, r;	adcs	 t9,  t9, pp
	umulh	pp,  m3, r;	adcs	t10, t10, pp
	umulh	pp,  m4, r;	adcs	t11, t11, pp
	umulh	pp,  m5, r;	adcs	t12, t12, pp
	umulh	pp,  m6, r;	adcs	t13, t13, pp
	umulh	pp,  m7, r;	adcs	t14, t14, pp
	umulh	pp,  m8, r;	adcs	t15, t15, pp
	umulh	pp,  m9, r;	adcs	t16, t16, pp
	umulh	pp, m10, r;	adcs	t17, t17, pp
	umulh	pp, m11, r;	adc	t18, t18, pp

	mul	r, inv, t7

	mul	pp,  m0, r;	adds	 t7,  t7, pp
	mul	pp,  m1, r;	adcs	 t8,  t8, pp
	mul	pp,  m2, r;	adcs	 t9,  t9, pp
	mul	pp,  m3, r;	adcs	t10, t10, pp
	mul	pp,  m4, r;	adcs	t11, t11, pp
	mul	pp,  m5, r;	adcs	t12, t12, pp
	mul	pp,  m6, r;	adcs	t13, t13, pp
	mul	pp,  m7, r;	adcs	t14, t14, pp
	mul	pp,  m8, r;	adcs	t15, t15, pp
	mul	pp,  m9, r;	adcs	t16, t16, pp
	mul	pp, m10, r;	adcs	t17, t17, pp
	mul	pp, m11, r;	adcs	t18, t18, pp
				adc	t19, XZR, XZR

	umulh	pp,  m0, r;	adds	 t8,  t8, pp
	umulh	pp,  m1, r;	adcs	 t9,  t9, pp
	umulh	pp,  m2, r;	adcs	t10, t10, pp
	umulh	pp,  m3, r;	adcs	t11, t11, pp
	umulh	pp,  m4, r;	adcs	t12, t12, pp
	umulh	pp,  m5, r;	adcs	t13, t13, pp
	umulh	pp,  m6, r;	adcs	t14, t14, pp
	umulh	pp,  m7, r;	adcs	t15, t15, pp
	umulh	pp,  m8, r;	adcs	t16, t16, pp
	umulh	pp,  m9, r;	adcs	t17, t17, pp
	umulh	pp, m10, r;	adcs	t18, t18, pp
	umulh	pp, m11, r;	adc	t19, t19, pp

	mul	r, inv, t8

	mul	pp,  m0, r;	adds	 t8,  t8, pp
	mul	pp,  m1, r;	adcs	 t9,  t9, pp
	mul	pp,  m2, r;	adcs	t10, t10, pp
	mul	pp,  m3, r;	adcs	t11, t11, pp
	mul	pp,  m4, r;	adcs	t12, t12, pp
	mul	pp,  m5, r;	adcs	t13, t13, pp
	mul	pp,  m6, r;	adcs	t14, t14, pp
	mul	pp,  m7, r;	adcs	t15, t15, pp
	mul	pp,  m8, r;	adcs	t16, t16, pp
	mul	pp,  m9, r;	adcs	t17, t17, pp
	mul	pp, m10, r;	adcs	t18, t18, pp
	mul	pp, m11, r;	adcs	t19, t19, pp
				adc	t20, XZR, XZR

	umulh	pp,  m0, r;	adds	 t9,  t9, pp
	umulh	pp,  m1, r;	adcs	t10, t10, pp
	umulh	pp,  m2, r;	adcs	t11, t11, pp
	umulh	pp,  m3, r;	adcs	t12, t12, pp
	umulh	pp,  m4, r;	adcs	t13, t13, pp
	umulh	pp,  m5, r;	adcs	t14, t14, pp
	umulh	pp,  m6, r;	adcs	t15, t15, pp
	umulh	pp,  m7, r;	adcs	t16, t16, pp
	umulh	pp,  m8, r;	adcs	t17, t17, pp
	umulh	pp,  m9, r;	adcs	t18, t18, pp
	umulh	pp, m10, r;	adcs	t19, t19, pp
	umulh	pp, m11, r;	adc	t20, t20, pp

	mul	r, inv, t9

	mul	pp,  m0, r;	adds	 t9,  t9, pp
	mul	pp,  m1, r;	adcs	t10, t10, pp
	mul	pp,  m2, r;	adcs	t11, t11, pp
	mul	pp,  m3, r;	adcs	t12, t12, pp
	mul	pp,  m4, r;	adcs	t13, t13, pp
	mul	pp,  m5, r;	adcs	t14, t14, pp
	mul	pp,  m6, r;	adcs	t15, t15, pp
	mul	pp,  m7, r;	adcs	t16, t16, pp
	mul	pp,  m8, r;	adcs	t17, t17, pp
	mul	pp,  m9, r;	adcs	t18, t18, pp
	mul	pp, m10, r;	adcs	t19, t19, pp
	mul	pp, m11, r;	adcs	t20, t20, pp
				adc	t21, XZR, XZR

	umulh	pp,  m0, r;	adds	t10, t10, pp
	umulh	pp,  m1, r;	adcs	t11, t11, pp
	umulh	pp,  m2, r;	adcs	t12, t12, pp
	umulh	pp,  m3, r;	adcs	t13, t13, pp
	umulh	pp,  m4, r;	adcs	t14, t14, pp
	umulh	pp,  m5, r;	adcs	t15, t15, pp
	umulh	pp,  m6, r;	adcs	t16, t16, pp
	umulh	pp,  m7, r;	adcs	t17, t17, pp
	umulh	pp,  m8, r;	adcs	t18, t18, pp
	umulh	pp,  m9, r;	adcs	t19, t19, pp
	umulh	pp, m10, r;	adcs	t20, t20, pp
	umulh	pp, m11, r;	adc	t21, t21, pp

	mul	r, inv, t10

	mul	pp,  m0, r;	adds	t10, t10, pp
	mul	pp,  m1, r;	adcs	t11, t11, pp
	mul	pp,  m2, r;	adcs	t12, t12, pp
	mul	pp,  m3, r;	adcs	t13, t13, pp
	mul	pp,  m4, r;	adcs	t14, t14, pp
	mul	pp,  m5, r;	adcs	t15, t15, pp
	mul	pp,  m6, r;	adcs	t16, t16, pp
	mul	pp,  m7, r;	adcs	t17, t17, pp
	mul	pp,  m8, r;	adcs	t18, t18, pp
	mul	pp,  m9, r;	adcs	t19, t19, pp
	mul	pp, m10, r;	adcs	t20, t20, pp
	mul	pp, m11, r;	adcs	t21, t21, pp
				adc	t22, XZR, XZR

	umulh	pp,  m0, r;	adds	t11, t11, pp
	umulh	pp,  m1, r;	adcs	t12, t12, pp
	umulh	pp,  m2, r;	adcs	t13, t13, pp
	umulh	pp,  m3, r;	adcs	t14, t14, pp
	umulh	pp,  m4, r;	adcs	t15, t15, pp
	umulh	pp,  m5, r;	adcs	t16, t16, pp
	umulh	pp,  m6, r;	adcs	t17, t17, pp
	umulh	pp,  m7, r;	adcs	t18, t18, pp
	umulh	pp,  m8, r;	adcs	t19, t19, pp
	umulh	pp,  m9, r;	adcs	t20, t20, pp
	umulh	pp, m10, r;	adcs	t21, t21, pp
	umulh	pp, m11, r;	adc	t22, t22, pp

	mul	r, inv, t11

	mul	pp,  m0, r;	adds	t11, t11, pp
	mul	pp,  m1, r;	adcs	t12, t12, pp
	mul	pp,  m2, r;	adcs	t13, t13, pp
	mul	pp,  m3, r;	adcs	t14, t14, pp
	mul	pp,  m4, r;	adcs	t15, t15, pp
	mul	pp,  m5, r;	adcs	t16, t16, pp
	mul	pp,  m6, r;	adcs	t17, t17, pp
	mul	pp,  m7, r;	adcs	t18, t18, pp
	mul	pp,  m8, r;	adcs	t19, t19, pp
	mul	pp,  m9, r;	adcs	t20, t20, pp
	mul	pp, m10, r;	adcs	t21, t21, pp
	mul	pp, m11, r;	adcs	t22, t22, pp
				adc	t23, XZR, XZR

	umulh	pp,  m0, r;	adds	t12, t12, pp
	umulh	pp,  m1, r;	adcs	t13, t13, pp
	umulh	pp,  m2, r;	adcs	t14, t14, pp
	umulh	pp,  m3, r;	adcs	t15, t15, pp
	umulh	pp,  m4, r;	adcs	t16, t16, pp
	umulh	pp,  m5, r;	adcs	t17, t17, pp
	umulh	pp,  m6, r;	adcs	t18, t18, pp
	umulh	pp,  m7, r;	adcs	t19, t19, pp
	umulh	pp,  m8, r;	adcs	t20, t20, pp
	umulh	pp,  m9, r;	adcs	t21, t21, pp
	umulh	pp, m10, r;	adcs	t22, t22, pp
	umulh	pp, m11, r;	adc	t23, t23, pp

	// Add top half of product

	ldp	a, b, [SP, #24*8];	adds	t12, t12, a;	adcs	t13, t13, b
	ldp	a, b, [SP, #12*8];	adcs	t14, t14, a;	adcs	t15, t15, b
	ldp	a, b, [SP, #14*8];	adcs	t16, t16, a;	adcs	t17, t17, b
	ldp	a, b, [SP, #16*8];	adcs	t18, t18, a;	adcs	t19, t19, b
	ldp	a, b, [SP, #18*8];	adcs	t20, t20, a;	adcs	t21, t21, b
	ldp	a, b, [SP, #20*8];	adcs	t22, t22, a;	adcs	t23, t23, b

	adc	t24, XZR, XZR

	// Conditionally subtract m:
	//
	// m=t-m
	// if (m >= 0)
	//   t=m
	// store t

	subs	 m0, t12,  m0
	sbcs	 m1, t13,  m1
	sbcs	 m2, t14,  m2
	sbcs	 m3, t15,  m3
	sbcs	 m4, t16,  m4
	sbcs	 m5, t17,  m5
	sbcs	 m6, t18,  m6
	sbcs	 m7, t19,  m7
	sbcs	 m8, t20,  m8
	sbcs	 m9, t21,  m9
	sbcs	m10, t22, m10
	sbcs	m11, t23, m11
	sbcs	XZR, t24, XZR

	csel	t12, t12,  m0, cc
	csel	t13, t13,  m1, cc
	csel	t14, t14,  m2, cc
	csel	t15, t15,  m3, cc
	csel	t16, t16,  m4, cc
	csel	t17, t17,  m5, cc
	csel	t18, t18,  m6, cc
	csel	t19, t19,  m7, cc
	csel	t20, t20,  m8, cc
	csel	t21, t21,  m9, cc
	csel	t22, t22, m10, cc
	csel	t23, t23, m11, cc

	stp	t12, t13, [c, # 0*8]
	stp	t14, t15, [c, # 2*8]
	stp	t16, t17, [c, # 4*8]
	stp	t18, t19, [c, # 6*8]
	stp	t20, t21, [c, # 8*8]
	stp	t22, t23, [c, #10*8]

	// Restore callee-save registers r19-r28

	ldp	x19, x20, [SP, #0*8]
	ldp	x21, x22, [SP, #2*8]
	ldp	x23, x24, [SP, #4*8]
	ldp	x25, x26, [SP, #6*8]
	ldp	x27, x28, [SP, #8*8]

	// Restore stack pointer

	add	SP, SP, #26*8

	ret
