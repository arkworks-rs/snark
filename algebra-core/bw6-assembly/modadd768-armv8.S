// void modadd768(const uint64_t a[12], const uint64_t b[12], const uint64_t m[12], uint64_t c[12])

// Function parameters

#define a   x0
#define b   x1
#define m   x2
#define c   x3

// Internal register roles

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
#define a12 a

#define b0  x16
#define b1  x17
#define b2  x19
#define b3  x20
#define b4  x21
#define b5  x22
#define b6  x23
#define b7  x24
#define b8  x25
#define b9  x26
#define b10 x27
#define b11 x28

#define m0  b0
#define m1  b1
#define m2  b2
#define m3  b3
#define m4  b4
#define m5  b5
#define m6  b6
#define m7  b7
#define m8  b8
#define m9  b9
#define m10 b10
#define m11 b11

#ifdef __APPLE__
#define modadd768 _modadd768
#endif

.global modadd768

.text
.p2align 6,,63

modadd768:
	// Allocate space for storing the 10 callee-save registers

	sub	SP, SP, #10*8

	// Preserve r19-r28		// Load a			// Load b			// a += b

	stp	x19, x20, [SP, # 0*8];	ldp	 a0,  a1, [a, # 0*8];	ldp	 b0,  b1, [b, # 0*8];	adds	 a0,  a0,  b0;	adcs	 a1,  a1,  b1
	stp	x21, x22, [SP, # 2*8];	ldp	 a2,  a3, [a, # 2*8];	ldp	 b2,  b3, [b, # 2*8];	adcs	 a2,  a2,  b2;	adcs	 a3,  a3,  b3
	stp	x23, x24, [SP, # 4*8];	ldp	 a4,  a5, [a, # 4*8];	ldp	 b4,  b5, [b, # 4*8];	adcs	 a4,  a4,  b4;	adcs	 a5,  a5,  b5
	stp	x25, x26, [SP, # 6*8];	ldp	 a6,  a7, [a, # 6*8];	ldp	 b6,  b7, [b, # 6*8];	adcs	 a6,  a6,  b6;	adcs	 a7,  a7,  b7
	stp	x27, x28, [SP, # 8*8];	ldp	 a8,  a9, [a, # 8*8];	ldp	 b8,  b9, [b, # 8*8];	adcs	 a8,  a8,  b8;	adcs	 a9,  a9,  b9
					ldp	a10, a11, [a, #10*8];	ldp	b10, b11, [b, #10*8];	adcs	a10, a10, b10;	adcs	a11, a11, b11
													adc	a12, XZR, XZR
	// Load m			// m = a-m

	ldp	 m0,  m1, [m, # 0*8];	subs	 m0,  a0,  m0;	sbcs	 m1,  a1,  m1;
	ldp	 m2,  m3, [m, # 2*8];	sbcs	 m2,  a2,  m2;	sbcs	 m3,  a3,  m3;
	ldp	 m4,  m5, [m, # 4*8];	sbcs	 m4,  a4,  m4;	sbcs	 m5,  a5,  m5;
	ldp	 m6,  m7, [m, # 6*8];	sbcs	 m6,  a6,  m6;	sbcs	 m7,  a7,  m7;
	ldp	 m8,  m9, [m, # 8*8];	sbcs	 m8,  a8,  m8;	sbcs	 m9,  a9,  m9;
	ldp	m10, m11, [m, #10*8];	sbcs	m10, a10, m10;	sbcs	m11, a11, m11;
					sbcs	XZR, a12, XZR

	// a = m<0 ? a : m

	csel	 a0,  a0,  m0, cc
	csel	 a1,  a1,  m1, cc
	csel	 a2,  a2,  m2, cc
	csel	 a3,  a3,  m3, cc
	csel	 a4,  a4,  m4, cc
	csel	 a5,  a5,  m5, cc
	csel	 a6,  a6,  m6, cc
	csel	 a7,  a7,  m7, cc
	csel	 a8,  a8,  m8, cc
	csel	 a9,  a9,  m9, cc
	csel	a10, a10, m10, cc
	csel	a11, a11, m11, cc

	// Store result in c		// Restore callee-save registers r19-r28, restore stack pointer

	stp	 a0,  a1, [c, # 0*8];	ldp	x19, x20, [SP, #0*8]
	stp	 a2,  a3, [c, # 2*8];	ldp	x21, x22, [SP, #2*8]
	stp	 a4,  a5, [c, # 4*8];	ldp	x23, x24, [SP, #4*8]
	stp	 a6,  a7, [c, # 6*8];	ldp	x25, x26, [SP, #6*8]
	stp	 a8,  a9, [c, # 8*8];	ldp	x27, x28, [SP, #8*8]
	stp	a10, a11, [c, #10*8];	add	SP, SP, #10*8

	ret
