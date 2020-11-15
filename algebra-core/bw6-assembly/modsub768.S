// void modsub768(const uint64_t x[12], const uint64_t y[12], const uint64_t m[13], uint64_t z[12])

#ifdef _WIN64
#	define x	%rcx
#	define y	%rdx
#	define m	%r8
#	define z	%r9

#	define t2	%rdi
#	define t3	%rsi
#else
#	define x	%rdi
#	define y	%rsi
#	define m	%rdx
#	define z	%rcx

#	define t2	%r9
#	define t3	%r8
#endif

#define t0	%r11
#define t1	%r10
#define t4	%r15
#define t5	%r14

#define t6	%r13
#define t7	%r12
#define t8	%rbx
#define t9	%rax
#define t10	%rbp
#define t11	x
#define t12	z

#define  y0	 0*8(y)
#define  y1	 1*8(y)
#define  y2	 2*8(y)
#define  y3	 3*8(y)
#define  y4	 4*8(y)
#define  y5	 5*8(y)
#define  y6	 6*8(y)
#define  y7	 7*8(y)
#define  y8	 8*8(y)
#define  y9	 9*8(y)
#define y10	10*8(y)
#define y11	11*8(y)

#define  m0	 0*8(m)
#define  m1	 1*8(m)
#define  m2	 2*8(m)
#define  m3	 3*8(m)
#define  m4	 4*8(m)
#define  m5	 5*8(m)
#define  m6	 6*8(m)
#define  m7	 7*8(m)
#define  m8	 8*8(m)
#define  m9	 9*8(m)
#define m10	10*8(m)
#define m11	11*8(m)

// We only use these after replacing y with z

#define  z0	 0*8(y)
#define  z1	 1*8(y)
#define  z2	 2*8(y)
#define  z3	 3*8(y)
#define  z4	 4*8(y)
#define  z5	 5*8(y)
#define  z6	 6*8(y)
#define  z7	 7*8(y)
#define  z8	 8*8(y)
#define  z9	 9*8(y)
#define z10	10*8(y)
#define z11	11*8(y)

.text

#ifdef __APPLE__
#define modsub768 _modsub768
#endif

.globl  modsub768
#ifndef __APPLE__
#ifndef _WIN64
.type   modsub768, @function
#endif
#endif

.p2align 6,,15
modsub768:

	// Callee-saves

#ifdef _WIN64
	mov	%rsi, 1*8(%rsp)
	mov	%rdi, 2*8(%rsp)
#endif
			// Load x
	push	%r15;	mov	 0*8(x),  t0;	mov	 1*8(x),  t1
	push	%r14;	mov	 2*8(x),  t2;	mov	 3*8(x),  t3
	push	%r13;	mov	 4*8(x),  t4;	mov	 5*8(x),  t5
	push	%r12;	mov	 6*8(x),  t6;	mov	 7*8(x),  t7
	push	%rbx;	mov	 8*8(x),  t8;	mov	 9*8(x),  t9
	push	%rbp;	mov	10*8(x), t10;	mov	11*8(x), t11
	push	z

	xor	t12, t12
	sub	 y0,  t0
	sbb	 y1,  t1
	sbb	 y2,  t2
	sbb	 y3,  t3
	sbb	 y4,  t4
	sbb	 y5,  t5
	sbb	 y6,  t6
	sbb	 y7,  t7
	sbb	 y8,  t8
	sbb	 y9,  t9
	sbb	y10, t10
	sbb	y11, t11
	sbb	 $0, t12	// -1 if y>x

	// no more need for y. load z in its place

	pop	y

	// Conditional addition of m

	mov	 t0,  z0;	add	 m0,  t0
	mov	 t1,  z1;	adc	 m1,  t1
	mov	 t2,  z2;	adc	 m2,  t2
	mov	 t3,  z3;	adc	 m3,  t3
	mov	 t4,  z4;	adc	 m4,  t4
	mov	 t5,  z5;	adc	 m5,  t5
	mov	 t6,  z6;	adc	 m6,  t6
	mov	 t7,  z7;	adc	 m7,  t7
	mov	 t8,  z8;	adc	 m8,  t8
	mov	 t9,  z9;	adc	 m9,  t9
	mov	t10, z10;	adc	m10, t10
	mov	t11, z11;	adc	m11, t11

	add	 $1, t12	// sets carry if adding m is needed

	cmovnc	 z0,  t0
	cmovnc	 z1,  t1
	cmovnc	 z2,  t2
	cmovnc	 z3,  t3
	cmovnc	 z4,  t4
	cmovnc	 z5,  t5
	cmovnc	 z6,  t6
	cmovnc	 z7,  t7
	cmovnc	 z8,  t8
	cmovnc	 z9,  t9
	cmovnc	z10, t10
	cmovnc	z11, t11

	mov	 t0,  z0
	mov	 t1,  z1
	mov	 t2,  z2
	mov	 t3,  z3
	mov	 t4,  z4
	mov	 t5,  z5
	mov	 t6,  z6
	mov	 t7,  z7
	mov	 t8,  z8
	mov	 t9,  z9
	mov	t10, z10
	mov	t11, z11

#ifdef _WIN64
	mov	7*8(%rsp), %rsi
	mov	8*8(%rsp), %rdi
#endif
			// Load x
	pop	%rbp
	pop	%rbx
	pop	%r12
	pop	%r13
	pop	%r14
	pop	%r15

	ret

// vim: noet ts=8 sw=8
