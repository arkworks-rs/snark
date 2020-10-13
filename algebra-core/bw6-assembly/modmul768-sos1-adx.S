// void modmul768(const uint64_t x[12], const uint64_t y[12], const uint64_t m[13], uint64_t z[12])

// m[12] contains the least significant word of the negated inverse of the modulus mod 2^768

#ifdef _WIN64
#	define x	%rcx
#	define y	%rdx
#	define m	%r8
#	define z	%r9
#else
#	define x	%rdi
#	define y	%rsi
#	define m	%rdx
#	define z	%rcx
#endif

#define l	%rax
#define h	%rbx

#define t0	%rcx
#define t1	%rbp
#define t2	%rsi
#define t3	%rdi
#define t4	%r8
#define t5	%r9
#define t6	%r10
#define t7	%r11
#define t8	%r12
#define t9	%r13
#define t10	%r14
#define t11	%r15
#define t12	t0
#define t13	t1
#define t14	t2
#define t15	t3
#define t16	t4
#define t17	t5
#define t18	t6
#define t19	t7
#define t20	t8
#define t21	t9
#define t22	t10
#define t23	t11
#define t24	h

#define x0	0*8(%rsp)
#define x1	1*8(%rsp)
#define x2	2*8(%rsp)
#define x3	3*8(%rsp)
#define x4	4*8(%rsp)
#define x5	5*8(%rsp)
#define x6	6*8(%rsp)
#define x7	7*8(%rsp)
#define x8	8*8(%rsp)
#define x9	9*8(%rsp)
#define x10	10*8(%rsp)
#define x11	11*8(%rsp)

#define m0	x0
#define m1	x1
#define m2	x2
#define m3	x3
#define m4	x4
#define m5	x5
#define m6	x6
#define m7	x7
#define m8	x8
#define m9	x9
#define m10	x10
#define m11	x11
#define inv	12*8(%rsp)

#define z0	13*8(%rsp)
#define z1	14*8(%rsp)
#define z2	15*8(%rsp)
#define z3	16*8(%rsp)
#define z4	17*8(%rsp)
#define z5	18*8(%rsp)
#define z6	19*8(%rsp)
#define z7	20*8(%rsp)
#define z8	21*8(%rsp)
#define z9	22*8(%rsp)
#define z10	23*8(%rsp)
#define z11	24*8(%rsp)
#define z12	25*8(%rsp)
#define z13	26*8(%rsp)
#define z14	27*8(%rsp)
#define z15	28*8(%rsp)
#define z16	29*8(%rsp)
#define z17	30*8(%rsp)
#define z18	31*8(%rsp)
#define z19	32*8(%rsp)
#define z20	33*8(%rsp)
#define z21	34*8(%rsp)
#define z22	35*8(%rsp)
#define z23	36*8(%rsp)

#define y1	z1
#define y2	z2
#define y3	z3
#define y4	z4
#define y5	z5
#define y6	z6
#define y7	z7
#define y8	z8
#define y9	z9
#define y10	z10
#define y11	z11

.text

#ifdef __APPLE__
#define modmul768 _modmul768
#endif

.globl  modmul768
#ifndef __APPLE__
.type   modmul768, @function
#endif

.p2align 6,,15
modmul768:

	// Allocate space on the stack:
	//   13 words for x or m
	//   24 words for x*y (used initially to store y)
	//    6 words for callee-saves
	//    2 words for the m and z pointers

	mov	%rsp, %rax
	sub	$8*(13+24+6+2), %rsp

	// Callee-saves

#ifdef _WIN64
	mov	%rsi, 1*8(%rax)
	mov	%rdi, 2*8(%rax)
#endif
	mov	%rbp, -6*8(%rax)
	mov	%rbx, -5*8(%rax)
	mov	%r12, -4*8(%rax)
	mov	%r13, -3*8(%rax)
	mov	%r14, -2*8(%rax)
	mov	%r15, -1*8(%rax)

	// m and z pointers

	mov	m, -8*8(%rax)
	mov	z, -7*8(%rax)

	// copy x to the stack

	mov	11*8(x), %rax;	mov	%rax, x11
	mov	10*8(x), %rax;	mov	%rax, x10
	mov	 9*8(x), %rax;	mov	%rax, x9
	mov	 8*8(x), %rax;	mov	%rax, x8
	mov	 7*8(x), %rax;	mov	%rax, x7
	mov	 6*8(x), %rax;	mov	%rax, x6
	mov	 5*8(x), %rax;	mov	%rax, x5
	mov	 4*8(x), %rax;	mov	%rax, x4
	mov	 3*8(x), %rax;	mov	%rax, x3
	mov	 2*8(x), %rax;	mov	%rax, x2
	mov	 1*8(x), %rax;	mov	%rax, x1
	mov	 0*8(x), %rax;	mov	%rax, x0

	// copy y to the stack

	mov	11*8(y), %rax;	mov	%rax, y11
	mov	10*8(y), %rax;	mov	%rax, y10
	mov	 9*8(y), %rax;	mov	%rax, y9
	mov	 8*8(y), %rax;	mov	%rax, y8
	mov	 7*8(y), %rax;	mov	%rax, y7
	mov	 6*8(y), %rax;	mov	%rax, y6
	mov	 5*8(y), %rax;	mov	%rax, y5
	mov	 4*8(y), %rax;	mov	%rax, y4
	mov	 3*8(y), %rax;	mov	%rax, y3
	mov	 2*8(y), %rax;	mov	%rax, y2
	mov	 1*8(y), %rax;	mov	%rax, y1
	mov	 0*8(y), %rdx

// y0

	// mov	y0, %rdx

	xor	l, l

	mulx	x0, t0, t1
	mov	t0, z0
	mulx	x2, t2, t3
	mulx	x1, l, h;	adcx	l, t1;	adcx	h, t2
	mulx	x4, t4, t5
	mulx	x3, l, h;	adcx	l, t3;	adcx	h, t4
	mulx	x6, t6, t7
	mulx	x5, l, h;	adcx	l, t5;	adcx	h, t6
	mulx	x8, t8, t9
	mulx	x7, l, h;	adcx	l, t7;	adcx	h, t8
	mulx	x10, t10, t11
	mulx	x9, l, h;	adcx	l, t9;	adcx	h, t10
	mulx	x11, l, t12;	adcx	l, t11;	adc	$0, t12

// y1

	mov	y1, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t1;	adox	h, t2
	mov	t1, z1
	mulx	x1, l, h;	adcx	l, t2;	adcx	h, t3
	mulx	x2, l, h;	adox	l, t3;	adox	h, t4
	mulx	x3, l, h;	adcx	l, t4;	adcx	h, t5
	mulx	x4, l, h;	adox	l, t5;	adox	h, t6
	mulx	x5, l, h;	adcx	l, t6;	adcx	h, t7
	mulx	x6, l, h;	adox	l, t7;	adox	h, t8
	mulx	x7, l, h;	adcx	l, t8;	adcx	h, t9
	mulx	x8, l, h;	adox	l, t9;	adox	h, t10
	mulx	x9, l, h;	adcx	l, t10;	adcx	h, t11
	mulx	x10, l, h;	adox	l, t11;	adox	h, t12
	mulx	x11, l, t13;	adcx	l, t12;	adc	$0, t13

// y2

	mov	y2, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t2;	adox	h, t3
	mov	t2, z2
	mulx	x1, l, h;	adcx	l, t3;	adcx	h, t4
	mulx	x2, l, h;	adox	l, t4;	adox	h, t5
	mulx	x3, l, h;	adcx	l, t5;	adcx	h, t6
	mulx	x4, l, h;	adox	l, t6;	adox	h, t7
	mulx	x5, l, h;	adcx	l, t7;	adcx	h, t8
	mulx	x6, l, h;	adox	l, t8;	adox	h, t9
	mulx	x7, l, h;	adcx	l, t9;	adcx	h, t10
	mulx	x8, l, h;	adox	l, t10;	adox	h, t11
	mulx	x9, l, h;	adcx	l, t11;	adcx	h, t12
	mulx	x10, l, h;	adox	l, t12;	adox	h, t13
	mulx	x11, l, t14;	adcx	l, t13;	adc	$0, t14

// y3

	mov	y3, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t3;	adox	h, t4
	mov	t3, z3
	mulx	x1, l, h;	adcx	l, t4;	adcx	h, t5
	mulx	x2, l, h;	adox	l, t5;	adox	h, t6
	mulx	x3, l, h;	adcx	l, t6;	adcx	h, t7
	mulx	x4, l, h;	adox	l, t7;	adox	h, t8
	mulx	x5, l, h;	adcx	l, t8;	adcx	h, t9
	mulx	x6, l, h;	adox	l, t9;	adox	h, t10
	mulx	x7, l, h;	adcx	l, t10;	adcx	h, t11
	mulx	x8, l, h;	adox	l, t11;	adox	h, t12
	mulx	x9, l, h;	adcx	l, t12;	adcx	h, t13
	mulx	x10, l, h;	adox	l, t13;	adox	h, t14
	mulx	x11, l, t15;	adcx	l, t14;	adc	$0, t15

// y4

	mov	y4, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t4;	adox	h, t5
	mov	t4, z4
	mulx	x1, l, h;	adcx	l, t5;	adcx	h, t6
	mulx	x2, l, h;	adox	l, t6;	adox	h, t7
	mulx	x3, l, h;	adcx	l, t7;	adcx	h, t8
	mulx	x4, l, h;	adox	l, t8;	adox	h, t9
	mulx	x5, l, h;	adcx	l, t9;	adcx	h, t10
	mulx	x6, l, h;	adox	l, t10;	adox	h, t11
	mulx	x7, l, h;	adcx	l, t11;	adcx	h, t12
	mulx	x8, l, h;	adox	l, t12;	adox	h, t13
	mulx	x9, l, h;	adcx	l, t13;	adcx	h, t14
	mulx	x10, l, h;	adox	l, t14;	adox	h, t15
	mulx	x11, l, t16;	adcx	l, t15;	adc	$0, t16

// y5

	mov	y5, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t5;	adox	h, t6
	mov	t5, z5
	mulx	x1, l, h;	adcx	l, t6;	adcx	h, t7
	mulx	x2, l, h;	adox	l, t7;	adox	h, t8
	mulx	x3, l, h;	adcx	l, t8;	adcx	h, t9
	mulx	x4, l, h;	adox	l, t9;	adox	h, t10
	mulx	x5, l, h;	adcx	l, t10;	adcx	h, t11
	mulx	x6, l, h;	adox	l, t11;	adox	h, t12
	mulx	x7, l, h;	adcx	l, t12;	adcx	h, t13
	mulx	x8, l, h;	adox	l, t13;	adox	h, t14
	mulx	x9, l, h;	adcx	l, t14;	adcx	h, t15
	mulx	x10, l, h;	adox	l, t15;	adox	h, t16
	mulx	x11, l, t17;	adcx	l, t16;	adc	$0, t17

// y6

	mov	y6, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t6;	adox	h, t7
	mov	t6, z6
	mulx	x1, l, h;	adcx	l, t7;	adcx	h, t8
	mulx	x2, l, h;	adox	l, t8;	adox	h, t9
	mulx	x3, l, h;	adcx	l, t9;	adcx	h, t10
	mulx	x4, l, h;	adox	l, t10;	adox	h, t11
	mulx	x5, l, h;	adcx	l, t11;	adcx	h, t12
	mulx	x6, l, h;	adox	l, t12;	adox	h, t13
	mulx	x7, l, h;	adcx	l, t13;	adcx	h, t14
	mulx	x8, l, h;	adox	l, t14;	adox	h, t15
	mulx	x9, l, h;	adcx	l, t15;	adcx	h, t16
	mulx	x10, l, h;	adox	l, t16;	adox	h, t17
	mulx	x11, l, t18;	adcx	l, t17;	adc	$0, t18

// y7

	mov	y7, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t7;	adox	h, t8
	mov	t7, z7
	mulx	x1, l, h;	adcx	l, t8;	adcx	h, t9
	mulx	x2, l, h;	adox	l, t9;	adox	h, t10
	mulx	x3, l, h;	adcx	l, t10;	adcx	h, t11
	mulx	x4, l, h;	adox	l, t11;	adox	h, t12
	mulx	x5, l, h;	adcx	l, t12;	adcx	h, t13
	mulx	x6, l, h;	adox	l, t13;	adox	h, t14
	mulx	x7, l, h;	adcx	l, t14;	adcx	h, t15
	mulx	x8, l, h;	adox	l, t15;	adox	h, t16
	mulx	x9, l, h;	adcx	l, t16;	adcx	h, t17
	mulx	x10, l, h;	adox	l, t17;	adox	h, t18
	mulx	x11, l, t19;	adcx	l, t18;	adc	$0, t19

// y8

	mov	y8, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t8;	adox	h, t9
	mov	t8, z8
	mulx	x1, l, h;	adcx	l, t9;	adcx	h, t10
	mulx	x2, l, h;	adox	l, t10;	adox	h, t11
	mulx	x3, l, h;	adcx	l, t11;	adcx	h, t12
	mulx	x4, l, h;	adox	l, t12;	adox	h, t13
	mulx	x5, l, h;	adcx	l, t13;	adcx	h, t14
	mulx	x6, l, h;	adox	l, t14;	adox	h, t15
	mulx	x7, l, h;	adcx	l, t15;	adcx	h, t16
	mulx	x8, l, h;	adox	l, t16;	adox	h, t17
	mulx	x9, l, h;	adcx	l, t17;	adcx	h, t18
	mulx	x10, l, h;	adox	l, t18;	adox	h, t19
	mulx	x11, l, t20;	adcx	l, t19;	adc	$0, t20

// y9

	mov	y9, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t9;	adox	h, t10
	mov	t9, z9
	mulx	x1, l, h;	adcx	l, t10;	adcx	h, t11
	mulx	x2, l, h;	adox	l, t11;	adox	h, t12
	mulx	x3, l, h;	adcx	l, t12;	adcx	h, t13
	mulx	x4, l, h;	adox	l, t13;	adox	h, t14
	mulx	x5, l, h;	adcx	l, t14;	adcx	h, t15
	mulx	x6, l, h;	adox	l, t15;	adox	h, t16
	mulx	x7, l, h;	adcx	l, t16;	adcx	h, t17
	mulx	x8, l, h;	adox	l, t17;	adox	h, t18
	mulx	x9, l, h;	adcx	l, t18;	adcx	h, t19
	mulx	x10, l, h;	adox	l, t19;	adox	h, t20
	mulx	x11, l, t21;	adcx	l, t20;	adc	$0, t21

// y10

	mov	y10, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t10;	adox	h, t11
	mov	t10, z10
	mulx	x1, l, h;	adcx	l, t11;	adcx	h, t12
	mulx	x2, l, h;	adox	l, t12;	adox	h, t13
	mulx	x3, l, h;	adcx	l, t13;	adcx	h, t14
	mulx	x4, l, h;	adox	l, t14;	adox	h, t15
	mulx	x5, l, h;	adcx	l, t15;	adcx	h, t16
	mulx	x6, l, h;	adox	l, t16;	adox	h, t17
	mulx	x7, l, h;	adcx	l, t17;	adcx	h, t18
	mulx	x8, l, h;	adox	l, t18;	adox	h, t19
	mulx	x9, l, h;	adcx	l, t19;	adcx	h, t20
	mulx	x10, l, h;	adox	l, t20;	adox	h, t21
	mulx	x11, l, t22;	adcx	l, t21;	adc	$0, t22

// y11

	mov	y11, %rdx
	xor	l, l

	mulx	x0, l, h;	adox	l, t11;	adox	h, t12
	mov	t11, z11
	mulx	x1, l, h;	adcx	l, t12;	adcx	h, t13
	mulx	x2, l, h;	adox	l, t13;	adox	h, t14
	mulx	x3, l, h;	adcx	l, t14;	adcx	h, t15
	mulx	x4, l, h;	adox	l, t15;	adox	h, t16
	mulx	x5, l, h;	adcx	l, t16;	adcx	h, t17
	mulx	x6, l, h;	adox	l, t17;	adox	h, t18
	mulx	x7, l, h;	adcx	l, t18;	adcx	h, t19
	mulx	x8, l, h;	adox	l, t19;	adox	h, t20
	mulx	x9, l, h;	adcx	l, t20;	adcx	h, t21
	mulx	x10, l, h;	adox	l, t21;	adox	h, t22
	mulx	x11, l, t23;	adcx	l, t22;	adc	$0, t23

	// Copy m to the stack, overwriting x

	mov	8*(13+24)(%rsp), h
	mov	12*8(h), %rdx;	mov	%rdx, inv
	mov	11*8(h), %rax;	mov	%rax, m11
	mov	10*8(h), %rax;	mov	%rax, m10
	mov	 9*8(h), %rax;	mov	%rax, m9
	mov	 8*8(h), %rax;	mov	%rax, m8
	mov	 7*8(h), %rax;	mov	%rax, m7
	mov	 6*8(h), %rax;	mov	%rax, m6
	mov	 5*8(h), %rax;	mov	%rax, m5
	mov	 4*8(h), %rax;	mov	%rax, m4
	mov	 3*8(h), %rax;	mov	%rax, m3
	mov	 2*8(h), %rax;	mov	%rax, m2
	mov	 1*8(h), %rax;	mov	%rax, m1
	mov	 0*8(h), %rax;	mov	%rax, m0

	// Write out the top half of x*y to the stack, load the low half back in

	mov	t12, z12;	mov	z0, t0
	mov	t13, z13;	mov	z1, t1
	mov	t14, z14;	mov	z2, t2
	mov	t15, z15;	mov	z3, t3
	mov	t16, z16;	mov	z4, t4
	mov	t17, z17;	mov	z5, t5
	mov	t18, z18;	mov	z6, t6
	mov	t19, z19;	mov	z7, t7
	mov	t20, z20;	mov	z8, t8
	mov	t21, z21;	mov	z9, t9
	mov	t22, z22;	mov	z10, t10
	mov	t23, z23;	mov	z11, t11

////////////////////////////////////////////////////////////////
// Reduction
////////////////////////////////////////////////////////////////

// z0

	//mov	inv, %rdx
	mulx	t0, %rdx, h

	mulx	m1, l, h;	adox	l, t1;	adox	h, t2
	mulx	m0, l, h;	adcx	l, t0;	adcx	h, t1
	mulx	m3, l, h;	adox	l, t3;	adox	h, t4
	mulx	m2, l, h;	adcx	l, t2;	adcx	h, t3
	mulx	m5, l, h;	adox	l, t5;	adox	h, t6
	mulx	m4, l, h;	adcx	l, t4;	adcx	h, t5
	mulx	m7, l, h;	adox	l, t7;	adox	h, t8
	mulx	m6, l, h;	adcx	l, t6;	adcx	h, t7
	mulx	m9, l, h;	adox	l, t9;	adox	h, t10
	mulx	m8, l, h;	adcx	l, t8;	adcx	h, t9
	mulx	m11, l, h;	adox	l, t11;	adox	h, t12
	mulx	m10, l, h;	adcx	l, t10;	adcx	h, t11;	adc	$0, t12

// z1

	mov	inv, %rdx
	mulx	t1, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t2;	adox	h, t3
	mulx	m0, l, h;	adcx	l, t1;	adcx	h, t2
	mulx	m3, l, h;	adox	l, t4;	adox	h, t5
	mulx	m2, l, h;	adcx	l, t3;	adcx	h, t4
	mulx	m5, l, h;	adox	l, t6;	adox	h, t7
	mulx	m4, l, h;	adcx	l, t5;	adcx	h, t6
	mulx	m7, l, h;	adox	l, t8;	adox	h, t9
	mulx	m6, l, h;	adcx	l, t7;	adcx	h, t8
	mulx	m9, l, h;	adox	l, t10;	adox	h, t11
	mulx	m8, l, h;	adcx	l, t9;	adcx	h, t10
	mulx	m11, l, h;	adox	l, t12;	adox	h, t13
	mulx	m10, l, h;	adcx	l, t11;	adcx	h, t12;	adc	$0, t13

// z2

	mov	inv, %rdx
	mulx	t2, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t3;	adox	h, t4
	mulx	m0, l, h;	adcx	l, t2;	adcx	h, t3
	mulx	m3, l, h;	adox	l, t5;	adox	h, t6
	mulx	m2, l, h;	adcx	l, t4;	adcx	h, t5
	mulx	m5, l, h;	adox	l, t7;	adox	h, t8
	mulx	m4, l, h;	adcx	l, t6;	adcx	h, t7
	mulx	m7, l, h;	adox	l, t9;	adox	h, t10
	mulx	m6, l, h;	adcx	l, t8;	adcx	h, t9
	mulx	m9, l, h;	adox	l, t11;	adox	h, t12
	mulx	m8, l, h;	adcx	l, t10;	adcx	h, t11
	mulx	m11, l, h;	adox	l, t13;	adox	h, t14
	mulx	m10, l, h;	adcx	l, t12;	adcx	h, t13;	adc	$0, t14

// z3

	mov	inv, %rdx
	mulx	t3, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t4;	adox	h, t5
	mulx	m0, l, h;	adcx	l, t3;	adcx	h, t4
	mulx	m3, l, h;	adox	l, t6;	adox	h, t7
	mulx	m2, l, h;	adcx	l, t5;	adcx	h, t6
	mulx	m5, l, h;	adox	l, t8;	adox	h, t9
	mulx	m4, l, h;	adcx	l, t7;	adcx	h, t8
	mulx	m7, l, h;	adox	l, t10;	adox	h, t11
	mulx	m6, l, h;	adcx	l, t9;	adcx	h, t10
	mulx	m9, l, h;	adox	l, t12;	adox	h, t13
	mulx	m8, l, h;	adcx	l, t11;	adcx	h, t12
	mulx	m11, l, h;	adox	l, t14;	adox	h, t15
	mulx	m10, l, h;	adcx	l, t13;	adcx	h, t14;	adc	$0, t15

// z4

	mov	inv, %rdx
	mulx	t4, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t5;	adox	h, t6
	mulx	m0, l, h;	adcx	l, t4;	adcx	h, t5
	mulx	m3, l, h;	adox	l, t7;	adox	h, t8
	mulx	m2, l, h;	adcx	l, t6;	adcx	h, t7
	mulx	m5, l, h;	adox	l, t9;	adox	h, t10
	mulx	m4, l, h;	adcx	l, t8;	adcx	h, t9
	mulx	m7, l, h;	adox	l, t11;	adox	h, t12
	mulx	m6, l, h;	adcx	l, t10;	adcx	h, t11
	mulx	m9, l, h;	adox	l, t13;	adox	h, t14
	mulx	m8, l, h;	adcx	l, t12;	adcx	h, t13
	mulx	m11, l, h;	adox	l, t15;	adox	h, t16
	mulx	m10, l, h;	adcx	l, t14;	adcx	h, t15;	adc	$0, t16

// z5

	mov	inv, %rdx
	mulx	t5, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t6;	adox	h, t7
	mulx	m0, l, h;	adcx	l, t5;	adcx	h, t6
	mulx	m3, l, h;	adox	l, t8;	adox	h, t9
	mulx	m2, l, h;	adcx	l, t7;	adcx	h, t8
	mulx	m5, l, h;	adox	l, t10;	adox	h, t11
	mulx	m4, l, h;	adcx	l, t9;	adcx	h, t10
	mulx	m7, l, h;	adox	l, t12;	adox	h, t13
	mulx	m6, l, h;	adcx	l, t11;	adcx	h, t12
	mulx	m9, l, h;	adox	l, t14;	adox	h, t15
	mulx	m8, l, h;	adcx	l, t13;	adcx	h, t14
	mulx	m11, l, h;	adox	l, t16;	adox	h, t17
	mulx	m10, l, h;	adcx	l, t15;	adcx	h, t16;	adc	$0, t17

// z6

	mov	inv, %rdx
	mulx	t6, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t7;	adox	h, t8
	mulx	m0, l, h;	adcx	l, t6;	adcx	h, t7
	mulx	m3, l, h;	adox	l, t9;	adox	h, t10
	mulx	m2, l, h;	adcx	l, t8;	adcx	h, t9
	mulx	m5, l, h;	adox	l, t11;	adox	h, t12
	mulx	m4, l, h;	adcx	l, t10;	adcx	h, t11
	mulx	m7, l, h;	adox	l, t13;	adox	h, t14
	mulx	m6, l, h;	adcx	l, t12;	adcx	h, t13
	mulx	m9, l, h;	adox	l, t15;	adox	h, t16
	mulx	m8, l, h;	adcx	l, t14;	adcx	h, t15
	mulx	m11, l, h;	adox	l, t17;	adox	h, t18
	mulx	m10, l, h;	adcx	l, t16;	adcx	h, t17;	adc	$0, t18

// z7

	mov	inv, %rdx
	mulx	t7, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t8;	adox	h, t9
	mulx	m0, l, h;	adcx	l, t7;	adcx	h, t8
	mulx	m3, l, h;	adox	l, t10;	adox	h, t11
	mulx	m2, l, h;	adcx	l, t9;	adcx	h, t10
	mulx	m5, l, h;	adox	l, t12;	adox	h, t13
	mulx	m4, l, h;	adcx	l, t11;	adcx	h, t12
	mulx	m7, l, h;	adox	l, t14;	adox	h, t15
	mulx	m6, l, h;	adcx	l, t13;	adcx	h, t14
	mulx	m9, l, h;	adox	l, t16;	adox	h, t17
	mulx	m8, l, h;	adcx	l, t15;	adcx	h, t16
	mulx	m11, l, h;	adox	l, t18;	adox	h, t19
	mulx	m10, l, h;	adcx	l, t17;	adcx	h, t18;	adc	$0, t19

// z8

	mov	inv, %rdx
	mulx	t8, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t9;	adox	h, t10
	mulx	m0, l, h;	adcx	l, t8;	adcx	h, t9
	mulx	m3, l, h;	adox	l, t11;	adox	h, t12
	mulx	m2, l, h;	adcx	l, t10;	adcx	h, t11
	mulx	m5, l, h;	adox	l, t13;	adox	h, t14
	mulx	m4, l, h;	adcx	l, t12;	adcx	h, t13
	mulx	m7, l, h;	adox	l, t15;	adox	h, t16
	mulx	m6, l, h;	adcx	l, t14;	adcx	h, t15
	mulx	m9, l, h;	adox	l, t17;	adox	h, t18
	mulx	m8, l, h;	adcx	l, t16;	adcx	h, t17
	mulx	m11, l, h;	adox	l, t19;	adox	h, t20
	mulx	m10, l, h;	adcx	l, t18;	adcx	h, t19;	adc	$0, t20

// z9

	mov	inv, %rdx
	mulx	t9, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t10;	adox	h, t11
	mulx	m0, l, h;	adcx	l, t9;	adcx	h, t10
	mulx	m3, l, h;	adox	l, t12;	adox	h, t13
	mulx	m2, l, h;	adcx	l, t11;	adcx	h, t12
	mulx	m5, l, h;	adox	l, t14;	adox	h, t15
	mulx	m4, l, h;	adcx	l, t13;	adcx	h, t14
	mulx	m7, l, h;	adox	l, t16;	adox	h, t17
	mulx	m6, l, h;	adcx	l, t15;	adcx	h, t16
	mulx	m9, l, h;	adox	l, t18;	adox	h, t19
	mulx	m8, l, h;	adcx	l, t17;	adcx	h, t18
	mulx	m11, l, h;	adox	l, t20;	adox	h, t21
	mulx	m10, l, h;	adcx	l, t19;	adcx	h, t20;	adc	$0, t21

// z10

	mov	inv, %rdx
	mulx	t10, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t11;	adox	h, t12
	mulx	m0, l, h;	adcx	l, t10;	adcx	h, t11
	mulx	m3, l, h;	adox	l, t13;	adox	h, t14
	mulx	m2, l, h;	adcx	l, t12;	adcx	h, t13
	mulx	m5, l, h;	adox	l, t15;	adox	h, t16
	mulx	m4, l, h;	adcx	l, t14;	adcx	h, t15
	mulx	m7, l, h;	adox	l, t17;	adox	h, t18
	mulx	m6, l, h;	adcx	l, t16;	adcx	h, t17
	mulx	m9, l, h;	adox	l, t19;	adox	h, t20
	mulx	m8, l, h;	adcx	l, t18;	adcx	h, t19
	mulx	m11, l, h;	adox	l, t21;	adox	h, t22
	mulx	m10, l, h;	adcx	l, t20;	adcx	h, t21;	adc	$0, t22

// z11

	mov	inv, %rdx
	mulx	t11, %rdx, h
	xor	h, h

	mulx	m1, l, h;	adox	l, t12;	adox	h, t13
	mulx	m0, l, h;	adcx	l, t11;	adcx	h, t12
	mulx	m3, l, h;	adox	l, t14;	adox	h, t15
	mulx	m2, l, h;	adcx	l, t13;	adcx	h, t14
	mulx	m5, l, h;	adox	l, t16;	adox	h, t17
	mulx	m4, l, h;	adcx	l, t15;	adcx	h, t16
	mulx	m7, l, h;	adox	l, t18;	adox	h, t19
	mulx	m6, l, h;	adcx	l, t17;	adcx	h, t18
	mulx	m9, l, h;	adox	l, t20;	adox	h, t21
	mulx	m8, l, h;	adcx	l, t19;	adcx	h, t20
	mulx	m11, l, h;	adox	l, t22;	adox	h, t23
	mulx	m10, l, h;	adcx	l, t21;	adcx	h, t22;	adc	$0, t23

	xor	t24, t24;	lea	8*(13+24+6+2)(%rsp), %rax

	add	z12, t12
	adc	z13, t13
	adc	z14, t14
	adc	z15, t15
	adc	z16, t16
	adc	z17, t17
	adc	z18, t18
	adc	z19, t19
	adc	z20, t20
	adc	z21, t21
	adc	z22, t22
	adc	z23, t23
	adc	$0, t24;	mov	-7*8(%rax), %rdx

	// Conditional subtraction of m

	mov	t12, z0;	sub	m0, t12
	mov	t13, z1;	sbb	m1, t13
	mov	t14, z2;	sbb	m2, t14
	mov	t15, z3;	sbb	m3, t15
	mov	t16, z4;	sbb	m4, t16
	mov	t17, z5;	sbb	m5, t17
	mov	t18, z6;	sbb	m6, t18
	mov	t19, z7;	sbb	m7, t19
	mov	t20, z8;	sbb	m8, t20
	mov	t21, z9;	sbb	m9, t21
	mov	t22, z10;	sbb	m10, t22
	mov	t23, z11;	sbb	m11, t23
				sbb	$0, t24

	cmovc	z0, t12
	cmovc	z1, t13
	cmovc	z2, t14
	cmovc	z3, t15
	cmovc	z4, t16
	cmovc	z5, t17
	cmovc	z6, t18
	cmovc	z7, t19
	cmovc	z8, t20
	cmovc	z9, t21
	cmovc	z10, t22
	cmovc	z11, t23

	mov	t12,  0*8(%rdx)
	mov	t13,  1*8(%rdx)
	mov	t14,  2*8(%rdx)
	mov	t15,  3*8(%rdx)
	mov	t16,  4*8(%rdx)
	mov	t17,  5*8(%rdx)
	mov	t18,  6*8(%rdx);	mov	-6*8(%rax), %rbp
	mov	t19,  7*8(%rdx);	mov	-5*8(%rax), %rbx
	mov	t20,  8*8(%rdx);	mov	-4*8(%rax), %r12
	mov	t21,  9*8(%rdx);	mov	-3*8(%rax), %r13
	mov	t22, 10*8(%rdx);	mov	-2*8(%rax), %r14
	mov	t23, 11*8(%rdx);	mov	-1*8(%rax), %r15

	add	$8*(13+24+6+2), %rsp 

#ifdef _WIN64
	mov	1*8(%rax), %rsi
	mov	2*8(%rax), %rdi
#endif
	ret

// vim: noet ts=8 sw=8
