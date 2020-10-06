// void modmul768(const uint64_t x[12], const uint64_t y[12], const uint64_t m[13], uint64_t z[24])

// m[12] contains the least significant word of the negated inverse of the modulus mod 2^768

/*
  For compact encoding of mov instructions:
  - prefer low registers (r[a-ds][ipx]) for pointers
  But:
  - reserve rdx for multiplication
  - prefer other low registers over rbp for pointers
  - prefer to keep function parameters in their original registers
*/

#ifdef _WIN64
#	 define x	%rcx
#	 define y	%rbx
#	 define m	%r8
#	 define z	%r9

#	 define r0	%rsi
#	 define t0	%rdi

#	 define l0	%rax
#	 define h0	%rbp
#	 define l1	%r10
#	 define h1	%r11
#else
#	 define x	%rdi
#	 define y	%rsi
#	 define m	%rax
#	 define z	%rcx

#	 define r0	%rbx
#	 define t0	%rbp

#	 define l0	%r8
#	 define h0	%r9
#	 define l1	%r10
#	 define h1	%r11
#endif

#define l2	l0
#define h2	h0
#define l3	l1
#define h3	h1

#define l4	l2
#define h4	h2
#define l5	l3
#define h5	h3

#define l6	l4
#define h6	h4
#define l7	l5
#define h7	h5

#define l8	l6
#define h8	h6
#define l9	l7
#define h9	h7

#define l10	l8
#define h10	h8
#define l11	l9
#define h11	h9

#define l12	l10
#define h12	h10

#define z0	%r12
#define z1	%r13
#define z2	%r14
#define z3	z1
#define z4	z2
#define z5	z3
#define z6	z4
#define z7	z5
#define z8	z6
#define z9	z7
#define z10	z8
#define z11	z9
#define z12	z10
#define z13	%r15

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
	push	%rbx

	// %rdx is used for multiplicands
#ifdef _WIN64
	mov	%rdx, y
#else
	mov	%rdx, m
#endif

	mov	0*8(x), %rdx

	push	%rbp

	mov	0*8(y), l0	// load y[0]	-> l0
	mov	1*8(y), l1	// load y[1]	-> l1

	push	%r12

	mulx	l0, z0, h0	// rdx * l0	-> h0:z0
	mov	2*8(y), l2	// load y[2]	-> l2
	mulx	l1, l1, h1	// rdx * l1	-> h1:l1
	add	h0, l1		// add h0,l0	-> l1
	mov	l1, 0*8(z)	// store l1	-> z[0]
	mov	3*8(y), l3	// load y[3]	-> l3
	mulx	l2, l2, h2	// rdx * l2	-> h2:l2
	adc	h1, l2		// add h1,l1	-> l2
	mov	l2, 1*8(z)	// store l2	-> z[1]

	mulx	l3, l3, h3	// rdx * l3	-> h3:l3
	adc	h2, l3		// add h2,l2	-> l3
	mov	l3, 2*8(z)	// store l3	-> z[2]

	mov	4*8(y), l4	// load y[4]	-> l4
	mulx	l4, l4, h4	// rdx * l4	-> h4:l4
	adc	h3, l4		// add h3,l3	-> l4
	mov	l4, 3*8(z)	// store l4	-> z[3]

#ifdef _WIN64
	push	%rsi
#endif

	mov	5*8(y), l5	// load y[5]	-> l5
	mulx	l5, l5, h5	// rdx * l5	-> h5:l5
	adc	h4, l5		// add h4,l4	-> l5
	mov	l5, 4*8(z)	// store l5	-> z[4]

#ifdef _WIN64
	push	%rdi
#endif

	mov	z0, %rdx	// mov z0	-> rdx
	mulx	12*8(m), r0, h4	// mulx inv	-> h4:r0 // h4 discarded
	mov	0*8(x), %rdx

	push	%r13

	mov	6*8(y), l6	// load y[6]	-> l6
	mulx	l6, l6, h6	// rdx * l6	-> h6:l6
	adc	h5, l6		// add h5,l5	-> l6
	mov	l6, 5*8(z)	// store l6	-> z[5]

	mov	7*8(y), l7	// load y[7]	-> l7
	mulx	l7, l7, h7	// rdx * l7	-> h7:l7
	adc	h6, l7		// add h6,l6	-> l7
	mov	l7, 6*8(z)	// store l7	-> z[6]

	push	%r14

	mov	8*8(y), l8	// load y[8]	-> l8
	mulx	l8, l8, h8	// rdx * l8	-> h8:l8
	adc	h7, l8		// add h7,l7	-> l8
	mov	l8, 7*8(z)	// store l8	-> z[7]

	mov	9*8(y), l9	// load y[9]	-> l9
	mulx	l9, l9, h9	// rdx * l9	-> h9:l9
	adc	h8, l9		// add h8,l8	-> l9
	mov	l9, 8*8(z)	// store l9	-> z[8]

	push	%r15

	mov	10*8(y), l10	// load y[10]	-> l10
	mulx	l10, l10, h10	// rdx * l10	-> h10:l10
	adc	h9, l10		// add h9,l9	-> l10
	mov	l10, 9*8(z)	// store l10	-> z[9]

	mov	11*8(y), l11	// load y[11]	-> l11
	mulx	l11, l11, h11	// rdx * l11	-> h11:l11
	adc	h10, l11	// add h10,l10	-> l11
	mov	l11, 10*8(z)	// store l11	-> z[10]

	adc	$0, h11		// adc 0,h11	-> h11
	mov	h11, 11*8(z)	// store h11	-> z[11]

	// Reduction

	mov	r0, %rdx	// r0		-> rdx
	xor	z13, z13	// clear flags

	mov	0*8(m), l0	// load m[0]	-> l0
	mulx	l0, l0, h0	// mulx l0	-> h0:l0
	adcx	l0, z0		// adcx l0,z0	-> z0	// 0, ignore

	mov	0*8(z), z1	// load z[0]	-> z1
	adox	h0, z1		// adox h0,z1	-> z1
	mulx	1*8(m), l1, h1	// mulx m[1]	-> h1:l1
	adcx	l1, z1		// adcx l1,z1	-> z1
	mov	z1, 0*8(z)	// store z1	-> z[0]

	mov	1*8(z), z2	// load z[1]	-> z2
	adox	h1, z2		// adox h1,z2	-> z2
	mulx	2*8(m), l2, h2	// mulx m[2]	-> h2:l2
	adcx	l2, z2		// adcx l2,z2	-> z2
	mov	z2, 1*8(z)	// store z2	-> z[1]

	mov	2*8(z), z3	// load z[2]	-> z3
	adox	h2, z3		// adox h2,z3	-> z3
	mulx	3*8(m), l3, h3	// mulx m[3]	-> h3:l3
	adcx	l3, z3		// adcx l3,z3	-> z3
	mov	z3, 2*8(z)	// store z3	-> z[2]

	mov	3*8(z), z4	// load z[3]	-> z4
	adox	h3, z4		// adox h3,z4	-> z4
	mulx	4*8(m), l4, h4	// mulx m[4]	-> h4:l4
	adcx	l4, z4		// adcx l4,z4	-> z4
	mov	z4, 3*8(z)	// store z4	-> z[3]

	mov	4*8(z), z5	// load z[4]	-> z5
	adox	h4, z5		// adox h4,z5	-> z5
	mulx	5*8(m), l5, h5	// mulx m[5]	-> h5:l5
	adcx	l5, z5		// adcx l5,z5	-> z5
	mov	z5, 4*8(z)	// store z5	-> z[4]

	mov	5*8(z), z6	// load z[5]	-> z6
	adox	h5, z6		// adox h5,z6	-> z6
	mulx	6*8(m), l6, h6	// mulx m[6]	-> h6:l6
	adcx	l6, z6		// adcx l6,z6	-> z6
	mov	z6, 5*8(z)	// store z6	-> z[5]

	mov	6*8(z), z7	// load z[6]	-> z7
	adox	h6, z7		// adox h6,z7	-> z7
	mulx	7*8(m), l7, h7	// mulx m[7]	-> h7:l7
	adcx	l7, z7		// adcx l7,z7	-> z7
	mov	z7, 6*8(z)	// store z7	-> z[6]

	mov	7*8(z), z8	// load z[7]	-> z8
	adox	h7, z8		// adox h7,z8	-> z8
	mulx	8*8(m), l8, h8	// mulx m[8]	-> h8:l8
	adcx	l8, z8		// adcx l8,z8	-> z8
	mov	z8, 7*8(z)	// store z8	-> z[7]

	mov	8*8(z), z9	// load z[8]	-> z9
	adox	h8, z9		// adox h8,z9	-> z9
	mulx	9*8(m), l9, h9	// mulx m[9]	-> h9:l9
	adcx	l9, z9		// adcx l9,z9	-> z9
	mov	z9, 8*8(z)	// store z9	-> z[8]

	mov	9*8(z), z10	// load z[9]	-> z10
	adox	h9, z10		// adox h9,z10	-> z10
	mulx	10*8(m), l10, h10	// mulx m[10]	-> h10:l10
	adcx	l10, z10	// adcx l10,z10	-> z10
	mov	z10, 9*8(z)	// store z10	-> z[9]

	mov	10*8(z), z11	// load z[10]	-> z11
	adox	h10, z11	// adox h10,z11	-> z11
	mulx	11*8(m), l11, h11	// mulx m[11]	-> h11:l11
	adcx	l11, z11	// adcx l11,z11	-> z11
	mov	z11, 10*8(z)	// store z11	-> z[10]

	mov	11*8(z), z12	// load z[11]	-> z12
	adox	h11, z12	// adox h11,z12	-> z12
	adcx	z13, z12	// adcx 0,z12	-> z12
	mov	z12, 11*8(z)	// store z12	-> z[11]

////////////////////////////////////////////////////////////////

	mov	1*8(x), %rdx	// load x[1]	-> rdx
	call step
	mov	2*8(x), %rdx	// load x[2]	-> rdx
	call step
	mov	3*8(x), %rdx	// load x[3]	-> rdx
	call step
	mov	4*8(x), %rdx	// load x[4]	-> rdx
	call step
	mov	5*8(x), %rdx	// load x[5]	-> rdx
	call step
	mov	6*8(x), %rdx	// load x[6]	-> rdx
	call step
	mov	7*8(x), %rdx	// load x[7]	-> rdx
	call step
	mov	8*8(x), %rdx	// load x[8]	-> rdx
	call step
	mov	9*8(x), %rdx	// load x[9]	-> rdx
	call step
	mov	10*8(x), %rdx	// load x[10]	-> rdx
	call step
	mov	11*8(x), %rdx	// load x[11]	-> rdx
	call step

	// Conditional subtraction of m

#ifdef _WIN64
	mov	0*8(z), %rax
	sub	0*8(m), %rax

	mov	1*8(z), %rcx
	sbb	1*8(m), %rcx
#else
	mov	0*8(z), %r8
	sub	0*8(m), %r8

	mov	1*8(z), %r9
	sbb	1*8(m), %r9
#endif
	mov	2*8(z), %r14
	sbb	2*8(m), %r14

	mov	3*8(z), %r13
	sbb	3*8(m), %r13

	mov	4*8(z), %rdi
	sbb	4*8(m), %rdi

	mov	5*8(z), %rsi
	sbb	5*8(m), %rsi

	mov	6*8(z), %r12
	sbb	6*8(m), %r12

	mov	7*8(z), %rbp
	sbb	7*8(m), %rbp

	mov	8*8(z), %rbx
	sbb	8*8(m), %rbx

	mov	9*8(z), %rdx
	sbb	9*8(m), %rdx

	mov	10*8(z), %r10
	sbb	10*8(m), %r10

	mov	11*8(z), %r11
	sbb	11*8(m), %r11

	sbb	$0, z13

#ifdef _WIN64
	cmovc	0*8(z), %rax
	cmovc	1*8(z), %rcx
#else
	cmovc	0*8(z), %r8
	cmovc	1*8(z), %r9
#endif
	cmovc	2*8(z), %r14
	cmovc	3*8(z), %r13
	cmovc	4*8(z), %rdi
	cmovc	5*8(z), %rsi
	cmovc	6*8(z), %r12
	cmovc	7*8(z), %rbp
	cmovc	8*8(z), %rbx
	cmovc	9*8(z), %rdx
	cmovc	10*8(z), %r10
	cmovc	11*8(z), %r11

#ifdef _WIN64
	mov	%rax, 0*8(z)
	mov	%rcx, 1*8(z)
#else
	mov	%r8, 0*8(z)
	mov	%r9, 1*8(z)
#endif
	pop	%r15
	mov	%r14, 2*8(z)
	pop	%r14
	mov	%r13, 3*8(z)
	pop	%r13
	mov	%rdi, 4*8(z)
#ifdef _WIN64
	pop	%rdi
#endif
	mov	%rsi, 5*8(z)
#ifdef _WIN64
	pop	%rsi
#endif
	mov	%r12, 6*8(z)
	pop	%r12
	mov	%rbp, 7*8(z)
	pop	%rbp
	mov	%rbx, 8*8(z)
	pop	%rbx

	mov	%rdx, 9*8(z)
	mov	%r10, 10*8(z)
	mov	%r11, 11*8(z)

	ret

////////////////////////////////////////////////////////////////

.p2align 6,,15
step:

	mulx	0*8(y), l0, h0	// rdx * y[0]	-> h0:l0
	mov	0*8(z), z0	// load z[0]	-> z0
	xor	t0, t0		// clear flags
	adcx	l0, z0		// adcx l0,z0	-> z0

	mov	1*8(z), z1	// load z[1]	-> z1
	mulx	1*8(y), l1, h1	// rdx * y[1]	-> h1:l1
	adox	h0, z1		// adox h0,z1	-> z1
	adcx	l1, z1		// adcx l1,z1	-> z1
	mov	z1, 0*8(z)	// store z1	-> z[0]

	mov	2*8(z), z2	// load z[2]	-> z2
	mulx	2*8(y), l2, h2	// rdx * y[2]	-> h2:l2
	adox	h1, z2		// adox h1,z2	-> z2
	adcx	l2, z2		// adcx l2,z2	-> z2
	mov	z2, 1*8(z)	// store z2	-> z[1]

	mov	3*8(z), z3	// load z[3]	-> z3
	mulx	3*8(y), l3, h3	// rdx * y[3]	-> h3:l3
	adox	h2, z3		// adox h2,z3	-> z3
	adcx	l3, z3		// adcx l3,z3	-> z3
	mov	z3, 2*8(z)	// store z3	-> z[2]

	mov	4*8(z), z4	// load z[4]	-> z4
	mulx	4*8(y), l4, h4	// rdx * y[4]	-> h4:l4
	adox	h3, z4		// adox h3,z4	-> z4
	adcx	l4, z4		// adcx l4,z4	-> z4
	mov	z4, 3*8(z)	// store z4	-> z[3]

	mov	5*8(z), z5	// load z[5]	-> z5
	mulx	5*8(y), l5, h5	// rdx * y[5]	-> h5:l5
	adox	h4, z5		// adox h4,z5	-> z5
	adcx	l5, z5		// adcx l5,z5	-> z5
	mov	z5, 4*8(z)	// store z5	-> z[4]

	mov	%rdx, l4
	mov	z0, %rdx	// mov z0	-> rdx
	mulx	12*8(m), r0, h4	// mulx inv	-> h4:r0 // h4 discarded
	mov	l4, %rdx	// load x[3]	-> rdx

	mov	6*8(z), z6	// load z[6]	-> z6
	mulx	6*8(y), l6, h6	// rdx * y[6]	-> h6:l6
	adox	h5, z6		// adox h5,z6	-> z6
	adcx	l6, z6		// adcx l6,z6	-> z6
	mov	z6, 5*8(z)	// store z6	-> z[5]

	mov	7*8(z), z7	// load z[7]	-> z7
	mulx	7*8(y), l7, h7	// rdx * y[7]	-> h7:l7
	adox	h6, z7		// adox h6,z7	-> z7
	adcx	l7, z7		// adcx l7,z7	-> z7
	mov	z7, 6*8(z)	// store z7	-> z[6]

	mov	8*8(z), z8	// load z[8]	-> z8
	mulx	8*8(y), l8, h8	// rdx * y[8]	-> h8:l8
	adox	h7, z8		// adox h7,z8	-> z8
	adcx	l8, z8		// adcx l8,z8	-> z8
	mov	z8, 7*8(z)	// store z8	-> z[7]

	mov	9*8(z), z9	// load z[9]	-> z9
	mulx	9*8(y), l9, h9	// rdx * y[9]	-> h9:l9
	adox	h8, z9		// adox h8,z9	-> z9
	adcx	l9, z9		// adcx l9,z9	-> z9
	mov	z9, 8*8(z)	// store z9	-> z[8]

	mov	10*8(z), z10	// load z[10]	-> z10
	mulx	10*8(y), l10, h10	// rdx * y[10]	-> h10:l10
	adox	h9, z10		// adox h9,z10	-> z10
	adcx	l10, z10	// adcx l10,z10	-> z10
	mov	z10, 9*8(z)	// store z10	-> z[9]

	mov	11*8(z), z11	// load z[11]	-> z11
	mulx	11*8(y), l11, z12	// rdx * y[11]	-> z12:l11
	adox	h10, z11	// adox h10,z11	-> z11
	adcx	l11, z11	// adcx l11,z11	-> z11
	mov	z11, 10*8(z)	// store z11	-> z[10]

	adox	z13, z12	// adox z13,z12	-> z12
	adcx	t0, z12		// adcx 0,z12	-> z12
	mov	t0, z13		// mov 0	-> z13
	mov	z12, 11*8(z)	// store z12	-> z[11]

	// Reduction

	mov	r0, %rdx	// r0		-> rdx
	xor	t0, t0		// clear flags

	mov	0*8(m), l0	// load m[0]	-> l0
	mulx	l0, l0, h0	// mulx l0	-> h0:l0
	adcx	l0, z0		// adcx l0,z0	-> z0	// 0, ignore

	mov	0*8(z), z1	// load z[0]	-> z1
	adox	h0, z1		// adox h0,z1	-> z1
	mulx	1*8(m), l1, h1	// mulx m[1]	-> h1:l1
	adcx	l1, z1		// adcx l1,z1	-> z1
	mov	z1, 0*8(z)	// store z1	-> z[0]

	mov	1*8(z), z2	// load z[1]	-> z2
	adox	h1, z2		// adox h1,z2	-> z2
	mulx	2*8(m), l2, h2	// mulx m[2]	-> h2:l2
	adcx	l2, z2		// adcx l2,z2	-> z2
	mov	z2, 1*8(z)	// store z2	-> z[1]

	mov	2*8(z), z3	// load z[2]	-> z3
	adox	h2, z3		// adox h2,z3	-> z3
	mulx	3*8(m), l3, h3	// mulx m[3]	-> h3:l3
	adcx	l3, z3		// adcx l3,z3	-> z3
	mov	z3, 2*8(z)	// store z3	-> z[2]

	mov	3*8(z), z4	// load z[3]	-> z4
	adox	h3, z4		// adox h3,z4	-> z4
	mulx	4*8(m), l4, h4	// mulx m[4]	-> h4:l4
	adcx	l4, z4		// adcx l4,z4	-> z4
	mov	z4, 3*8(z)	// store z4	-> z[3]

	mov	4*8(z), z5	// load z[4]	-> z5
	adox	h4, z5		// adox h4,z5	-> z5
	mulx	5*8(m), l5, h5	// mulx m[5]	-> h5:l5
	adcx	l5, z5		// adcx l5,z5	-> z5
	mov	z5, 4*8(z)	// store z5	-> z[4]

	mov	5*8(z), z6	// load z[5]	-> z6
	adox	h5, z6		// adox h5,z6	-> z6
	mulx	6*8(m), l6, h6	// mulx m[6]	-> h6:l6
	adcx	l6, z6		// adcx l6,z6	-> z6
	mov	z6, 5*8(z)	// store z6	-> z[5]

	mov	6*8(z), z7	// load z[6]	-> z7
	adox	h6, z7		// adox h6,z7	-> z7
	mulx	7*8(m), l7, h7	// mulx m[7]	-> h7:l7
	adcx	l7, z7		// adcx l7,z7	-> z7
	mov	z7, 6*8(z)	// store z7	-> z[6]

	mov	7*8(z), z8	// load z[7]	-> z8
	adox	h7, z8		// adox h7,z8	-> z8
	mulx	8*8(m), l8, h8	// mulx m[8]	-> h8:l8
	adcx	l8, z8		// adcx l8,z8	-> z8
	mov	z8, 7*8(z)	// store z8	-> z[7]

	mov	8*8(z), z9	// load z[8]	-> z9
	adox	h8, z9		// adox h8,z9	-> z9
	mulx	9*8(m), l9, h9	// mulx m[9]	-> h9:l9
	adcx	l9, z9		// adcx l9,z9	-> z9
	mov	z9, 8*8(z)	// store z9	-> z[8]

	mov	9*8(z), z10	// load z[9]	-> z10
	adox	h9, z10		// adox h9,z10	-> z10
	mulx	10*8(m), l10, h10	// mulx m[10]	-> h10:l10
	adcx	l10, z10	// adcx l10,z10	-> z10
	mov	z10, 9*8(z)	// store z10	-> z[9]

	mov	10*8(z), z11	// load z[10]	-> z11
	adox	h10, z11	// adox h10,z11	-> z11
	mulx	11*8(m), l11, h11	// mulx m[11]	-> h11:l11
	adcx	l11, z11	// adcx l11,z11	-> z11
	mov	z11, 10*8(z)	// store z11	-> z[10]

	mov	11*8(z), z12	// load z[11]	-> z12
	adox	h11, z12	// adox h11,z12	-> z12
	adcx	t0, z12		// adcx 0,z12	-> z12
	mov	z12, 11*8(z)	// store z12	-> z[11]

	ret

#ifndef __APPLE__
.size   modmul768, .-modmul768
#endif

// vim: noet ts=8 sw=8
