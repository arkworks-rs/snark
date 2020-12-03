// void modsqr768(const uint64_t a[12], const uint64_t m[13], uint64_t c[12])

// m[12] contains the least significant word of the negated inverse of the modulus mod 2^768

// Function parameters

#define a	x0
#define m	x1
#define c	x2

#define cl0	x4
#define cl1	x5
#define cl2	x6
#define cl3	x7
#define cl4	x8
#define cl5	x9
#define cl6	x10
#define cl7	x11
#define cl8	x12
#define cl9	x13
#define cl10	x14
#define cl11	x15

#define cm0	x16
#define cm1	x17
#define cm2	x19
#define cm3	x20
#define cm4	x21
#define cm5	x22
#define cm6	x23
#define cm7	x24
#define cm8	x25
#define cm9	x26
#define cm10	x27
#define cm11	x28

#define ch0	cl0
#define ch1	cl1
#define ch2	cl2
#define ch3	cl3
#define ch4	cl4
#define ch5	cl5
#define ch6	cl6
#define ch7	cl7
#define ch8	cl8
#define ch9	cl9
#define ch10	cl10
#define ch11	cl11

#define a0	cm0
#define a1	cm1
#define a2	cm2
#define a3	cm3
#define a4	cm4
#define a5	cm5

#define b0	cm6
#define b1	cm7
#define b2	cm8
#define b3	cm9
#define b4	cm10
#define b5	cm11

#define t0	cl0
#define t1	cl1
#define t2	cl2
#define t3	cl3
#define t4	cl4
#define t5	cl5

#define u0	cl6
#define u1	cl7
#define u2	cl8
#define u3	cl9
#define u4	cl10
#define u5	cl11

#define pp	c
#define sa	m
#define sb	c
#define sc	x3

#define inv	a
#define r	x3

#define m0	cm0
#define m1	cm1
#define m2	cm2
#define m3	cm3
#define m4	cm4
#define m5	cm5
#define m6	cm6
#define m7	cm7
#define m8	cm8
#define m9	cm9
#define m10	cm10
#define m11	cm11

#define cl12	cl0
#define cl13	cl1
#define cl14	cl2
#define cl15	cl3
#define cl16	cl4
#define cl17	cl5
#define cl18	cl6
#define cl19	cl7
#define cl20	cl8
#define cl21	cl9
#define cl22	cl10
#define cl23	cl11
#define cl24	inv

// 384-bit squaring

// Off-center diagonals: 3+5+7+9 additions
// Doubling: 11 additions
// Center diagonal: 11 additions
// Total: 46 additions

#define	sqr384(a, c) \
\
	/* Off-center diagonals */\
\
	mul	c##5, a##5, a##0;\
	umulh	c##6, a##5, a##0;\
\
	mul	c##4, a##4, a##0;\
	umulh	pp, a##4, a##0;	adds	 c##5, c##5, pp;\
	mul	pp, a##5, a##1;	adcs	 c##6, c##6, pp;\
	umulh	pp, a##5, a##1;	adc	 c##7,  XZR, pp;\
\
	mul	c##3, a##3, a##0;\
	umulh	pp, a##3, a##0;	adds	 c##4, c##4, pp;\
	mul	pp, a##4, a##1;	adcs	 c##5, c##5, pp;\
	umulh	pp, a##4, a##1;	adcs	 c##6, c##6, pp;\
	mul	pp, a##5, a##2;	adcs	 c##7, c##7, pp;\
	umulh	pp, a##5, a##2;	adc	 c##8,  XZR, pp;\
\
	mul	c##2, a##2, a##0;\
	umulh	pp, a##2, a##0;	adds	 c##3, c##3, pp;\
	mul	pp, a##3, a##1;	adcs	 c##4, c##4, pp;\
	umulh	pp, a##3, a##1;	adcs	 c##5, c##5, pp;\
	mul	pp, a##4, a##2;	adcs	 c##6, c##6, pp;\
	umulh	pp, a##4, a##2;	adcs	 c##7, c##7, pp;\
	mul	pp, a##5, a##3;	adcs	 c##8, c##8, pp;\
	umulh	pp, a##5, a##3;	adc	 c##9,  XZR, pp;\
\
	mul	c##1, a##1, a##0;\
	umulh	pp, a##1, a##0;	adds	 c##2, c##2, pp;\
	mul	pp, a##2, a##1;	adcs	 c##3, c##3, pp;\
	umulh	pp, a##2, a##1;	adcs	 c##4, c##4, pp;\
	mul	pp, a##3, a##2;	adcs	 c##5, c##5, pp;\
	umulh	pp, a##3, a##2;	adcs	 c##6, c##6, pp;\
	mul	pp, a##4, a##3;	adcs	 c##7, c##7, pp;\
	umulh	pp, a##4, a##3;	adcs	 c##8, c##8, pp;\
	mul	pp, a##5, a##4;	adcs	 c##9, c##9, pp;\
	umulh	pp, a##5, a##4;	adc	c##10,  XZR, pp;\
\
	/* Double */\
\
	adds	 c##1,  c##1,  c##1;\
	adcs	 c##2,  c##2,  c##2;\
	adcs	 c##3,  c##3,  c##3;\
	adcs	 c##4,  c##4,  c##4;\
	adcs	 c##5,  c##5,  c##5;\
	adcs	 c##6,  c##6,  c##6;\
	adcs	 c##7,  c##7,  c##7;\
	adcs	 c##8,  c##8,  c##8;\
	adcs	 c##9,  c##9,  c##9;\
	adcs	c##10, c##10, c##10;\
	adc	c##11, XZR, XZR;\
\
	/* Main diagonal */\
\
	mul	c##0, a##0, a##0;\
	umulh	pp, a##0, a##0;	adds	 c##1,  c##1, pp;\
	mul	pp, a##1, a##1;	adcs	 c##2,  c##2, pp;\
	umulh	pp, a##1, a##1;	adcs	 c##3,  c##3, pp;\
	mul	pp, a##2, a##2;	adcs	 c##4,  c##4, pp;\
	umulh	pp, a##2, a##2;	adcs	 c##5,  c##5, pp;\
	mul	pp, a##3, a##3;	adcs	 c##6,  c##6, pp;\
	umulh	pp, a##3, a##3;	adcs	 c##7,  c##7, pp;\
	mul	pp, a##4, a##4;	adcs	 c##8,  c##8, pp;\
	umulh	pp, a##4, a##4;	adcs	 c##9,  c##9, pp;\
	mul	pp, a##5, a##5;	adcs	c##10, c##10, pp;\
	umulh	pp, a##5, a##5;	adc	c##11, c##11, pp;

#ifdef __APPLE__
#define modsqr768 _modsqr768
#endif

.global modsqr768

.text
.p2align 6,,63

modsqr768:
	// Allocate space for storing 36 registers

	sub	SP, SP, #36*8

	// Preserve m, c, and callee-save registers r19-r28    

	stp	x19, x20, [SP, # 0*8]
	stp	x21, x22, [SP, # 2*8]
	stp	x23, x24, [SP, # 4*8]
	stp	x25, x26, [SP, # 6*8]
	stp	x27, x28, [SP, # 8*8]
	stp	  m,   c, [SP, #10*8]

	//// Karatsuba Squaring ////

	// Load a

	ldp	a0, a1, [a, # 0*8]
	ldp	a2, a3, [a, # 2*8]
	ldp	a4, a5, [a, # 4*8]
	ldp	b0, b1, [a, # 6*8]
	ldp	b2, b3, [a, # 8*8]
	ldp	b4, b5, [a, #10*8]

	// Compute al^2

	sqr384(a, cl)

	// Save cl to stack

	stp	 cl0,  cl1, [SP, #12*8]
	stp	 cl2,  cl3, [SP, #14*8]
	stp	 cl4,  cl5, [SP, #16*8]
	stp	 cl6,  cl7, [SP, #18*8]
	stp	 cl8,  cl9, [SP, #20*8]
	stp	cl10, cl11, [SP, #22*8]

	// Compute ah^2

	sqr384(b, ch)

	// Save ch to stack

	stp	 ch0,  ch1, [SP, #24*8]
	stp	 ch2,  ch3, [SP, #26*8]
	stp	 ch4,  ch5, [SP, #28*8]
	stp	 ch6,  ch7, [SP, #30*8]
	stp	 ch8,  ch9, [SP, #32*8]
	stp	ch10, ch11, [SP, #34*8]

	// Compute t = |al-ah|

	subs	t0, a0, b0
	sbcs	t1, a1, b1
	sbcs	t2, a2, b2
	sbcs	t3, a3, b3
	sbcs	t4, a4, b4
	sbcs	t5, a5, b5
	sbc	sa, XZR, XZR	// -1 iff t<0, else 0

	eor	t0, t0, sa
	eor	t1, t1, sa
	eor	t2, t2, sa
	eor	t3, t3, sa
	eor	t4, t4, sa
	eor	t5, t5, sa

	subs	t0, t0, sa
	sbcs	t1, t1, sa
	sbcs	t2, t2, sa
	sbcs	t3, t3, sa
	sbcs	t4, t4, sa
	sbc	t5, t5, sa

	// Compute (al-ah)^2

	sqr384(t, cm)

	// Load cl

	ldp	 cl0,  cl1, [SP, #12*8]
	ldp	 cl2,  cl3, [SP, #14*8]
	ldp	 cl4,  cl5, [SP, #16*8]
	ldp	 cl6,  cl7, [SP, #18*8]
	ldp	 cl8,  cl9, [SP, #20*8]
	ldp	cl10, cl11, [SP, #22*8]

	// cm = cl - cm

	subs	 cm0,  cl0,  cm0
	sbcs	 cm1,  cl1,  cm1
	sbcs	 cm2,  cl2,  cm2
	sbcs	 cm3,  cl3,  cm3
	sbcs	 cm4,  cl4,  cm4
	sbcs	 cm5,  cl5,  cm5
	sbcs	 cm6,  cl6,  cm6
	sbcs	 cm7,  cl7,  cm7
	sbcs	 cm8,  cl8,  cm8
	sbcs	 cm9,  cl9,  cm9
	sbcs	cm10, cl10, cm10
	sbcs	cm11, cl11, cm11
	sbc	  sc,  XZR,  XZR

	// Load ch from stack

	ldp	 ch0,  ch1, [SP, #24*8]
	ldp	 ch2,  ch3, [SP, #26*8]
	ldp	 ch4,  ch5, [SP, #28*8]
	ldp	 ch6,  ch7, [SP, #30*8]
	ldp	 ch8,  ch9, [SP, #32*8]
	ldp	ch10, ch11, [SP, #34*8]

	// cm = cm + ch

	adds	 cm0,  cm0,  ch0
	adcs	 cm1,  cm1,  ch1
	adcs	 cm2,  cm2,  ch2
	adcs	 cm3,  cm3,  ch3
	adcs	 cm4,  cm4,  ch4
	adcs	 cm5,  cm5,  ch5
	adcs	 cm6,  cm6,  ch6
	adcs	 cm7,  cm7,  ch7
	adcs	 cm8,  cm8,  ch8
	adcs	 cm9,  cm9,  ch9
	adcs	cm10, cm10, ch10
	adcs	cm11, cm11, ch11
	adc	  sc,   sc,  XZR

	// Compute high tree quarters of a^2

	ldp	sa, sb, [SP, #18*8];	adds	cm0, cm0, sa;	adcs	cm1, cm1, sb
	ldp	sa, sb, [SP, #20*8];	adcs	cm2, cm2, sa;	adcs	cm3, cm3, sb
	ldp	sa, sb, [SP, #22*8];	adcs	cm4, cm4, sa;	adcs	cm5, cm5, sb

	adcs	 ch0,  ch0, cm6
	adcs	 ch1,  ch1, cm7
	adcs	 ch2,  ch2, cm8
	adcs	 ch3,  ch3, cm9
	adcs	 ch4,  ch4, cm10
	adcs	 ch5,  ch5, cm11
	adcs	 ch6,  ch6,  sc
	adcs	 ch7,  ch7, XZR
	adcs	 ch8,  ch8, XZR
	adcs	 ch9,  ch9, XZR
	adcs	ch10, ch10, XZR
	adc	ch11, ch11, XZR

	// Save ch to stack

	stp	 ch0,  ch1, [SP, #24*8]
	stp	 ch2,  ch3, [SP, #26*8]
	stp	 ch4,  ch5, [SP, #28*8]
	stp	 ch6,  ch7, [SP, #30*8]
	stp	 ch8,  ch9, [SP, #32*8]
	stp	ch10, ch11, [SP, #34*8]

	// Load low half of cl

	ldp	 cl0,  cl1, [SP, #12*8]
	ldp	 cl2,  cl3, [SP, #14*8]
	ldp	 cl4,  cl5, [SP, #16*8]

	// Move bottom half of cm to top half of cl

	mov	 cl6, cm0
	mov	 cl7, cm1
	mov	 cl8, cm2
	mov	 cl9, cm3
	mov	cl10, cm4
	mov	cl11, cm5

	//// Reduction (Operand Scanning) ////

	// Restore m

	ldr	m, [SP, #10*8]

	// Load m[]

	ldp	 m0,  m1, [m, # 0*8]
	ldp	 m2,  m3, [m, # 2*8]
	ldp	 m4,  m5, [m, # 4*8]
	ldp	 m6,  m7, [m, # 6*8]
	ldp	 m8,  m9, [m, # 8*8]
	ldp	m10, m11, [m, #10*8]
	ldr	inv, [m, #12*8]

#undef pp
#define pp m

	mul	r, inv, cl0

	mul	pp,  m0, r;	adds	 cl0,  cl0, pp
	mul	pp,  m1, r;	adcs	 cl1,  cl1, pp
	mul	pp,  m2, r;	adcs	 cl2,  cl2, pp
	mul	pp,  m3, r;	adcs	 cl3,  cl3, pp
	mul	pp,  m4, r;	adcs	 cl4,  cl4, pp
	mul	pp,  m5, r;	adcs	 cl5,  cl5, pp
	mul	pp,  m6, r;	adcs	 cl6,  cl6, pp
	mul	pp,  m7, r;	adcs	 cl7,  cl7, pp
	mul	pp,  m8, r;	adcs	 cl8,  cl8, pp
	mul	pp,  m9, r;	adcs	 cl9,  cl9, pp
	mul	pp, m10, r;	adcs	cl10, cl10, pp
	mul	pp, m11, r;	adcs	cl11, cl11, pp
				adc	cl12, XZR, XZR

	umulh	pp,  m0, r;	adds	 cl1,  cl1, pp
	umulh	pp,  m1, r;	adcs	 cl2,  cl2, pp
	umulh	pp,  m2, r;	adcs	 cl3,  cl3, pp
	umulh	pp,  m3, r;	adcs	 cl4,  cl4, pp
	umulh	pp,  m4, r;	adcs	 cl5,  cl5, pp
	umulh	pp,  m5, r;	adcs	 cl6,  cl6, pp
	umulh	pp,  m6, r;	adcs	 cl7,  cl7, pp
	umulh	pp,  m7, r;	adcs	 cl8,  cl8, pp
	umulh	pp,  m8, r;	adcs	 cl9,  cl9, pp
	umulh	pp,  m9, r;	adcs	cl10, cl10, pp
	umulh	pp, m10, r;	adcs	cl11, cl11, pp
	umulh	pp, m11, r;	adc	cl12, cl12, pp

	mul	r, inv, cl1

	mul	pp,  m0, r;	adds	 cl1,  cl1, pp
	mul	pp,  m1, r;	adcs	 cl2,  cl2, pp
	mul	pp,  m2, r;	adcs	 cl3,  cl3, pp
	mul	pp,  m3, r;	adcs	 cl4,  cl4, pp
	mul	pp,  m4, r;	adcs	 cl5,  cl5, pp
	mul	pp,  m5, r;	adcs	 cl6,  cl6, pp
	mul	pp,  m6, r;	adcs	 cl7,  cl7, pp
	mul	pp,  m7, r;	adcs	 cl8,  cl8, pp
	mul	pp,  m8, r;	adcs	 cl9,  cl9, pp
	mul	pp,  m9, r;	adcs	cl10, cl10, pp
	mul	pp, m10, r;	adcs	cl11, cl11, pp
	mul	pp, m11, r;	adcs	cl12, cl12, pp
				adc	cl13, XZR, XZR

	umulh	pp,  m0, r;	adds	 cl2,  cl2, pp
	umulh	pp,  m1, r;	adcs	 cl3,  cl3, pp
	umulh	pp,  m2, r;	adcs	 cl4,  cl4, pp
	umulh	pp,  m3, r;	adcs	 cl5,  cl5, pp
	umulh	pp,  m4, r;	adcs	 cl6,  cl6, pp
	umulh	pp,  m5, r;	adcs	 cl7,  cl7, pp
	umulh	pp,  m6, r;	adcs	 cl8,  cl8, pp
	umulh	pp,  m7, r;	adcs	 cl9,  cl9, pp
	umulh	pp,  m8, r;	adcs	cl10, cl10, pp
	umulh	pp,  m9, r;	adcs	cl11, cl11, pp
	umulh	pp, m10, r;	adcs	cl12, cl12, pp
	umulh	pp, m11, r;	adc	cl13, cl13, pp

	mul	r, inv, cl2

	mul	pp,  m0, r;	adds	 cl2,  cl2, pp
	mul	pp,  m1, r;	adcs	 cl3,  cl3, pp
	mul	pp,  m2, r;	adcs	 cl4,  cl4, pp
	mul	pp,  m3, r;	adcs	 cl5,  cl5, pp
	mul	pp,  m4, r;	adcs	 cl6,  cl6, pp
	mul	pp,  m5, r;	adcs	 cl7,  cl7, pp
	mul	pp,  m6, r;	adcs	 cl8,  cl8, pp
	mul	pp,  m7, r;	adcs	 cl9,  cl9, pp
	mul	pp,  m8, r;	adcs	cl10, cl10, pp
	mul	pp,  m9, r;	adcs	cl11, cl11, pp
	mul	pp, m10, r;	adcs	cl12, cl12, pp
	mul	pp, m11, r;	adcs	cl13, cl13, pp
				adc	cl14, XZR, XZR

	umulh	pp,  m0, r;	adds	 cl3,  cl3, pp
	umulh	pp,  m1, r;	adcs	 cl4,  cl4, pp
	umulh	pp,  m2, r;	adcs	 cl5,  cl5, pp
	umulh	pp,  m3, r;	adcs	 cl6,  cl6, pp
	umulh	pp,  m4, r;	adcs	 cl7,  cl7, pp
	umulh	pp,  m5, r;	adcs	 cl8,  cl8, pp
	umulh	pp,  m6, r;	adcs	 cl9,  cl9, pp
	umulh	pp,  m7, r;	adcs	cl10, cl10, pp
	umulh	pp,  m8, r;	adcs	cl11, cl11, pp
	umulh	pp,  m9, r;	adcs	cl12, cl12, pp
	umulh	pp, m10, r;	adcs	cl13, cl13, pp
	umulh	pp, m11, r;	adc	cl14, cl14, pp

	mul	r, inv, cl3

	mul	pp,  m0, r;	adds	 cl3,  cl3, pp
	mul	pp,  m1, r;	adcs	 cl4,  cl4, pp
	mul	pp,  m2, r;	adcs	 cl5,  cl5, pp
	mul	pp,  m3, r;	adcs	 cl6,  cl6, pp
	mul	pp,  m4, r;	adcs	 cl7,  cl7, pp
	mul	pp,  m5, r;	adcs	 cl8,  cl8, pp
	mul	pp,  m6, r;	adcs	 cl9,  cl9, pp
	mul	pp,  m7, r;	adcs	cl10, cl10, pp
	mul	pp,  m8, r;	adcs	cl11, cl11, pp
	mul	pp,  m9, r;	adcs	cl12, cl12, pp
	mul	pp, m10, r;	adcs	cl13, cl13, pp
	mul	pp, m11, r;	adcs	cl14, cl14, pp
				adc	cl15, XZR, XZR

	umulh	pp,  m0, r;	adds	 cl4,  cl4, pp
	umulh	pp,  m1, r;	adcs	 cl5,  cl5, pp
	umulh	pp,  m2, r;	adcs	 cl6,  cl6, pp
	umulh	pp,  m3, r;	adcs	 cl7,  cl7, pp
	umulh	pp,  m4, r;	adcs	 cl8,  cl8, pp
	umulh	pp,  m5, r;	adcs	 cl9,  cl9, pp
	umulh	pp,  m6, r;	adcs	cl10, cl10, pp
	umulh	pp,  m7, r;	adcs	cl11, cl11, pp
	umulh	pp,  m8, r;	adcs	cl12, cl12, pp
	umulh	pp,  m9, r;	adcs	cl13, cl13, pp
	umulh	pp, m10, r;	adcs	cl14, cl14, pp
	umulh	pp, m11, r;	adc	cl15, cl15, pp

	mul	r, inv, cl4

	mul	pp,  m0, r;	adds	 cl4,  cl4, pp
	mul	pp,  m1, r;	adcs	 cl5,  cl5, pp
	mul	pp,  m2, r;	adcs	 cl6,  cl6, pp
	mul	pp,  m3, r;	adcs	 cl7,  cl7, pp
	mul	pp,  m4, r;	adcs	 cl8,  cl8, pp
	mul	pp,  m5, r;	adcs	 cl9,  cl9, pp
	mul	pp,  m6, r;	adcs	cl10, cl10, pp
	mul	pp,  m7, r;	adcs	cl11, cl11, pp
	mul	pp,  m8, r;	adcs	cl12, cl12, pp
	mul	pp,  m9, r;	adcs	cl13, cl13, pp
	mul	pp, m10, r;	adcs	cl14, cl14, pp
	mul	pp, m11, r;	adcs	cl15, cl15, pp
				adc	cl16, XZR, XZR

	umulh	pp,  m0, r;	adds	 cl5,  cl5, pp
	umulh	pp,  m1, r;	adcs	 cl6,  cl6, pp
	umulh	pp,  m2, r;	adcs	 cl7,  cl7, pp
	umulh	pp,  m3, r;	adcs	 cl8,  cl8, pp
	umulh	pp,  m4, r;	adcs	 cl9,  cl9, pp
	umulh	pp,  m5, r;	adcs	cl10, cl10, pp
	umulh	pp,  m6, r;	adcs	cl11, cl11, pp
	umulh	pp,  m7, r;	adcs	cl12, cl12, pp
	umulh	pp,  m8, r;	adcs	cl13, cl13, pp
	umulh	pp,  m9, r;	adcs	cl14, cl14, pp
	umulh	pp, m10, r;	adcs	cl15, cl15, pp
	umulh	pp, m11, r;	adc	cl16, cl16, pp

	mul	r, inv, cl5

	mul	pp,  m0, r;	adds	 cl5,  cl5, pp
	mul	pp,  m1, r;	adcs	 cl6,  cl6, pp
	mul	pp,  m2, r;	adcs	 cl7,  cl7, pp
	mul	pp,  m3, r;	adcs	 cl8,  cl8, pp
	mul	pp,  m4, r;	adcs	 cl9,  cl9, pp
	mul	pp,  m5, r;	adcs	cl10, cl10, pp
	mul	pp,  m6, r;	adcs	cl11, cl11, pp
	mul	pp,  m7, r;	adcs	cl12, cl12, pp
	mul	pp,  m8, r;	adcs	cl13, cl13, pp
	mul	pp,  m9, r;	adcs	cl14, cl14, pp
	mul	pp, m10, r;	adcs	cl15, cl15, pp
	mul	pp, m11, r;	adcs	cl16, cl16, pp
				adc	cl17, XZR, XZR

	umulh	pp,  m0, r;	adds	 cl6,  cl6, pp
	umulh	pp,  m1, r;	adcs	 cl7,  cl7, pp
	umulh	pp,  m2, r;	adcs	 cl8,  cl8, pp
	umulh	pp,  m3, r;	adcs	 cl9,  cl9, pp
	umulh	pp,  m4, r;	adcs	cl10, cl10, pp
	umulh	pp,  m5, r;	adcs	cl11, cl11, pp
	umulh	pp,  m6, r;	adcs	cl12, cl12, pp
	umulh	pp,  m7, r;	adcs	cl13, cl13, pp
	umulh	pp,  m8, r;	adcs	cl14, cl14, pp
	umulh	pp,  m9, r;	adcs	cl15, cl15, pp
	umulh	pp, m10, r;	adcs	cl16, cl16, pp
	umulh	pp, m11, r;	adc	cl17, cl17, pp

	mul	r, inv, cl6

	mul	pp,  m0, r;	adds	 cl6,  cl6, pp
	mul	pp,  m1, r;	adcs	 cl7,  cl7, pp
	mul	pp,  m2, r;	adcs	 cl8,  cl8, pp
	mul	pp,  m3, r;	adcs	 cl9,  cl9, pp
	mul	pp,  m4, r;	adcs	cl10, cl10, pp
	mul	pp,  m5, r;	adcs	cl11, cl11, pp
	mul	pp,  m6, r;	adcs	cl12, cl12, pp
	mul	pp,  m7, r;	adcs	cl13, cl13, pp
	mul	pp,  m8, r;	adcs	cl14, cl14, pp
	mul	pp,  m9, r;	adcs	cl15, cl15, pp
	mul	pp, m10, r;	adcs	cl16, cl16, pp
	mul	pp, m11, r;	adcs	cl17, cl17, pp
				adc	cl18, XZR, XZR

	umulh	pp,  m0, r;	adds	 cl7,  cl7, pp
	umulh	pp,  m1, r;	adcs	 cl8,  cl8, pp
	umulh	pp,  m2, r;	adcs	 cl9,  cl9, pp
	umulh	pp,  m3, r;	adcs	cl10, cl10, pp
	umulh	pp,  m4, r;	adcs	cl11, cl11, pp
	umulh	pp,  m5, r;	adcs	cl12, cl12, pp
	umulh	pp,  m6, r;	adcs	cl13, cl13, pp
	umulh	pp,  m7, r;	adcs	cl14, cl14, pp
	umulh	pp,  m8, r;	adcs	cl15, cl15, pp
	umulh	pp,  m9, r;	adcs	cl16, cl16, pp
	umulh	pp, m10, r;	adcs	cl17, cl17, pp
	umulh	pp, m11, r;	adc	cl18, cl18, pp

	mul	r, inv, cl7

	mul	pp,  m0, r;	adds	 cl7,  cl7, pp
	mul	pp,  m1, r;	adcs	 cl8,  cl8, pp
	mul	pp,  m2, r;	adcs	 cl9,  cl9, pp
	mul	pp,  m3, r;	adcs	cl10, cl10, pp
	mul	pp,  m4, r;	adcs	cl11, cl11, pp
	mul	pp,  m5, r;	adcs	cl12, cl12, pp
	mul	pp,  m6, r;	adcs	cl13, cl13, pp
	mul	pp,  m7, r;	adcs	cl14, cl14, pp
	mul	pp,  m8, r;	adcs	cl15, cl15, pp
	mul	pp,  m9, r;	adcs	cl16, cl16, pp
	mul	pp, m10, r;	adcs	cl17, cl17, pp
	mul	pp, m11, r;	adcs	cl18, cl18, pp
				adc	cl19, XZR, XZR

	umulh	pp,  m0, r;	adds	 cl8,  cl8, pp
	umulh	pp,  m1, r;	adcs	 cl9,  cl9, pp
	umulh	pp,  m2, r;	adcs	cl10, cl10, pp
	umulh	pp,  m3, r;	adcs	cl11, cl11, pp
	umulh	pp,  m4, r;	adcs	cl12, cl12, pp
	umulh	pp,  m5, r;	adcs	cl13, cl13, pp
	umulh	pp,  m6, r;	adcs	cl14, cl14, pp
	umulh	pp,  m7, r;	adcs	cl15, cl15, pp
	umulh	pp,  m8, r;	adcs	cl16, cl16, pp
	umulh	pp,  m9, r;	adcs	cl17, cl17, pp
	umulh	pp, m10, r;	adcs	cl18, cl18, pp
	umulh	pp, m11, r;	adc	cl19, cl19, pp

	mul	r, inv, cl8

	mul	pp,  m0, r;	adds	 cl8,  cl8, pp
	mul	pp,  m1, r;	adcs	 cl9,  cl9, pp
	mul	pp,  m2, r;	adcs	cl10, cl10, pp
	mul	pp,  m3, r;	adcs	cl11, cl11, pp
	mul	pp,  m4, r;	adcs	cl12, cl12, pp
	mul	pp,  m5, r;	adcs	cl13, cl13, pp
	mul	pp,  m6, r;	adcs	cl14, cl14, pp
	mul	pp,  m7, r;	adcs	cl15, cl15, pp
	mul	pp,  m8, r;	adcs	cl16, cl16, pp
	mul	pp,  m9, r;	adcs	cl17, cl17, pp
	mul	pp, m10, r;	adcs	cl18, cl18, pp
	mul	pp, m11, r;	adcs	cl19, cl19, pp
				adc	cl20, XZR, XZR

	umulh	pp,  m0, r;	adds	 cl9,  cl9, pp
	umulh	pp,  m1, r;	adcs	cl10, cl10, pp
	umulh	pp,  m2, r;	adcs	cl11, cl11, pp
	umulh	pp,  m3, r;	adcs	cl12, cl12, pp
	umulh	pp,  m4, r;	adcs	cl13, cl13, pp
	umulh	pp,  m5, r;	adcs	cl14, cl14, pp
	umulh	pp,  m6, r;	adcs	cl15, cl15, pp
	umulh	pp,  m7, r;	adcs	cl16, cl16, pp
	umulh	pp,  m8, r;	adcs	cl17, cl17, pp
	umulh	pp,  m9, r;	adcs	cl18, cl18, pp
	umulh	pp, m10, r;	adcs	cl19, cl19, pp
	umulh	pp, m11, r;	adc	cl20, cl20, pp

	mul	r, inv, cl9

	mul	pp,  m0, r;	adds	 cl9,  cl9, pp
	mul	pp,  m1, r;	adcs	cl10, cl10, pp
	mul	pp,  m2, r;	adcs	cl11, cl11, pp
	mul	pp,  m3, r;	adcs	cl12, cl12, pp
	mul	pp,  m4, r;	adcs	cl13, cl13, pp
	mul	pp,  m5, r;	adcs	cl14, cl14, pp
	mul	pp,  m6, r;	adcs	cl15, cl15, pp
	mul	pp,  m7, r;	adcs	cl16, cl16, pp
	mul	pp,  m8, r;	adcs	cl17, cl17, pp
	mul	pp,  m9, r;	adcs	cl18, cl18, pp
	mul	pp, m10, r;	adcs	cl19, cl19, pp
	mul	pp, m11, r;	adcs	cl20, cl20, pp
				adc	cl21, XZR, XZR

	umulh	pp,  m0, r;	adds	cl10, cl10, pp
	umulh	pp,  m1, r;	adcs	cl11, cl11, pp
	umulh	pp,  m2, r;	adcs	cl12, cl12, pp
	umulh	pp,  m3, r;	adcs	cl13, cl13, pp
	umulh	pp,  m4, r;	adcs	cl14, cl14, pp
	umulh	pp,  m5, r;	adcs	cl15, cl15, pp
	umulh	pp,  m6, r;	adcs	cl16, cl16, pp
	umulh	pp,  m7, r;	adcs	cl17, cl17, pp
	umulh	pp,  m8, r;	adcs	cl18, cl18, pp
	umulh	pp,  m9, r;	adcs	cl19, cl19, pp
	umulh	pp, m10, r;	adcs	cl20, cl20, pp
	umulh	pp, m11, r;	adc	cl21, cl21, pp

	mul	r, inv, cl10

	mul	pp,  m0, r;	adds	cl10, cl10, pp
	mul	pp,  m1, r;	adcs	cl11, cl11, pp
	mul	pp,  m2, r;	adcs	cl12, cl12, pp
	mul	pp,  m3, r;	adcs	cl13, cl13, pp
	mul	pp,  m4, r;	adcs	cl14, cl14, pp
	mul	pp,  m5, r;	adcs	cl15, cl15, pp
	mul	pp,  m6, r;	adcs	cl16, cl16, pp
	mul	pp,  m7, r;	adcs	cl17, cl17, pp
	mul	pp,  m8, r;	adcs	cl18, cl18, pp
	mul	pp,  m9, r;	adcs	cl19, cl19, pp
	mul	pp, m10, r;	adcs	cl20, cl20, pp
	mul	pp, m11, r;	adcs	cl21, cl21, pp
				adc	cl22, XZR, XZR

	umulh	pp,  m0, r;	adds	cl11, cl11, pp
	umulh	pp,  m1, r;	adcs	cl12, cl12, pp
	umulh	pp,  m2, r;	adcs	cl13, cl13, pp
	umulh	pp,  m3, r;	adcs	cl14, cl14, pp
	umulh	pp,  m4, r;	adcs	cl15, cl15, pp
	umulh	pp,  m5, r;	adcs	cl16, cl16, pp
	umulh	pp,  m6, r;	adcs	cl17, cl17, pp
	umulh	pp,  m7, r;	adcs	cl18, cl18, pp
	umulh	pp,  m8, r;	adcs	cl19, cl19, pp
	umulh	pp,  m9, r;	adcs	cl20, cl20, pp
	umulh	pp, m10, r;	adcs	cl21, cl21, pp
	umulh	pp, m11, r;	adc	cl22, cl22, pp

	mul	r, inv, cl11

	mul	pp,  m0, r;	adds	cl11, cl11, pp
	mul	pp,  m1, r;	adcs	cl12, cl12, pp
	mul	pp,  m2, r;	adcs	cl13, cl13, pp
	mul	pp,  m3, r;	adcs	cl14, cl14, pp
	mul	pp,  m4, r;	adcs	cl15, cl15, pp
	mul	pp,  m5, r;	adcs	cl16, cl16, pp
	mul	pp,  m6, r;	adcs	cl17, cl17, pp
	mul	pp,  m7, r;	adcs	cl18, cl18, pp
	mul	pp,  m8, r;	adcs	cl19, cl19, pp
	mul	pp,  m9, r;	adcs	cl20, cl20, pp
	mul	pp, m10, r;	adcs	cl21, cl21, pp
	mul	pp, m11, r;	adcs	cl22, cl22, pp
				adc	cl23, XZR, XZR

	umulh	pp,  m0, r;	adds	cl12, cl12, pp
	umulh	pp,  m1, r;	adcs	cl13, cl13, pp
	umulh	pp,  m2, r;	adcs	cl14, cl14, pp
	umulh	pp,  m3, r;	adcs	cl15, cl15, pp
	umulh	pp,  m4, r;	adcs	cl16, cl16, pp
	umulh	pp,  m5, r;	adcs	cl17, cl17, pp
	umulh	pp,  m6, r;	adcs	cl18, cl18, pp
	umulh	pp,  m7, r;	adcs	cl19, cl19, pp
	umulh	pp,  m8, r;	adcs	cl20, cl20, pp
	umulh	pp,  m9, r;	adcs	cl21, cl21, pp
	umulh	pp, m10, r;	adcs	cl22, cl22, pp
	umulh	pp, m11, r;	adc	cl23, cl23, pp

	// Add top half of product

	ldp	a, m, [SP, #24*8];	adds	cl12, cl12, a;	adcs	cl13, cl13, m
	ldp	a, m, [SP, #26*8];	adcs	cl14, cl14, a;	adcs	cl15, cl15, m
	ldp	a, m, [SP, #28*8];	adcs	cl16, cl16, a;	adcs	cl17, cl17, m
	ldp	a, m, [SP, #30*8];	adcs	cl18, cl18, a;	adcs	cl19, cl19, m
	ldp	a, m, [SP, #32*8];	adcs	cl20, cl20, a;	adcs	cl21, cl21, m
	ldp	a, m, [SP, #34*8];	adcs	cl22, cl22, a;	adcs	cl23, cl23, m

	adc	cl24, XZR, XZR

	// Conditionally subtract m:
	//
	// m=cl-m
	// if (m >= 0)
	//   cl=m
	// store cl

	subs	 m0, cl12,  m0
	sbcs	 m1, cl13,  m1
	sbcs	 m2, cl14,  m2
	sbcs	 m3, cl15,  m3
	sbcs	 m4, cl16,  m4
	sbcs	 m5, cl17,  m5
	sbcs	 m6, cl18,  m6
	sbcs	 m7, cl19,  m7
	sbcs	 m8, cl20,  m8
	sbcs	 m9, cl21,  m9
	sbcs	m10, cl22, m10
	sbcs	m11, cl23, m11
	sbcs	XZR, cl24, XZR

	csel	cl12, cl12,  m0, cc
	csel	cl13, cl13,  m1, cc
	csel	cl14, cl14,  m2, cc
	csel	cl15, cl15,  m3, cc
	csel	cl16, cl16,  m4, cc
	csel	cl17, cl17,  m5, cc
	csel	cl18, cl18,  m6, cc
	csel	cl19, cl19,  m7, cc
	csel	cl20, cl20,  m8, cc
	csel	cl21, cl21,  m9, cc
	csel	cl22, cl22, m10, cc
	csel	cl23, cl23, m11, cc

	// Restore c

	ldr	c, [SP, #11*8]

	// Save result

	stp	cl12, cl13, [c, # 0*8]
	stp	cl14, cl15, [c, # 2*8]
	stp	cl16, cl17, [c, # 4*8]
	stp	cl18, cl19, [c, # 6*8]
	stp	cl20, cl21, [c, # 8*8]
	stp	cl22, cl23, [c, #10*8]

	//// Cleanup ////

	// Restore callee-save registers r19-r28

	ldp	x19, x20, [SP, #0*8]
	ldp	x21, x22, [SP, #2*8]
	ldp	x23, x24, [SP, #4*8]
	ldp	x25, x26, [SP, #6*8]
	ldp	x27, x28, [SP, #8*8]

	// Restore stack pointer

	add	SP, SP, #36*8

	ret
