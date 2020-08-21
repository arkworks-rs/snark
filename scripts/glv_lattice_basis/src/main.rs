use algebra_core::bigint::BigInteger;
// We work on arrays of size 3
// We assume that |E(F_q)| < R = 2^{ceil(limbs/2) * 64}
fn extended_euclidean<BigInt: BigInteger>(n: BigInt, lambda: BigInt) -> ((BigInt, BigInt), (BigInt, BigInt)) {
    let mut r = [n, lambda, n];
    let one = iBigInteger::<BigInt>{ value: BigInt::from(1), neg: false };
    let zero = iBigInteger::<BigInt>{ value: BigInt::from(0), neg: false };
    let mut s = [one, zero, zero];
    let mut t = [zero, one, zero];

    let sqrt_n = as_f64(n.0).sqrt();

    let mut i = 0;
    // While r_i >= sqrt(n), we then return the vectors (r_i, t_i), (r_i+1, t_i+1)
    while as_f64(r[(i + 1) % 3].0) >= sqrt_n {
        let (q, r): (BigInt, BigInt) = div_with_remainder::<BigInt>(r[i % 3], r[(i + 1) % 3]);
        r[(i + 2) % 3] = r;
        let int_q = iBigInteger::<BigInt>::from(q);
        s[(i + 2) % 3] = s[i % 3] - int_q * (s[(i + 1) % 3]);
        t[(i + 2) % 3] = t[i % 3] - int_q * (t[(i + 1) % 3]);

    }
    i += 1;

    vec_1 = (r[(i + 1) % 3], t[(i + 2) % 3].value)
}

fn as_f64(bigint_ref: &[u64]) -> f64 {
    let mut n_float: f64 = 0;
    for (i, limb) in n.iter().enumerate() {
        n_float += (limb as f64) * 2.pow((i as i32) * 64i32)
    }
    n_float
}

struct iBigInteger<BigInt: BigInteger> {
    value: BigInt,
    neg: bool,
}

impl iBigInteger {}

impl<BigInt: BigInteger> Mul for iBigInteger<BigInt> {
    fn mul_assign(&mut self, other: &Self) {
        self.value *= other.value;
        match (self.neg, other.neg) {
            (true, true) => self.neg(),
            (false, true) => self.neg(),
            _ => (),
        }
    }
}

impl<BigInt: BigInteger> Neg for iBigInteger<BigInt> {
    fn neg(&mut self) {
        if self.neg {
            self.neg = false;
        } else {
            self.neg = true;
        }
    }
}

impl<BigInt: BigInteger> Sub for iBigInteger<BigInt> {
    fn sub_assign(&mut self, other: &Self) {
        self.add_nocarry(other.neg());
    }
}

impl<BigInt: BigInteger> Add for iBigInteger<BigInt> {
    fn add_assign(&mut self, other: &Self) {
        // If operators have the same sign, just add the values
        if self.neg + other.neg == false {
            self.value += other.value;
        } else {
            if self.value > other.value {
                self.sub_noborrow(other);
            } else {
                let mut tmp = other.clone();
                tmp.sub_noborrow(self.value);
                self.value = tmp;
                self.neg();
            }
        }
    }
}

impl<BigInt: BigInteger> From<BigInt> for iBigInteger<BigInt> {
    #[inline]
    fn from(val: BigInt) -> iBigInteger<BigInt> {
        iBigInteger::<BigInt>{
            value: val,
            neg: false,
        }
    }
}
