use crate::{LinearCombination, Variable, Vec};
use algebra_core::Field;
use core::ops::{Add, AddAssign, Deref, DerefMut, Mul, MulAssign, Neg, Sub};

impl<F: Field> LinearCombination<F> {
    /// Create a new empty linear combination.
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// Create a new empty linear combination.
    pub fn zero() -> Self {
        Self(Vec::new())
    }

    /// Deduplicate entries in `self`.
    pub fn compactify(&mut self) {
        self.0.sort_by_key(|e| e.1);
        let mut current_var = None;
        let mut current_var_first_index = 0;
        for i in 0..self.0.len() {
            let (f, v) = self.0[i];
            if Some(v) == current_var {
                self.0[current_var_first_index].0 += &f;
            } else {
                current_var = Some(v);
                current_var_first_index = i;
            }
        }
        self.0.dedup_by_key(|e| e.1);
    }
}

impl<'a, F: Field> Deref for LinearCombination<F> {
    type Target = Vec<(F, Variable)>;

    #[inline]
    fn deref(&self) -> &Vec<(F, Variable)> {
        &self.0
    }
}

impl<F: Field> DerefMut for LinearCombination<F> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<F: Field> From<(F, Variable)> for LinearCombination<F> {
    #[inline]
    fn from(input: (F, Variable)) -> Self {
        LinearCombination(crate::vec![input])
    }
}

impl<F: Field> From<Variable> for LinearCombination<F> {
    #[inline]
    fn from(var: Variable) -> Self {
        LinearCombination(crate::vec![(F::one(), var)])
    }
}

impl<F: Field> LinearCombination<F> {
    /// Negate the coefficients of all variables in `self`.
    #[inline]
    pub fn negate_in_place(&mut self) {
        self.0.iter_mut().for_each(|(coeff, _)| *coeff = -(*coeff));
    }

    /// Get the location of a variable in `self`.
    #[inline]
    pub fn get_var_loc(&self, search_var: &Variable) -> Result<usize, usize> {
        if self.0.len() < 6 {
            let mut found_index = 0;
            for (i, (_, var)) in self.iter().enumerate() {
                if var >= search_var {
                    found_index = i;
                    break;
                } else {
                    found_index += 1;
                }
            }
            Err(found_index)
        } else {
            self.0
                .binary_search_by_key(search_var, |&(_, cur_var)| cur_var)
        }
    }
}

impl<F: Field> Add<(F, Variable)> for LinearCombination<F> {
    type Output = Self;

    #[inline]
    fn add(mut self, coeff_var: (F, Variable)) -> Self {
        self += coeff_var;
        self
    }
}

impl<F: Field> AddAssign<(F, Variable)> for LinearCombination<F> {
    #[inline]
    fn add_assign(&mut self, (coeff, var): (F, Variable)) {
        match self.get_var_loc(&var) {
            Ok(found) => self.0[found].0 += &coeff,
            Err(not_found) => self.0.insert(not_found, (coeff, var)),
        }
    }
}

impl<F: Field> Sub<(F, Variable)> for LinearCombination<F> {
    type Output = Self;

    #[inline]
    fn sub(self, (coeff, var): (F, Variable)) -> Self {
        self + (-coeff, var)
    }
}

impl<F: Field> Neg for LinearCombination<F> {
    type Output = Self;

    #[inline]
    fn neg(mut self) -> Self {
        self.negate_in_place();
        self
    }
}

impl<F: Field> Mul<F> for LinearCombination<F> {
    type Output = Self;

    #[inline]
    fn mul(mut self, scalar: F) -> Self {
        self *= scalar;
        self
    }
}

impl<'a, F: Field> Mul<F> for &'a LinearCombination<F> {
    type Output = LinearCombination<F>;

    #[inline]
    fn mul(self, scalar: F) -> LinearCombination<F> {
        let mut cur = self.clone();
        cur *= scalar;
        cur
    }
}

impl<F: Field> MulAssign<F> for LinearCombination<F> {
    #[inline]
    fn mul_assign(&mut self, scalar: F) {
        self.0.iter_mut().for_each(|(coeff, _)| *coeff *= &scalar);
    }
}

impl<F: Field> Add<Variable> for LinearCombination<F> {
    type Output = Self;

    #[inline]
    fn add(self, other: Variable) -> LinearCombination<F> {
        self + (F::one(), other)
    }
}

impl<'a, F: Field> Add<&'a Variable> for LinearCombination<F> {
    type Output = Self;

    #[inline]
    fn add(self, other: &'a Variable) -> LinearCombination<F> {
        self + *other
    }
}

impl<'a, F: Field> Sub<&'a Variable> for LinearCombination<F> {
    type Output = Self;

    #[inline]
    fn sub(self, other: &'a Variable) -> LinearCombination<F> {
        self - *other
    }
}

impl<F: Field> Sub<Variable> for LinearCombination<F> {
    type Output = LinearCombination<F>;

    #[inline]
    fn sub(self, other: Variable) -> LinearCombination<F> {
        self - (F::one(), other)
    }
}

fn op_impl<F: Field, F1, F2>(
    cur: &LinearCombination<F>,
    other: &LinearCombination<F>,
    push_fn: F1,
    combine_fn: F2,
) -> LinearCombination<F>
where
    F1: Fn(F) -> F,
    F2: Fn(F, F) -> F,
{
    let mut new_vec = Vec::new();
    let mut i = 0;
    let mut j = 0;
    while i < cur.len() && j < other.len() {
        let self_cur = &cur[i];
        let other_cur = &other[j];
        if self_cur.1 > other_cur.1 {
            new_vec.push((push_fn(other[j].0), other[j].1));
            j += 1;
        } else if self_cur.1 < other_cur.1 {
            new_vec.push(*self_cur);
            i += 1;
        } else {
            new_vec.push((combine_fn(self_cur.0, other_cur.0), self_cur.1));
            i += 1;
            j += 1;
        }
    }
    new_vec.extend_from_slice(&cur[i..]);
    while j < other.0.len() {
        new_vec.push((push_fn(other[j].0), other[j].1));
        j += 1;
    }
    LinearCombination(new_vec)
}

impl<F: Field> Add<&LinearCombination<F>> for &LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn add(self, other: &LinearCombination<F>) -> LinearCombination<F> {
        if other.0.is_empty() {
            return self.clone();
        } else if self.0.is_empty() {
            return other.clone();
        }
        op_impl(
            self,
            other,
            |coeff| coeff,
            |cur_coeff, other_coeff| cur_coeff + &other_coeff,
        )
    }
}

impl<F: Field> Add<LinearCombination<F>> for &LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn add(self, other: LinearCombination<F>) -> LinearCombination<F> {
        if self.0.is_empty() {
            return other;
        } else if other.0.is_empty() {
            return self.clone();
        }
        op_impl(
            self,
            &other,
            |coeff| coeff,
            |cur_coeff, other_coeff| cur_coeff + &other_coeff,
        )
    }
}

impl<'a, F: Field> Add<&'a LinearCombination<F>> for LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn add(self, other: &'a LinearCombination<F>) -> LinearCombination<F> {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            return other.clone();
        }
        op_impl(
            &self,
            other,
            |coeff| coeff,
            |cur_coeff, other_coeff| cur_coeff + &other_coeff,
        )
    }
}

impl<F: Field> Add<LinearCombination<F>> for LinearCombination<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            return other;
        }
        op_impl(
            &self,
            &other,
            |coeff| coeff,
            |cur_coeff, other_coeff| cur_coeff + &other_coeff,
        )
    }
}

impl<F: Field> Sub<&LinearCombination<F>> for &LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn sub(self, other: &LinearCombination<F>) -> LinearCombination<F> {
        if other.0.is_empty() {
            let cur = self.clone();
            return cur;
        } else if self.0.is_empty() {
            let mut other = other.clone();
            other.negate_in_place();
            return other;
        }

        op_impl(
            self,
            other,
            |coeff| -coeff,
            |cur_coeff, other_coeff| cur_coeff - &other_coeff,
        )
    }
}

impl<'a, F: Field> Sub<&'a LinearCombination<F>> for LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn sub(self, other: &'a LinearCombination<F>) -> LinearCombination<F> {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            let mut other = other.clone();
            other.negate_in_place();
            return other;
        }
        op_impl(
            &self,
            other,
            |coeff| -coeff,
            |cur_coeff, other_coeff| cur_coeff - &other_coeff,
        )
    }
}

impl<F: Field> Sub<LinearCombination<F>> for &LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn sub(self, mut other: LinearCombination<F>) -> LinearCombination<F> {
        if self.0.is_empty() {
            other.negate_in_place();
            return other;
        } else if other.0.is_empty() {
            return self.clone();
        }

        op_impl(
            self,
            &other,
            |coeff| -coeff,
            |cur_coeff, other_coeff| cur_coeff - &other_coeff,
        )
    }
}

impl<F: Field> Sub<LinearCombination<F>> for LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn sub(self, mut other: LinearCombination<F>) -> LinearCombination<F> {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            other.negate_in_place();
            return other;
        }
        op_impl(
            &self,
            &other,
            |coeff| -coeff,
            |cur_coeff, other_coeff| cur_coeff - &other_coeff,
        )
    }
}

impl<F: Field> Add<(F, &LinearCombination<F>)> for &LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn add(self, (mul_coeff, other): (F, &LinearCombination<F>)) -> LinearCombination<F> {
        if other.0.is_empty() {
            return self.clone();
        } else if self.0.is_empty() {
            let mut other = other.clone();
            other.mul_assign(mul_coeff);
            return other;
        }
        op_impl(
            self,
            other,
            |coeff| mul_coeff * &coeff,
            |cur_coeff, other_coeff| cur_coeff + &(mul_coeff * &other_coeff),
        )
    }
}

impl<'a, F: Field> Add<(F, &'a LinearCombination<F>)> for LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn add(self, (mul_coeff, other): (F, &'a LinearCombination<F>)) -> LinearCombination<F> {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            let mut other = other.clone();
            other.mul_assign(mul_coeff);
            return other;
        }
        op_impl(
            &self,
            other,
            |coeff| mul_coeff * &coeff,
            |cur_coeff, other_coeff| cur_coeff + &(mul_coeff * &other_coeff),
        )
    }
}

impl<F: Field> Add<(F, LinearCombination<F>)> for &LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn add(self, (mul_coeff, mut other): (F, LinearCombination<F>)) -> LinearCombination<F> {
        if other.0.is_empty() {
            return self.clone();
        } else if self.0.is_empty() {
            other.mul_assign(mul_coeff);
            return other;
        }
        op_impl(
            self,
            &other,
            |coeff| mul_coeff * &coeff,
            |cur_coeff, other_coeff| cur_coeff + &(mul_coeff * &other_coeff),
        )
    }
}

impl<F: Field> Add<(F, Self)> for LinearCombination<F> {
    type Output = Self;

    fn add(self, (mul_coeff, other): (F, Self)) -> Self {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            let mut other = other;
            other.mul_assign(mul_coeff);
            return other;
        }
        op_impl(
            &self,
            &other,
            |coeff| mul_coeff * &coeff,
            |cur_coeff, other_coeff| cur_coeff + &(mul_coeff * &other_coeff),
        )
    }
}

impl<F: Field> Sub<(F, &LinearCombination<F>)> for &LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn sub(self, (coeff, other): (F, &LinearCombination<F>)) -> LinearCombination<F> {
        self + (-coeff, other)
    }
}

impl<'a, F: Field> Sub<(F, &'a LinearCombination<F>)> for LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn sub(self, (coeff, other): (F, &'a LinearCombination<F>)) -> LinearCombination<F> {
        self + (-coeff, other)
    }
}

impl<F: Field> Sub<(F, LinearCombination<F>)> for &LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn sub(self, (coeff, other): (F, LinearCombination<F>)) -> LinearCombination<F> {
        self + (-coeff, other)
    }
}

impl<'a, F: Field> Sub<(F, LinearCombination<F>)> for LinearCombination<F> {
    type Output = LinearCombination<F>;

    fn sub(self, (coeff, other): (F, LinearCombination<F>)) -> LinearCombination<F> {
        self + (-coeff, other)
    }
}
