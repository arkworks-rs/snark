#![allow(unsafe_code)]

use ark_ff::Field;
use ark_std::vec::*;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::gr1cs::{
    field_interner::{FieldInterner, InternedField},
    Variable,
};

/// The following invariant must always be maintained:
/// 1. `self.offsets` is a non-empty vector of length `n + 1`, where
///    `n` is the number of linear combinations.
/// 2. `self.offsets[i]` is the index in `self.coeffs` and `self.vars`
///    where the `i`-th linear combination starts.
/// 3. `self.offsets[i + 1]` is the index in `self.coeffs` and `self.vars`
///    where the `i`-th linear combination ends.
/// 4. `self.coeffs` and `self.vars` are of the same length,
///    and they contain the coefficients and variables of the linear combinations respectively.
/// 5. `self.coeffs` and `self.vars` are interleaved such that for each linear combination `i`,
///    the coefficients and variables are stored in the same order, i.e.,
///    letting `start = self.offsets[i]` and `end = self.offsets[i + 1]`, then
///    `self.vars[start..end]` corresponds to `self.coeffs[start..end]
///
/// Invariants 2 and 3 imply the following lemma:
/// Lemma 1. `self.offsets[i] <= self.offsets[i + 1]` for all `i` in `0..n`.
/// Proof:
///   By invariant 2, `self.offsets[i]` is the index in `self.coeffs` and `self.vars`
///   where the `i`-th linear combination starts.
///   By invariant 3, `self.offsets[i + 1]` is the index in `self.coeffs` and `self.vars`
///   where the `i`-th linear combination ends.
///   Since each linear combination is of length at least 0, we have the lemma.
///
/// Invariants 2, 3, 4, 5 imply the following lemma:
/// Lemma 2. For all `i` in `0..n`,
///   * `self.offsets[i + 1] - self.offsets[i] <= self.vars.len()`
///   * `self.offsets[i + 1] - self.offsets[i] <= self.coeffs.len()`
///   
/// Proof:
///   Assume not. Then there exists an `i` such that either
///   * `self.offsets[i + 1] - self.offsets[i] > self.vars.len()`
///   * `self.offsets[i + 1] - self.offsets[i] > self.coeffs.len()`
///
///   In both cases, this would imply that the `i`-th linear combination
///   has more variables or coefficients than the total number of variables or coefficients,
///   which contradicts invariant 5.
#[derive(Debug, Clone, Default)]
pub struct LcMap<F: Field> {
    vars: Vec<Variable>,
    coeffs: Vec<InternedField>,
    offsets: Vec<usize>,
    _f: core::marker::PhantomData<F>,
}

pub(crate) fn to_non_interned_lc<'a, F: Field>(
    lc: impl Iterator<Item = (&'a InternedField, &'a Variable)> + 'a,
    f_interner: &'a FieldInterner<F>,
) -> impl Iterator<Item = (F, Variable)> + 'a {
    lc.map(|(&c, &v)| (f_interner.value(c).unwrap(), v))
}

type LcMapIterItem<'a> =
    core::iter::Zip<core::slice::Iter<'a, InternedField>, core::slice::Iter<'a, Variable>>;
type LcVarsIterMutItem<'a> = core::slice::IterMut<'a, Variable>;

impl<F: Field> LcMap<F> {
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            vars: Vec::new(),
            coeffs: Vec::new(),
            // We preserve invariant 1., i.e. that `self.offsets` has length `n + 1` where `n` is the number of linear combinations.
            // Initially, `n = 0`, so `self.offsets` has length 1.
            offsets: vec![0],
            _f: core::marker::PhantomData,
        }
    }

    #[inline(always)]
    pub fn with_capacity(expected_num_lcs: usize, expected_total_lc_size: usize) -> Self {
        let mut result = Self::new();
        result.vars.reserve(expected_total_lc_size);
        result.coeffs.reserve(expected_total_lc_size);
        result.offsets.reserve(expected_num_lcs + 1);
        result
    }

    #[inline(always)]
    pub fn push(
        &mut self,
        v: impl IntoIterator<Item = (F, Variable)>,
        f_interner: &mut FieldInterner<F>,
    ) {
        // This loop preserves invariants 4 and 5, i.e. that
        // * `self.coeffs.len() == self.vars.len()`
        // * `self.coeffs` and `self.vars` are interleaved such that for each linear combination `i`, the coefficients and variables are stored in the same order.
        for (coeff, var) in v {
            self.coeffs.push(f_interner.get_or_intern(coeff));
            self.vars.push(var);
        }
        // This step preserves invariants 1, 2, and 3, i.e. that
        // * `self.offsets` has length `n + 1` where `n` is the number of linear combinations.
        // * `self.offsets[i]` is the index in `self.coeffs` and `self.vars` where the `i`-th linear combination starts.
        // * `self.offsets[i + 1]` is the index in `self.coeffs` and `self.vars` where the `i`-th linear combination ends.
        self.offsets.push(self.coeffs.len());
    }

    #[inline(always)]
    pub fn push_by_ref<'a>(
        &mut self,
        v: impl IntoIterator<Item = &'a (F, Variable)>,
        f_interner: &mut FieldInterner<F>,
    ) {
        // See `push` for why the invariants are preserved.
        for (coeff, var) in v {
            self.coeffs.push(f_interner.get_or_intern(*coeff));
            self.vars.push(*var);
        }
        self.offsets.push(self.coeffs.len());
    }

    #[inline(always)]
    pub fn iter(&self) -> impl Iterator<Item = LcMapIterItem<'_>> {
        self.offsets.windows(2).map(|w|
                // SAFETY:
                // Precondition 1 (`w.len() == 2`) is satisfied for the following reason:
                //   By invariant 1, `self.offsets` is of length `n + 1`.
                //   Clearly, when `n == 0`, `self.offsets.next()` will return `None`, 
                //.  so this closure will never be called.
                //   Next, when `n > 0`, `self.offsets.next()` will return slices of length 2 due to
                //   the definition of windows(2) .
                //
                // Precondition 2 (`w[0] <= w[1]`) is satisfied due to Lemma 1.
                // Precondition 3 (`w[1] <= coeffs.len() && w[1] <= vars.len()`) is satisfied due to Invariant 3 and Invariant 5.
                unsafe { windowed_access(w, &self.coeffs, &self.vars) })
    }

    #[cfg(feature = "parallel")]
    #[inline(always)]
    pub fn par_iter(&self) -> impl ParallelIterator<Item = LcMapIterItem<'_>> {
        self.offsets.par_windows(2).map(|w|
                // SAFETY:
                // Precondition 1 (`w.len() == 2`) is satisfied for the following reason:
                //   By invariant 1, `self.offsets` is of length `n + 1`.
                //   Clearly, when `n == 0`, `self.offsets.next()` will return `None`, 
                //.  so this closure will never be called.
                //   Next, when `n > 0`, `self.offsets.next()` will return slices of length 2 due to
                //   the definition of windows(2) .
                //
                // Precondition 2 (`w[0] <= w[1]`) is satisfied due to Lemma 1.
                // Precondition 3 (`w[1] <= coeffs.len() && w[1] <= vars.len()`) is satisfied due to Invariant 3 and Invariant 5.
                unsafe { windowed_access(w, &self.coeffs, &self.vars) })
    }

    #[inline(always)]
    pub fn lc_vars_iter_mut(&mut self) -> impl Iterator<Item = LcVarsIterMutItem<'_>> {
        LcVarsIterMut::new(self)
    }

    #[cfg(feature = "parallel")]
    #[inline(always)]
    pub fn lc_vars_par_iter_mut(&mut self) -> LcVarsParIterMut<'_> {
        LcVarsParIterMut::new(self)
    }

    #[inline(always)]
    pub fn num_lcs(&self) -> usize {
        self.offsets.len() - 1
    }

    #[inline(always)]
    pub fn total_lc_size(&self) -> usize {
        self.vars.len()
    }

    #[allow(unsafe_code)]
    #[inline(always)]
    pub fn get(&self, idx: usize) -> Option<LcMapIterItem<'_>> {
        if idx >= self.num_lcs() || self.offsets.len() < 2 {
            cold()
        } else {
            unsafe {
                // SAFETY:
                // Precondition 1 (`self.offsets.len() >= 2`) is satisfied due to the check
                // in the `if` condition.
                //
                // Precondition 2 (`offsets[idx] <= offsets[idx + 1]`) is satisfied due to Lemma 1.
                // Precondition 3 (`offsets[idx + 1] <= coeffs.len() && offsets[idx + 1] <= vars.len()`)
                // is satisfied due to Invariant 3 and Invariant 5.
                Some(windowed_access(
                    // SAFETY:
                    // `idx < self.len()` implies that `idx < self.offsets.len() - 1`,
                    // and so `self.offsets.get_unchecked(idx..=(idx + 1))` is a valid slice of length 2.
                    self.offsets.get_unchecked(idx..=(idx + 1)),
                    &self.coeffs,
                    &self.vars,
                ))
            }
        }
    }
}

/// Preconditions:
/// 1. `w` is a slice of length 2.
/// 2. `w[0] <= w[1]`.
/// 3. `w[1]` is within bounds of `coeffs` and `vars`,
///    i.e. `w[1] <= self.coeffs.len()` and `w[1] <= self.vars.len()`.
#[inline(always)]
unsafe fn windowed_access<'a>(
    w: &'a [usize],
    coeffs: &'a [InternedField],
    vars: &'a [Variable],
) -> LcMapIterItem<'a> {
    debug_assert!(w.len() == 2, "Expected a slice of length 2");
    debug_assert!(w[0] <= w[1], "Expected w[0] <= w[1]");
    debug_assert!(w[1] <= coeffs.len(), "Expected w[1] <= coeffs.len()");
    unsafe {
        // SAFETY:
        // By precondition 1, `w` is a slice of length 2, so these accesses are within bounds.
        let start = *w.get_unchecked(0);
        let end = *w.get_unchecked(1);
        // By precondition 2, `start..end` is not an empty range.
        // By precondition 3, `end` is within bounds of `coeffs` and `vars`.
        coeffs
            .get_unchecked(start..end)
            .iter()
            .zip(vars.get_unchecked(start..end))
    }
}

/// Preconditions:
/// 1. `w` is a slice of length 2.
/// 2. `w[0] <= w[1]`.
/// 3. `w[1] - w[0] < vars.len()`.
///
/// Note that precondition 3 here differs from that in `windowed_access`;
/// instead of requiring that `w[1] < vars.len()`,
/// we require that `w[1] - w[0]` is within bounds of `vars`.
///
/// Note that `windowed_access_mut` mutates `vars` by splitting it into two parts:
/// the first part is the linear combination variables corresponding to the range `w[0]..w[1]`,
/// and the second part is the remaining variables.
/// `vars` is updated to point to the second part.
#[inline(always)]
unsafe fn windowed_access_mut<'b>(
    w: &[usize],
    vars: &mut &'b mut [Variable],
) -> LcVarsIterMutItem<'b> {
    debug_assert!(w.len() == 2, "Expected a slice of length 2");
    debug_assert!(w[0] <= w[1], "Expected w[0] <= w[1]");
    debug_assert!(
        w[1] - w[0] < vars.len(),
        "Expected w[1] - w[0] < vars.len()"
    );
    #[allow(unsafe_code)]
    unsafe {
        let start = *w.get_unchecked(0);
        let end = *w.get_unchecked(1);
        let len = end - start;
        let (v_head, v_tail) = core::mem::take(vars).split_at_mut_unchecked(len);
        *vars = v_tail;
        v_head.iter_mut()
    }
}

#[cold]
fn cold<'a>() -> Option<LcMapIterItem<'a>> {
    None
}

pub struct LcVarsIterMut<'a> {
    vars: &'a mut [Variable],
    offsets: core::slice::Windows<'a, usize>,
}

impl<'a> LcVarsIterMut<'a> {
    #[inline(always)]
    pub fn new<F: Field>(lc_map: &'a mut LcMap<F>) -> Self {
        Self {
            vars: &mut lc_map.vars,
            offsets: lc_map.offsets.windows(2),
        }
    }
}

impl<'a> Iterator for LcVarsIterMut<'a> {
    type Item = LcVarsIterMutItem<'a>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.offsets.next().map(move |w| {
            // SAFETY:
            // Precondition 1 (`w.len() == 2`) is satisfied for the following reason:
            //   By invariant 1, `self.offsets` is of length `n + 1`.
            //   Clearly, when `n == 0`, `self.offsets.next()` will return `None`,
            //.  so this closure will never be called.
            //   Next, when `n > 0`, `self.offsets.next()` will return slices of length 2 due to
            //   the definition of windows(2) .
            //
            // Precondition 2 (`w[0] <= w[1]`) is satisfied due to Lemma 1.
            // Precondition 3 (`w[1] - w[0] < vars.len()`) is satisfied due to Lemma 2.
            unsafe { windowed_access_mut(w, &mut self.vars) }
        })
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.offsets.len() - 1;
        (len, Some(len))
    }
}

#[cfg(feature = "parallel")]
pub struct LcVarsParIterMut<'a> {
    vars: &'a mut [Variable],
    offsets: &'a [usize],
}

#[cfg(feature = "parallel")]
impl<'a> LcVarsParIterMut<'a> {
    #[inline(always)]
    pub fn new<F: Field>(lc_map: &'a mut LcMap<F>) -> Self {
        Self {
            vars: &mut lc_map.vars,
            offsets: &lc_map.offsets,
        }
    }
}

#[cfg(feature = "parallel")]
impl<'a> ParallelIterator for LcVarsParIterMut<'a> {
    type Item = LcVarsIterMutItem<'a>;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: rayon::iter::plumbing::UnindexedConsumer<Self::Item>,
    {
        rayon::iter::plumbing::bridge(self, consumer)
    }
}

#[cfg(feature = "parallel")]
impl<'a> IndexedParallelIterator for LcVarsParIterMut<'a> {
    fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    fn drive<C>(self, consumer: C) -> C::Result
    where
        C: rayon::iter::plumbing::Consumer<Self::Item>,
    {
        rayon::iter::plumbing::bridge(self, consumer)
    }

    fn with_producer<CB>(self, callback: CB) -> CB::Output
    where
        CB: rayon::iter::plumbing::ProducerCallback<Self::Item>,
    {
        // At all times, the contents of `Producer` must preserve similar invariants
        // as those in `LcMap`:
        // 1. `self.offsets` is a non-empty vector of length `n + 1`, where `n` is the number of linear combinations being iterated over.
        // 2. `self.offsets[i]` is the index in `self.vars` where the `i`-th linear combination starts.
        // 3. `self.offsets[i + 1]` is the index in `self.vars` where the `i`-th linear combination ends.
        // 4. `self.vars` is a mutable slice of variables that contains the variables of the linear combinations being iterated over.
        //
        // (Note that invariant 5 does not apply here, as we are not interleaving coefficients and variables in this iterator.)
        struct Producer<'a> {
            vars: &'a mut [Variable],
            offsets: &'a [usize],
        }

        impl<'a> rayon::iter::plumbing::Producer for Producer<'a> {
            type Item = LcVarsIterMutItem<'a>;
            type IntoIter = std::vec::IntoIter<LcVarsIterMutItem<'a>>;

            fn into_iter(mut self) -> Self::IntoIter {
                // SAFETY:
                // Precondition 1 (`w.len() == 2`) is satisfied for the following reason:
                //   By invariant 1, `self.offsets` is of length `n + 1`.
                //   Clearly, when `n == 0`, `self.offsets.next()` will return `None`,
                //.  so this closure will never be called.
                //   Next, when `n > 0`, `self.offsets.next()` will return slices of length 2 due to
                //   the definition of windows(2) .
                //
                // Precondition 2 (`w[0] <= w[1]`) is satisfied due to (an analog of) Lemma 1.
                // Precondition 3 (`w[1] - w[0] < vars.len()`) is satisfied due to (an analog of) Lemma 2.
                self.offsets
                    .windows(2)
                    .map(|w| unsafe { windowed_access_mut(w, &mut self.vars) })
                    .collect::<Vec<_>>()
                    .into_iter()
            }

            // The freshly produced `left` and `right` `Producer`s
            // preserve the `Producer` invariants.
            //
            //
            // First, the definition of `split_at` guarantees that
            // `index < n = self.offsets.len() - 1`.
            //
            // We then establish some common facts we will use in the reasoning below.
            //
            // Fact 1. By invariant 2 on `self`, `split_point` is the index in `self.vars`
            // where the `index`-th linear combination starts.
            //
            // The latter implies the following:
            //
            // Fact 2. `left_vars` will contain the variables of all linear combinations
            // from `0` to `index - 1`, (i.e. the first `index` linear combinations).
            // Fact 3. `right_vars` will contain the variables of all linear combinations
            // from `index` to `n - 1`, (i.e. the last `n - index` linear combinations).
            //
            // ## Invariant 1:
            //
            // If `self` contains `n` linear combinations, then
            // * `left_offsets` will contain `index + 1` elements,
            //   thus representing `index`-many linear combinations,
            // * `right_offsets` will contain `self.offsets.len() - index + 1` elements,
            //   thus representing `n - index`-many linear combinations.
            //
            //
            // ## Invariant 2:
            //
            // Fact 2 and Fact 3 imply that `left_vars` and `right_vars`
            // contain the correct variables.
            //
            // We are left to show that `left_offsets` and `right_offsets` contain
            // the correct starting indices of the linear combinations in `left_vars` and `right_vars`.
            // This is guaranteed by invariant 2 on `self`:
            // * `self.offsets[..=index]` contains the starting indices of the first `index` linear combinations,
            //   which are the linear combinations in `left_vars`.
            // * `self.offsets[index..]` contains the starting indices of the last `n - index` linear combinations,
            //   which are the linear combinations in `right_vars`.
            //
            // ## Invariant 3:
            //
            // Follows similarly to the reasoning above, but just replacing starting indices with ending indices.
            //
            // ## Invariant 4:
            //
            // Follows from Fact 2 for `left_vars` and Fact 3 for `right_vars`.
            fn split_at(mut self, index: usize) -> (Self, Self) {
                let left_offsets = &self.offsets[..=index];
                let right_offsets = &self.offsets[index..];

                let split_point = self.offsets[index];
                let (left_vars, right_vars) =
                    std::mem::take(&mut self.vars).split_at_mut(split_point);
                let left = Producer {
                    vars: left_vars,
                    offsets: left_offsets,
                };
                let right = Producer {
                    vars: right_vars,
                    offsets: right_offsets,
                };
                (left, right)
            }
        }

        // This construction of `Producer` ensures that the invariants are preserved
        // because `self.vars` and `self.offsets` are used directly from `LcMap`,
        // which enforces the relevant `LcVarsParIterMut` invariants.
        callback.callback(Producer {
            vars: self.vars,
            offsets: self.offsets,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_test_curves::bls12_381::Fr;

    #[test]
    #[cfg(feature = "parallel")]
    fn test_lc_map_par_iter_mut() {
        let mut interner = FieldInterner::<Fr>::new();

        let mut lcmap = LcMap::<Fr>::new();

        lcmap.push(
            [
                (1u8.into(), Variable::One),
                (2u8.into(), Variable::instance(2)),
            ],
            &mut interner,
        );
        lcmap.push(
            [
                (3u8.into(), Variable::witness(4)),
                (4u8.into(), Variable::instance(4)),
            ],
            &mut interner,
        );

        // Parallel mutation: double coefficients, increment vars
        lcmap.lc_vars_par_iter_mut().for_each(|chunk| {
            for v in chunk {
                if v.is_instance() {
                    *v = Variable::instance(v.index().unwrap() + 1);
                }
            }
        });

        // Convert back to (F, Variable) for assertions
        let flattened: Vec<_> = lcmap
            .iter()
            .flat_map(|chunk| to_non_interned_lc(chunk, &interner))
            .collect();

        let expected = vec![
            (1u8.into(), Variable::One),
            (2u8.into(), Variable::instance(3)),
            (3u8.into(), Variable::witness(4)),
            (4u8.into(), Variable::instance(5)),
        ];

        assert_eq!(flattened, expected);
    }

    #[test]
    fn test_lc_map_iter_mut() {
        let mut interner = FieldInterner::<Fr>::new();

        let mut lcmap = LcMap::<Fr>::new();

        lcmap.push(
            [
                (1u8.into(), Variable::One),
                (2u8.into(), Variable::instance(2)),
            ],
            &mut interner,
        );
        lcmap.push(
            [
                (3u8.into(), Variable::witness(4)),
                (4u8.into(), Variable::instance(4)),
            ],
            &mut interner,
        );

        // Parallel mutation: double coefficients, increment vars
        lcmap.lc_vars_iter_mut().for_each(|chunk| {
            for v in chunk {
                if v.is_instance() {
                    *v = Variable::instance(v.index().unwrap() + 1);
                }
            }
        });

        // Convert back to (F, Variable) for assertions
        let flattened: Vec<_> = lcmap
            .iter()
            .flat_map(|chunk| to_non_interned_lc(chunk, &interner))
            .collect();

        let expected = vec![
            (1u8.into(), Variable::One),
            (2u8.into(), Variable::instance(3)),
            (3u8.into(), Variable::witness(4)),
            (4u8.into(), Variable::instance(5)),
        ];

        assert_eq!(flattened, expected);
    }
}
