#![allow(unsafe_code)]

use ark_ff::Field;
use ark_std::vec::*;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::gr1cs::{
    field_interner::{FieldInterner, InternedField},
    Variable,
};

#[derive(Debug, Clone, Default)]
pub struct LcMap<F: Field> {
    vars: Vec<Variable>,
    coeffs: Vec<InternedField>,
    offsets: Vec<usize>, // Starting indices of each inner slice in `vars` and `coeffs``
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
type LcMapIterMutItem<'a> =
    core::iter::Zip<core::slice::IterMut<'a, InternedField>, core::slice::IterMut<'a, Variable>>;

impl<F: Field> LcMap<F> {
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            vars: Vec::new(),
            coeffs: Vec::new(),
            offsets: vec![0],
            _f: core::marker::PhantomData,
        }
    }

    #[inline(always)]
    pub fn with_capacity(capacity: usize) -> Self {
        let mut result = Self::new();
        result.vars.reserve(capacity * 2);
        result.coeffs.reserve(capacity * 2);
        result.offsets.reserve(capacity + 1);
        result
    }

    #[inline(always)]
    pub fn push(
        &mut self,
        v: impl IntoIterator<Item = (F, Variable)>,
        f_interner: &mut FieldInterner<F>,
    ) {
        for (coeff, var) in v {
            self.coeffs.push(f_interner.get_or_intern(coeff));
            self.vars.push(var);
        }
        self.offsets.push(self.coeffs.len());
    }

    #[inline(always)]
    pub fn push_by_ref<'a>(
        &mut self,
        v: impl IntoIterator<Item = &'a (F, Variable)>,
        f_interner: &mut FieldInterner<F>,
    ) {
        for (coeff, var) in v {
            self.coeffs.push(f_interner.get_or_intern(*coeff));
            self.vars.push(*var);
        }
        self.offsets.push(self.coeffs.len());
    }

    #[inline(always)]
    pub fn iter(&self) -> impl Iterator<Item = LcMapIterItem<'_>> {
        self.offsets
            .windows(2)
            .map(|w| windowed_access(w, &self.coeffs, &self.vars))
    }

    #[cfg(feature = "parallel")]
    #[inline(always)]
    pub fn par_iter(&self) -> impl ParallelIterator<Item = LcMapIterItem<'_>> {
        self.offsets
            .par_windows(2)
            .map(|w| windowed_access(w, &self.coeffs, &self.vars))
    }

    #[inline(always)]
    pub fn iter_mut(&mut self) -> impl Iterator<Item = LcMapIterMutItem<'_>> {
        LcMapIterMut {
            coeffs: &mut self.coeffs,
            vars: &mut self.vars,
            offsets: self.offsets.windows(2),
        }
    }

    #[cfg(feature = "parallel")]
    #[inline(always)]
    pub fn par_iter_mut(&mut self) -> LcMapParIterMut<'_> {
        LcMapParIterMut {
            coeffs: self.coeffs.as_mut_ptr(),
            vars: self.vars.as_mut_ptr(),
            offsets: &self.offsets,
        }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    #[allow(unsafe_code)]
    #[inline(always)]
    pub fn get(&self, idx: usize) -> Option<LcMapIterItem<'_>> {
        if idx >= self.len() {
            return cold();
        } else {
            Some(windowed_access(
                unsafe { self.offsets.get_unchecked(idx..=(idx + 1)) },
                &self.coeffs,
                &self.vars,
            ))
        }
    }
}

///
fn windowed_access<'a>(
    w: &'a [usize],
    coeffs: &'a [InternedField],
    vars: &'a [Variable],
) -> LcMapIterItem<'a> {
    unsafe {
        // SAFETY:
        // by construction, `self.offsets` always has an odd number of elements,
        // and so `w` always has two elements
        let start = *w.get_unchecked(0);
        let end = *w.get_unchecked(1);
        // `start` and `end` are guaranteed to be within bounds of `self.coeffs` and `self.vars`
        coeffs
            .get_unchecked(start..end)
            .iter()
            .zip(vars.get_unchecked(start..end))
    }
}

#[cold]
fn cold<'a>() -> Option<LcMapIterItem<'a>> {
    None
}

pub struct LcMapIterMut<'a> {
    coeffs: &'a mut [InternedField],
    vars: &'a mut [Variable],
    offsets: core::slice::Windows<'a, usize>,
}

impl<'a> Iterator for LcMapIterMut<'a> {
    type Item = LcMapIterMutItem<'a>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let Some([start, end]) = self.offsets.next() else {
            return None;
        };

        let len = end - start;
        // SAFETY:
        // * By construction, `start <= end` (we only append something larger than `start` to `offsets`)
        // * `since `end = start + length_of_lc`, `len` is guaranteed to be within bounds
        //    of `self.coeffs` and `self.vars`.
        //    This is because the latter two are always appended to together.
        #[allow(unsafe_code)]
        unsafe {
            let (c_head, c_tail) = core::mem::take(&mut self.coeffs).split_at_mut_unchecked(len);
            let (v_head, v_tail) = core::mem::take(&mut self.vars).split_at_mut_unchecked(len);
            self.coeffs = c_tail;
            self.vars = v_tail;
            Some(c_head.iter_mut().zip(v_head))
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.offsets.len() - 1;
        (len, Some(len))
    }
}

#[cfg(feature = "parallel")]
pub struct LcMapParIterMut<'a> {
    coeffs: *mut InternedField,
    vars: *mut Variable,
    offsets: &'a [usize],
}

#[cfg(feature = "parallel")]
#[allow(unsafe_code)]
unsafe impl<'a> Send for LcMapParIterMut<'a> {}

#[cfg(feature = "parallel")]
#[allow(unsafe_code)]
unsafe impl<'a> Sync for LcMapParIterMut<'a> {}

#[cfg(feature = "parallel")]
impl<'a> ParallelIterator for LcMapParIterMut<'a> {
    type Item = LcMapIterMutItem<'a>;

    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: rayon::iter::plumbing::UnindexedConsumer<Self::Item>,
    {
        rayon::iter::plumbing::bridge(self, consumer)
    }
}

#[cfg(feature = "parallel")]
impl<'a> IndexedParallelIterator for LcMapParIterMut<'a> {
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
        struct Producer<'a> {
            coeffs: *mut InternedField,
            vars: *mut Variable,
            offsets: &'a [usize],
        }

        #[allow(unsafe_code)]
        unsafe impl<'a> Send for Producer<'a> {}
        impl<'a> rayon::iter::plumbing::Producer for Producer<'a> {
            type Item = LcMapIterMutItem<'a>;
            type IntoIter = IntoIter<LcMapIterMutItem<'a>>;

            fn into_iter(self) -> Self::IntoIter {
                self.offsets
                    .windows(2)
                    .map(|window| {
                        // SAFETY: See logic from `LcMapIterMut::next`
                        let start = unsafe { *window.get_unchecked(0) };
                        let end = unsafe { *window.get_unchecked(1) };
                        let len = end - start;
                        #[allow(unsafe_code)]
                        unsafe {
                            let coeffs_ptr = self.coeffs.add(start);
                            let coeffs_slice = core::slice::from_raw_parts_mut(coeffs_ptr, len);
                            let vars_ptr = self.vars.add(start);
                            let vars_slice = core::slice::from_raw_parts_mut(vars_ptr, len);
                            coeffs_slice.iter_mut().zip(vars_slice)
                        }
                    })
                    .collect::<Vec<_>>()
                    .into_iter()
            }

            fn split_at(self, index: usize) -> (Self, Self) {
                let (left, right) = self.offsets.split_at(index + 1);
                (
                    Producer {
                        coeffs: self.coeffs,
                        vars: self.vars,
                        offsets: left,
                    },
                    Producer {
                        coeffs: self.coeffs,
                        vars: self.vars,
                        offsets: right,
                    },
                )
            }
        }

        callback.callback(Producer {
            coeffs: self.coeffs,
            vars: self.vars,
            offsets: self.offsets,
        })
    }
}
