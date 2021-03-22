use crate::{BasicRadix2Domain, EvaluationDomainImpl, MixedRadix2Domain};
use crate::{
    PrimeField, FpParameters,
};

/// Wraps EvaluationDomain variants, suited for
/// different field parameters and domain sizes.
#[derive(Copy, Clone, Hash, Eq, PartialEq, Debug)]
pub enum EvaluationDomain<F: PrimeField> {
    #[doc(hidden)]
    BasicRadix2Domain(BasicRadix2Domain<F>),
    #[doc(hidden)]
    MixedRadix2Domain(MixedRadix2Domain<F>),
}

impl<F: PrimeField> EvaluationDomainImpl<F> for EvaluationDomain<F> {

    fn new(num_coeffs: usize) -> Option<Self> {

        // Let's assign an index to each domain:
        // -1: No suitable domain found
        // 0: BasicRadix2Domain
        // 1: MixedRadix2Domain
        // 2: ...

        let mut index = -1;
        let mut domain_size = std::usize::MAX;

        match BasicRadix2Domain::<F>::compute_size_of_domain(num_coeffs) {
            Some(size) => {
                if size < domain_size {
                    index = 0;
                    domain_size = size;
                }
            },
            None => {}
        };

        if F::Params::SMALL_SUBGROUP_DEFINED {
            match MixedRadix2Domain::<F>::compute_size_of_domain(num_coeffs) {
                Some(size) => {
                    if size < domain_size {
                        index = 1;
                        //domain_size = size;
                    }
                },
                None => {}
            };
        }

        //Return best domain or None if no suitable domain has been found
        match index {
            0 => Some(EvaluationDomain::<F>::BasicRadix2Domain(BasicRadix2Domain::<F>::new(num_coeffs).unwrap())),
            1 => Some(EvaluationDomain::<F>::MixedRadix2Domain(MixedRadix2Domain::<F>::new(num_coeffs).unwrap())),
            _ => None,
        }
    }

    //Returns the smallest domain size among the domains able to handle at least num_coeffs
    fn compute_size_of_domain(num_coeffs: usize) -> Option<usize> {
        let mut compatible_domains_size: Vec<usize> =  vec![];

        let brd_size = BasicRadix2Domain::<F>::compute_size_of_domain(num_coeffs);
        if brd_size.is_some(){
            compatible_domains_size.push(brd_size.unwrap())
        }

        let mrd_size = MixedRadix2Domain::<F>::compute_size_of_domain(num_coeffs);
        if mrd_size.is_some(){
            compatible_domains_size.push(mrd_size.unwrap())
        }

        match compatible_domains_size.len() {
            0 => None,
            _ => {
                let mut curr_min = std::usize::MAX;
                compatible_domains_size.iter().for_each(|size| curr_min = curr_min.min(*size));
                Some(curr_min)
            }
        }
    }

    #[inline]
    fn size(&self) -> usize {
        match self {
            EvaluationDomain::BasicRadix2Domain(brd) => brd.size(),
            EvaluationDomain::MixedRadix2Domain(mrd) => mrd.size(),

        }
    }

    #[inline]
    fn size_inv(&self) -> F {
        match self {
            EvaluationDomain::BasicRadix2Domain(brd) => brd.size_inv(),
            EvaluationDomain::MixedRadix2Domain(mrd) => mrd.size_inv(),
        }
    }

    #[inline]
    fn group_gen(&self) -> F {
        match self {
            EvaluationDomain::BasicRadix2Domain(brd) => brd.group_gen(),
            EvaluationDomain::MixedRadix2Domain(mrd) => mrd.group_gen(),
        }
    }

    fn fft_in_place(&self, coeffs: &mut Vec<F>) {
        match self {
            EvaluationDomain::BasicRadix2Domain(brd) => brd.fft_in_place(coeffs),
            EvaluationDomain::MixedRadix2Domain(mrd) => mrd.fft_in_place(coeffs),
        }
    }

    fn coset_fft_in_place(&self, coeffs: &mut Vec<F>){
        match self {
            EvaluationDomain::BasicRadix2Domain(brd) => brd.coset_fft_in_place(coeffs),
            EvaluationDomain::MixedRadix2Domain(mrd) => mrd.coset_fft_in_place(coeffs),
        }
    }

    #[inline]
    fn ifft_in_place(&self, evals: &mut Vec<F>){
        match self {
            EvaluationDomain::BasicRadix2Domain(brd) => brd.ifft_in_place(evals),
            EvaluationDomain::MixedRadix2Domain(mrd) => mrd.ifft_in_place(evals),
        }
    }

    fn coset_ifft_in_place(&self, evals: &mut Vec<F>){
        match self {
            EvaluationDomain::BasicRadix2Domain(brd) => brd.coset_ifft_in_place(evals),
            EvaluationDomain::MixedRadix2Domain(mrd) => mrd.coset_ifft_in_place(evals)
        }
    }
}


impl<F: PrimeField> Default for EvaluationDomain<F> {
    fn default() -> Self {
        EvaluationDomain::BasicRadix2Domain(BasicRadix2Domain::<F>::default())
    }
}

#[cfg(test)]
mod test {
    use crate::fields::mnt6753::fr::Fr;
    use crate::domain::*;

    #[test]
    fn test_mnt6753_best_evaluation_domain() {
        //The basic domain size increases with 2^k, with k <= 15, while
        //the mixed domain increases with 2^k * 5^s, with k <= 15 and s <= 2

        let mut domain_size = 2048;
        let mut domain = EvaluationDomain::<Fr>::new(domain_size).unwrap();

        match domain {
            EvaluationDomain::BasicRadix2Domain(_) => assert_eq!(domain.size(), 2048, "Unexpected domain size"),
            _ => panic!("Unexpected chosen domain"),
        };

        domain_size = 5000;
        domain = EvaluationDomain::<Fr>::new(domain_size).unwrap();
        match domain {
            EvaluationDomain::MixedRadix2Domain(_) => assert_eq!(domain.size(), 5120, "Unexpected domain size"),
            _ => panic!("Unexpected chosen domain"),
        };

        //Limit for the basic radix2 domain support
        domain_size = 32768;
        domain = EvaluationDomain::<Fr>::new(domain_size).unwrap();
        match domain {
            EvaluationDomain::BasicRadix2Domain(_) => assert_eq!(domain.size(), 32768, "Unexpected domain size"),
            _ => panic!("Unexpected chosen domain"),
        };

        domain_size = 32769;
        domain = EvaluationDomain::<Fr>::new(domain_size).unwrap();
        match domain {
            EvaluationDomain::MixedRadix2Domain(_) => assert_eq!(domain.size(), 40960, "Unexpected domain size"),
            _ => panic!("Unexpected chosen domain"),
        };

        //Limit for the mixed radix2 domain support
        domain_size = 819200;
        domain = EvaluationDomain::<Fr>::new(domain_size).unwrap();
        match domain {
            EvaluationDomain::MixedRadix2Domain(_) => assert_eq!(domain.size(), 819200, "Unexpected domain size"),
            _ => panic!("Unexpected chosen domain"),
        };

        //No supported domain for this size should exist
        domain_size = 819201;
        match EvaluationDomain::<Fr>::new(domain_size) {
            None => {},
           _ => panic!("No domain should exists for this size")
        }
    }
}