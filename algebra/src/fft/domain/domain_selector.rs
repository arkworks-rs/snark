use crate::{BasicRadix2Domain, EvaluationDomain, MixedRadix2Domain};
use crate::{FpParameters, PrimeField};

/// Return the smallest sized and most efficient Evaluation Domain able to support `num_coeffs` size
pub fn get_best_evaluation_domain<F: PrimeField>(
    num_coeffs: usize,
) -> Option<Box<dyn EvaluationDomain<F>>> {
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
        }
        None => {}
    };

    if F::Params::SMALL_SUBGROUP_DEFINED {
        match MixedRadix2Domain::<F>::compute_size_of_domain(num_coeffs) {
            Some(size) => {
                if size < domain_size {
                    index = 1;
                    //domain_size = size;
                }
            }
            None => {}
        };
    }

    //Return best domain or None if no suitable domain has been found
    match index {
        0 => Some(Box::new(BasicRadix2Domain::<F>::new(num_coeffs).unwrap())),
        1 => Some(Box::new(MixedRadix2Domain::<F>::new(num_coeffs).unwrap())),
        _ => None,
    }
}

#[cfg(test)]
mod test {
    use crate::domain::*;
    use crate::fields::mnt6753::fr::Fr;

    #[test]
    fn test_mnt6753_best_evaluation_domain() {
        //The basic domain size increases with 2^k, with k <= 15, while
        //the mixed domain increases with 2^k * 5^s, with k <= 15 and s <= 2

        let mut domain_size = 2048;
        let mut domain = get_best_evaluation_domain::<Fr>(domain_size).unwrap();

        //Expected Basic to be chosen
        assert_eq!(domain.size(), 2048, "Unexpected domain size");

        domain_size = 5000;
        domain = get_best_evaluation_domain::<Fr>(domain_size).unwrap();
        //Expected Mixed to be chosen
        assert_eq!(domain.size(), 5120, "Unexpected domain size");

        //Limit for the basic radix2 domain support
        domain_size = 32768;
        domain = get_best_evaluation_domain::<Fr>(domain_size).unwrap();
        assert_eq!(domain.size(), 32768, "Unexpected domain size");

        domain_size = 32769;
        //Expected Mixed to be chosen
        domain = get_best_evaluation_domain::<Fr>(domain_size).unwrap();
        assert_eq!(domain.size(), 40960, "Unexpected domain size");

        //Limit for the mixed radix2 domain support
        domain_size = 819200;
        domain = get_best_evaluation_domain::<Fr>(domain_size).unwrap();
        assert_eq!(domain.size(), 819200, "Unexpected domain size");

        //No supported domain for this size should exist
        domain_size = 819201;
        match get_best_evaluation_domain::<Fr>(domain_size) {
            None => {}
            _ => panic!("No domain should exists for this size"),
        }
    }
}
