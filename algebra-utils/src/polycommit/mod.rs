use algebra::AffineCurve;

#[cfg(not(feature = "gpu"))]
use rayon::prelude::*;
#[cfg(feature = "gpu")]
use algebra_kernels::polycommit::get_kernels;

pub fn polycommit_round_reduce<
    G: AffineCurve
>(
    round_challenge: G::ScalarField,
    round_challenge_inv: G::ScalarField,
    c_l: &mut [G::ScalarField],
    c_r: &[G::ScalarField],
    z_l: &mut [G::ScalarField],
    z_r: &[G::ScalarField],
    k_l: &mut [G::Projective],
    k_r: &[G],
) {
    #[cfg(not(feature = "gpu"))]
    {
        c_l.par_iter_mut()
            .zip(c_r)
            .for_each(|(c_l, c_r)| *c_l += &(round_challenge_inv * &c_r));

        z_l.par_iter_mut()
            .zip(z_r)
            .for_each(|(z_l, z_r)| *z_l += &(round_challenge * &z_r));

        k_l.par_iter_mut()
            .zip(k_r)
            .for_each(|(k_l, k_r)| *k_l += &(k_r.mul(round_challenge)));
    }

    #[cfg(feature = "gpu")]
    match get_kernels() {
        Ok(kernels) => {
            match kernels[0].polycommit_round_reduce(
                round_challenge,
                round_challenge_inv,
                c_l,
                c_r,
                z_l,
                z_r,
                k_l,
                k_r        
            ) {
                Ok(_) => {},
                Err(error) => { panic!("{}", error); }
            }
        },
        Err(error) => {
            panic!("{}", error);
        }
    }    
}

#[cfg(test)]
mod tests {

    use super::*;
    use algebra::curves::bls12_381::G1Projective;
    use algebra::fields::bls12_381::Fr;
    use algebra::{Field, UniformRand, ProjectiveCurve};

    #[test]
    fn test_polycommit_round_reduce() {
        
        use rayon::prelude::*;
        
        let mut rng = &mut rand::thread_rng();

        let round_challenge = Fr::rand(&mut rng);
        let round_challenge_inv = round_challenge.inverse().unwrap();

        let samples = 1 << 10;

        let mut coeffs_l = (0..samples)
            .map(|_| Fr::rand(&mut rng))
            .collect::<Vec<_>>();

        let coeffs_r = (0..samples)
            .map(|_| Fr::rand(&mut rng))
            .collect::<Vec<_>>();

        let mut z_l = (0..samples)
            .map(|_| Fr::rand(&mut rng))
            .collect::<Vec<_>>();

        let z_r= (0..samples)
            .map(|_| Fr::rand(&mut rng))
            .collect::<Vec<_>>();

        let mut key_proj_l= (0..samples)
            .map(|_| G1Projective::rand(&mut rng))
            .collect::<Vec<_>>();

        let key_r= (0..samples)
            .map(|_| G1Projective::rand(&mut rng).into_affine())
            .collect::<Vec<_>>();

        let mut gpu_coeffs_l = coeffs_l.clone();
        let gpu_coeffs_r = coeffs_r.clone();
        let mut gpu_z_l = z_l.clone();
        let gpu_z_r = z_r.clone();
        let mut gpu_key_proj_l = key_proj_l.clone();
        let gpu_key_r = key_r.clone();

        coeffs_l.par_iter_mut()
            .zip(coeffs_r)
            .for_each(|(c_l, c_r)| *c_l += &(round_challenge_inv * &c_r));

        z_l.par_iter_mut()
            .zip(z_r)
            .for_each(|(z_l, z_r)| *z_l += &(round_challenge * &z_r));

        key_proj_l.par_iter_mut()
            .zip(key_r)
            .for_each(|(k_l, k_r)| *k_l += &(k_r.mul(round_challenge)));

        polycommit_round_reduce(
            round_challenge,
            round_challenge_inv,
            &mut gpu_coeffs_l,
            &gpu_coeffs_r,
            &mut gpu_z_l,
            &gpu_z_r,
            &mut gpu_key_proj_l,
            &gpu_key_r
        );

        assert_eq!(coeffs_l, gpu_coeffs_l);
        assert_eq!(z_l, gpu_z_l);
        assert_eq!(key_proj_l, gpu_key_proj_l);
    }
}