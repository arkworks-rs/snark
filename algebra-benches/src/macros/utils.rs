const N: usize = 1000;

macro_rules! tmp_dot_func {
    ($tmp:ident, $v:ident, $func:ident, $count:ident) => {
        #[cfg(not(feature = "n_fold"))]
        $tmp.$func(&$v[$count].1);
        #[cfg(feature = "n_fold")]
        for _ in 0..N { $tmp.$func(&$v[$count].1); }
    }
}

macro_rules! prepared_v {
    ($v:ident, $rng:ident) => {
        let $v: Vec<(G1Prepared<Parameters>, G2Prepared<Parameters>)> = (0..SAMPLES)
            .map(|_| {
                (
                    G1Affine::from(G1::rand(&mut $rng)).into(),
                    G2Affine::from(G2::rand(&mut $rng)).into(),
                )
            })
            .collect();
    }
}

macro_rules! affine_v {
    ($v:ident, $rng:ident) => {
        let $v: Vec<(G1Affine, G2Affine)> = (0..SAMPLES)
            .map(|_| {
                (
                    G1Affine::from(G1::rand(&mut $rng)).into(),
                    G2Affine::from(G2::rand(&mut $rng)).into(),
                )
            })
            .collect();
    }
}
