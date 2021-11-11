macro_rules! n_fold {
    ($tmp:ident, $v:ident, $func:ident, $count:ident) => {
        #[cfg(not(feature = "n_fold"))]
        $tmp.$func(&$v[$count].1);
        #[cfg(feature = "n_fold")]
        for _ in 0..1000 {
            $tmp.$func(&$v[$count].1);
        }
    };

    ($tmp:ident, $func:ident) => {
        #[cfg(not(feature = "n_fold"))]
        $tmp.$func();
        #[cfg(feature = "n_fold")]
        for _ in 0..1000 {
            $tmp.$func();
        }
    };
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
    };
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
    };
}
