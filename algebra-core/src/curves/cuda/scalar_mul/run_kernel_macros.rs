#[macro_export]
macro_rules! impl_run_kernel {
    () => {
        // We drop a lock only after the parallel portion has been handled
        fn par_run_kernel_sync<T>(
            ctx: &Context,
            bases_h: &[<Self as ProjectiveCurve>::Affine],
            exps_h: &[<<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt],
            cuda_group_size: usize,
            lock: T,
        ) -> DeviceMemory<Self> {
            assert_eq!(bases_h.len(), exps_h.len());
            let n = bases_h.len();

            let mut tables_h = vec![Self::zero(); n * Self::table_size()];
            let mut exps_recode_h = vec![0u8; n * Self::num_u8()];

            let now = std::time::Instant::now();
            Self::generate_tables_and_recoding(
                bases_h,
                &mut tables_h[..],
                exps_h,
                &mut exps_recode_h[..],
            );
            drop(lock);
            println!(
                "Generated tables and recoding: {}us",
                now.elapsed().as_micros()
            );

            let now = std::time::Instant::now();
            let mut out = DeviceMemory::<Self>::zeros(&ctx, n);
            let mut tables = DeviceMemory::<Self>::zeros(&ctx, n * Self::table_size());
            let mut exps = DeviceMemory::<u8>::zeros(&ctx, n * Self::num_u8());
            println!("Allocated device memory: {}us", now.elapsed().as_micros());

            let now = std::time::Instant::now();
            tables.copy_from_slice(&tables_h);
            exps.copy_from_slice(&exps_recode_h);
            println!("Copied data to device: {}us", now.elapsed().as_micros());

            let now = std::time::Instant::now();
            P::scalar_mul_kernel(
                &ctx,
                n / cuda_group_size, // grid
                cuda_group_size,     // block
                tables.as_ptr(),
                exps.as_ptr(),
                out.as_mut_ptr(),
                n as isize,
            )
            .expect("Kernel call failed");

            println!("Ran kernel: {}us", now.elapsed().as_micros());
            out
        }

        fn par_run_kernel(
            ctx: &Context,
            bases_h: &[<Self as ProjectiveCurve>::Affine],
            exps_h: &[<<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt],
            cuda_group_size: usize,
        ) -> DeviceMemory<Self> {
            assert_eq!(bases_h.len(), exps_h.len());
            let n = bases_h.len();

            let now = std::time::Instant::now();
            let mut tables = DeviceMemory::<Self>::zeros(&ctx, n * Self::table_size());
            let mut exps = DeviceMemory::<u8>::zeros(&ctx, n * Self::num_u8());
            let mut out = DeviceMemory::<Self>::zeros(&ctx, n);
            println!("Allocated device memory: {}us", now.elapsed().as_micros());

            let now = std::time::Instant::now();
            Self::generate_tables_and_recoding(bases_h, &mut tables[..], exps_h, &mut exps[..]);
            println!(
                "Generated tables and recoding: {}us",
                now.elapsed().as_micros()
            );
            P::scalar_mul_kernel(
                &ctx,
                n / cuda_group_size, // grid
                cuda_group_size,     // block
                tables.as_ptr(),
                exps.as_ptr(),
                out.as_mut_ptr(),
                n as isize,
            )
            .expect("Kernel call failed");
            out
        }
    };
}
