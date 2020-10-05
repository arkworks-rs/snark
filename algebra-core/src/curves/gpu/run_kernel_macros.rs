

macro_rules! impl_run_kernel {
    () => {
        // We drop a lock only after the parallel portion has been handled
        fn par_run_kernel_sync<T>(
            ctx: &Context,
            bases_h: &[<G as ProjectiveCurve>::Affine],
            exps_h: &[BigInt],
            cuda_group_size: usize,
            lock: T,
        ) -> DeviceMemory<G> {
            assert_eq!(bases_h.len(), exps_h.len());
            let n = bases_h.len();

            let mut tables_h = vec![G::zero(); n * TABLE_SIZE];
            let mut exps_recode_h = vec![0u8; n * NUM_U8];

            let now = std::time::Instant::now();
            generate_tables_and_recoding(
                bases_h,
                &mut tables_h[..],
                exps_h,
                &mut exps_recode_h[..],
                true,
            );
            drop(lock);
            println!(
                "Generated tables and recoding: {}us",
                now.elapsed().as_micros()
            );

            let now = std::time::Instant::now();
            let mut out = DeviceMemory::<G>::zeros(&ctx, n);
            let mut tables = DeviceMemory::<G>::zeros(&ctx, n * TABLE_SIZE);
            let mut exps = DeviceMemory::<u8>::zeros(&ctx, n * NUM_U8);
            println!("Allocated device memory: {}us", now.elapsed().as_micros());

            let now = std::time::Instant::now();
            tables.copy_from_slice(&tables_h);
            exps.copy_from_slice(&exps_recode_h);
            println!("Copied data to device: {}us", now.elapsed().as_micros());

            let now = std::time::Instant::now();
            scalar_mul_kernel::scalar_mul(
                &ctx,
                n / cuda_group_size, // grid
                cuda_group_size,     // block
                (tables.as_ptr(), exps.as_ptr(), out.as_mut_ptr(), n as isize),
            )
            .expect("Kernel call failed");

            println!("Ran kernel: {}us", now.elapsed().as_micros());
            out
        }

        pub fn par_run_kernel(
            ctx: &Context,
            bases_h: &[<G as ProjectiveCurve>::Affine],
            exps_h: &[BigInt],
            cuda_group_size: usize,
        ) -> DeviceMemory<G> {
            assert_eq!(bases_h.len(), exps_h.len());
            let n = bases_h.len();

            let now = std::time::Instant::now();
            let mut tables = DeviceMemory::<G>::zeros(&ctx, n * TABLE_SIZE);
            let mut exps = DeviceMemory::<u8>::zeros(&ctx, n * NUM_U8);
            let mut out = DeviceMemory::<G>::zeros(&ctx, n);
            println!("Allocated device memory: {}us", now.elapsed().as_micros());

            let now = std::time::Instant::now();
            generate_tables_and_recoding(bases_h, &mut tables[..], exps_h, &mut exps[..], true);
            println!(
                "Generated tables and recoding: {}us",
                now.elapsed().as_micros()
            );
            // Accessible from CPU as usual Rust slice (though this will be slow)
            // Can this be changed to a memcpy?
            P::scalar_mul_kernel(
                &ctx,
                n / cuda_group_size, // grid
                cuda_group_size,     // block
                (tables.as_ptr(), exps.as_ptr(), out.as_mut_ptr(), n as isize),
            )
            .expect("Kernel call failed");
            out
        }
    };
}
