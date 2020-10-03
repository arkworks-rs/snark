// TODO: make this more generic
#[macro_export]
macro_rules! impl_gpu_cpu_run_kernel {
    ($STATIC_MICROBENCH: ident) =>  {
        use peekmore::PeekMore;
        use closure::closure;

        pub fn cpu_gpu_static_partition_run_kernel(
            ctx: &Context,
            bases_h: &[<G as ProjectiveCurve>::Affine],
            exps_h: &[BigInt],
            cuda_group_size: usize,
            // size of the batch for cpu scalar mul
            cpu_chunk_size: usize,
        ) -> Vec<<G as ProjectiveCurve>::Affine> {
            let mut bases_res = bases_h.to_vec();

            let now = std::time::Instant::now();
            let mut profile_data = $STATIC_MICROBENCH.lock().unwrap();
            let ratio = profile_data.0;

            let n = bases_res.len();
            let n_cpu = (ratio * (n as f64)).round() as usize;
            let n_gpu = n - n_cpu;

            let (bases_cpu, bases_gpu) = bases_res.split_at_mut(n_cpu);
            let (exps_cpu, exps_gpu) = exps_h.split_at(n_cpu);

            let mut tables = DeviceMemory::<G>::zeros(&ctx, n_gpu * TABLE_SIZE);
            let mut exps = DeviceMemory::<u8>::zeros(&ctx, n_gpu * NUM_U8);

            let (mut time_cpu, mut time_gpu) = (0, 0);

            println!("Split statically and allocated device: {}us", now.elapsed().as_micros());

            generate_tables_and_recoding(bases_gpu, &mut tables[..], exps_gpu, &mut exps[..], true);

            rayon::scope(|s| {
                // Here, we should write directly to device
                s.spawn(|_| {
                    let now = std::time::Instant::now();
                    let mut out = DeviceMemory::<G>::zeros(&ctx, n_gpu);
                    scalar_mul_kernel::scalar_mul(
                        &ctx,
                        (n_gpu - 1) / cuda_group_size + 1, // grid
                        cuda_group_size,     // block
                        (tables.as_ptr(), exps.as_ptr(), out.as_mut_ptr(), n_gpu as isize),
                    )
                    .expect("Kernel call failed");
                    G::batch_normalization(&mut out[..]);
                    bases_gpu.clone_from_slice(&out.par_iter().map(|p| p.into_affine()).collect::<Vec<_>>()[..]);
                    time_gpu = now.elapsed().as_micros();
                    println!("GPU finish");
                });

                s.spawn(|_| {
                    let now = std::time::Instant::now();
                    let exps_mut = &mut exps_cpu.to_vec()[..];
                    rayon::scope(|t| {
                        for (b, s) in bases_cpu.chunks_mut(cpu_chunk_size).zip(exps_mut.chunks_mut(cpu_chunk_size)) {
                            t.spawn(move |_| b[..].batch_scalar_mul_in_place(&mut s[..], 4));
                        }
                    });
                    time_cpu = now.elapsed().as_micros();
                    println!("CPU finish");
                });
            });

            // Update global microbenchmarking state
            println!("old profile_data: {:?}", profile_data);
            let cpu_throughput = n_cpu as f64 / time_cpu as f64;
            let gpu_throughput = n_gpu as f64 / time_gpu as f64;
            let new_ratio = cpu_throughput / (cpu_throughput + gpu_throughput);
            println!("new ratio: {:?}", new_ratio);
            let n_data_points = profile_data.1 as f64;
            profile_data.1 += 1;
            profile_data.0 = (new_ratio + n_data_points * profile_data.0) / profile_data.1 as f64;
            println!("new profile_data: {:?}", profile_data);

            bases_res
        }

        pub fn cpu_gpu_load_balance_run_kernel(
            ctx: &Context,
            bases_h: &[<G as ProjectiveCurve>::Affine],
            exps_h: &[BigInt],
            cuda_group_size: usize,
            // size of a single job in the queue e.g. 2 << 14
            job_size: usize,
            // size of the batch for cpu scalar mul
            cpu_chunk_size: usize,
        ) -> Vec<<G as ProjectiveCurve>::Affine> {
            let mut bases_res = bases_h.to_vec();
            let queue = Mutex::new(bases_res.chunks_mut(job_size).zip(exps_h.chunks(job_size)).peekmore());

            rayon::scope(|s| {
                // We launch two concurrent GPU threads that block on waiting for GPU to hide latency
                for i in 0..2 {
                    s.spawn(closure!(move i, ref queue, |_| {
                        std::thread::sleep(std::time::Duration::from_millis(i * 500));
                        let mut iter = queue.lock().unwrap();
                        while let Some((bases, exps)) = iter.next() {
                            iter.peek();
                            if iter.peek().is_none() { break; }
                            let mut proj_res = par_run_kernel_sync(ctx, bases, exps, cuda_group_size, iter);
                            G::batch_normalization(&mut proj_res[..]);
                            bases.clone_from_slice(&proj_res.par_iter().map(|p| p.into_affine()).collect::<Vec<_>>()[..]);
                            iter = queue.lock().unwrap();
                        }
                    }));
                }

                s.spawn(|_| {
                    std::thread::sleep(std::time::Duration::from_millis(20));
                    let mut iter = queue.lock().unwrap();
                    println!("acquired cpu");
                    while let Some((bases, exps)) = iter.next() {
                        let exps_mut = &mut exps.to_vec()[..];
                        rayon::scope(|t| {
                            for (b, s) in bases.chunks_mut(cpu_chunk_size).zip(exps_mut.chunks_mut(cpu_chunk_size)) {
                                t.spawn(move |_| b[..].batch_scalar_mul_in_place(&mut s[..], 4));
                            }
                        });
                        // Sleep to allow other threads to unlock
                        drop(iter);
                        println!("unlocked cpu");
                        std::thread::sleep(std::time::Duration::from_millis(20));
                        iter = queue.lock().unwrap();
                        println!("acquired cpu");
                    }
                    println!("CPU FINISH");
                });
            });
            drop(queue);
            bases_res
        }
    }
}
