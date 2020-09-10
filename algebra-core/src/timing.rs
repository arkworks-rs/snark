// This instrumentation should only be used for functions
// which run on the order of >> 1ms, as the time for processing
// and printing is on the order of 20-70 us with whitelists and
// blacklists (due to needing to unwind the backtrace at runtime)
// and around 3 us without.

#[macro_export]
macro_rules! timer_println {
    ($now: ident, $string: expr) => {
        #[cfg(any(feature = "timing", feature = "timing_detailed"))]
        {
            use backtrace::Backtrace;
            #[cfg(feature = "timing_thread_id")]
            use thread_id;

            const MAX_CALL_DEPTH: usize = 10;

            let elapsed = $now.1.elapsed().as_micros();

            // This is for reference
            let _instrumented_functions: Vec<&'static str> = vec![
                "batch_bucketed_add",
                "verify_points",
                "batch_scalar_mul_in_place",
            ];

            let whitelisted_functions: Vec<&'static str> = vec!["verify_points"];

            let blacklisted_parent_functions: Vec<&'static str> = vec![];
            let whitelisted_parent_functions: Vec<&'static str> = vec![];

            macro_rules! function {
                () => {{
                    fn f() {}
                    fn type_name_of<T>(_: T) -> &'static str {
                        core::any::type_name::<T>()
                    }
                    let name = type_name_of(f);
                    &name[..name.len() - 3]
                }};
            }
            let func_string = function!();
            let mut fs_vec = func_string.split("::").collect::<Vec<_>>();
            while *fs_vec.last().unwrap() == "{{closure}}" {
                fs_vec.pop();
            }
            let func_name = *fs_vec.last().unwrap();
            let whitelisted = whitelisted_functions.iter().any(|&w| w == func_name);

            if whitelisted {
                let (blacklisted, whitelisted_parents) = if whitelisted_parent_functions.len() == 0
                    && blacklisted_parent_functions.len() == 0
                {
                    (false, true)
                } else {
                    let bt = Backtrace::new();
                    let mut bt_iter = bt.frames().iter().flat_map(|x| x.symbols());

                    let mut b = !(blacklisted_parent_functions.len() == 0);
                    let mut wp = whitelisted_parent_functions.len() == 0;

                    for _ in 0..MAX_CALL_DEPTH {
                        if b == true {
                            break;
                        }
                        if let Some(symbol) = bt_iter.next() {
                            let calling_func_string = format!("{}", symbol.name().unwrap());
                            let mut vec = calling_func_string.split("::").collect::<Vec<_>>();

                            vec.pop();
                            if let Some(func) = vec.last() {
                                if whitelisted_parent_functions.iter().any(|&x| x == *func) {
                                    wp = true;
                                }
                                if blacklisted_parent_functions.iter().any(|&x| x == *func) {
                                    b = true;
                                }
                            }
                        } else {
                            break;
                        }
                    }
                    (b, wp)
                };

                if !blacklisted && whitelisted_parents {
                    if cfg!(feature = "timing") {
                        let std_info = format!("[{:^28}] {} us", $string, elapsed);
                        #[cfg(feature = "timing_thread_id")]
                        let std_info =
                            format!("{:25} {}", format!("(tid: {})", thread_id::get()), std_info);
                        println!("{}", std_info);
                    }

                    if cfg!(feature = "timing_detailed") {
                        let std_info = format!(
                            "{:30} {:26} [{:^28}] {} us",
                            format!(
                                "{} {}:{}",
                                String::from(file!()).split("/").last().unwrap(),
                                $now.0,
                                line!()
                            ),
                            func_name,
                            $string,
                            elapsed
                        );
                        #[cfg(feature = "timing_thread_id")]
                        let std_info =
                            format!("{:25} {}", format!("(tid: {})", thread_id::get()), std_info);
                        println!("{}", std_info);
                    }
                }
            }
        }
    };
}

#[macro_export]
macro_rules! timer {
    () => {{
        #[cfg(any(feature = "timing", feature = "timing_detailed"))]
        let now = (line!(), std::time::Instant::now());

        #[cfg(not(any(feature = "timing", feature = "timing_detailed")))]
        let now = ();
        now
    }};
}
