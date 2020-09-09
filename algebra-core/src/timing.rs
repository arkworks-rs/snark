// This instrumentation should only be used for functions
// which run on the order of >> 1ms, as the time for processing
// and printing is on the order of 1-3 us
#[macro_export]
macro_rules! timer_println {
    ($now: ident, $string: expr) => {
        #[cfg(any(feature = "timing", feature = "timing_detailed"))]
        {
            #[cfg(feature = "timing_thread_id")]
            use thread_id;
            // This is for reference
            const _INSTRUMENTED_FUNCTIONS: [&'static str; 3] = [
                "batch_bucketed_add",
                "verify_points",
                "batch_scalar_mul_in_place",
            ];

            const WHITELISTED_FUNCTIONS: [&'static str; 1] = ["verify_points"];

            const BLACKLISTED_PARENT_FUNCTIONS: [&'static str; 0] = [];

            // If not empty, we only run the instrumentation if
            // one of the parents of the function is contained here
            const WHITELISTED_PARENT_FUNCTIONS: [&'static str; 0] = [];

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

            let whitelisted_parents = if WHITELISTED_PARENT_FUNCTIONS.len() == 0 {
                true
            } else {
                func_string
                    .split("::")
                    .any(|func| WHITELISTED_PARENT_FUNCTIONS.iter().any(|&x| x == func))
            };
            // Note this has n^2 complexity, please be cautious.
            let blacklisted = func_string
                .split("::")
                .any(|func| BLACKLISTED_PARENT_FUNCTIONS.iter().any(|&x| x == func));

            if !blacklisted && whitelisted_parents {
                let mut fs_vec = func_string.split("::").collect::<Vec<_>>();
                while *fs_vec.last().unwrap() == "{{closure}}" {
                    fs_vec.pop();
                }

                let func_name = *fs_vec.last().unwrap();
                let whitelisted = WHITELISTED_FUNCTIONS.iter().any(|&w| w == func_name);

                if cfg!(feature = "timing") && whitelisted {
                    let std_info = format!("[{:^28}] {} us", $string, $now.1.elapsed().as_micros());
                    #[cfg(feature = "timing_thread_id")]
                    let std_info =
                        format!("{:25} {}", format!("(tid: {})", thread_id::get()), std_info);
                    println!("{}", std_info);
                }

                if cfg!(feature = "timing_detailed") && whitelisted {
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
                        $now.1.elapsed().as_micros()
                    );
                    #[cfg(feature = "timing_thread_id")]
                    let std_info =
                        format!("{:25} {}", format!("(tid: {})", thread_id::get()), std_info);
                    println!("{}", std_info);
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
